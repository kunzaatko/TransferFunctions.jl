using ImageFiltering: padarray, Fill
using SpecialFunctions, FourierTools # , IterTools

# TODO: Add examples with transfer functions. For example an model OTF, when `ifft`ed with padding of the same parity
# should be real...
@doc """
    padtosize(a::AbstractArray{T,N}, size...; fourier=false, padvalue=0)
    padtosize(a::AbstractArray{T,N}, size; fourier=false, padvalue=0)

Pad `a` to a `size`

Padding is done in a way which ensures that if `a` is symmetric and decreases to zero on the edges (!!) such that
`eltype(fft(a)) <: Real`, then, if possible, output is symmetric under the DFT and `elyptype(fft(a_padded)) <: Real`).
If `a` is in the Fourier domain, then `fourier` should de set to `true`.

# Arguments:
- `size`: size of the output array. If `size isa Integer` then size along all the dimensions is `size`
- `fourier`: `a` is in the Fourier domain (default: `false`).

# Examples:
```jldoctest; setup = :(using TransferFunctions: padtosize; using FFTW)
julia> padtosize(reshape(1:4, 2,2), 3, 3)
3×3 Matrix{Int64}:
 0  0  0
 0  1  3
 0  2  4

julia> padtosize(reshape(1:4, 2,2), 4, 4)
4×4 Matrix{Int64}:
 0  0  0  0
 0  1  3  0
 0  2  4  0
 0  0  0  0

julia> padtosize(reshape(1:4, 2,2), 3, 3; fourier=true)
3×3 Matrix{Int64}:
 4  2  0
 3  1  0
 0  0  0

julia> padtosize(reshape(1:4, 2,2), 4, 4; fourier=true)
4×4 Matrix{Int64}:
 1  0  0  3
 0  0  0  0
 0  0  0  0
 2  0  0  4

julia> ones(1,1) |> fft
1×1 Matrix{ComplexF64}:
 1.0 + 0.0im

julia> padtosize(ones(1,1), 2,2) |> fft
2×2 Matrix{ComplexF64}:
  1.0+0.0im  -1.0+0.0im
 -1.0+0.0im   1.0+0.0im

julia> ones(2,2) |> fft
2×2 Matrix{ComplexF64}:
 4.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im

julia> padtosize(ones(2,2), 3, 3) |> fft # can be made symmetric
3×3 Matrix{ComplexF64}:
  4.0+0.0im  -2.0+0.0im  -2.0+0.0im
 -2.0+0.0im   1.0+0.0im   1.0+0.0im
 -2.0+0.0im   1.0+0.0im   1.0+0.0im

julia> padtosize(ones(2,2), 4, 4) |> fft # cannot be symmetric
4×4 Matrix{ComplexF64}:
  4.0+0.0im  -2.0-2.0im  0.0+0.0im  -2.0+2.0im
 -2.0-2.0im   0.0+2.0im  0.0+0.0im   2.0+0.0im
  0.0+0.0im   0.0+0.0im  0.0+0.0im   0.0+0.0im
 -2.0+2.0im   2.0+0.0im  0.0+0.0im   0.0-2.0im

julia> @assert all(-10^-5 .< imag(fft(padtosize(ones(2,2), 5, 5))) .< 10^-5) # .== 0 (numerical errors)

```
"""
function padtosize(
    a::AbstractArray{T,N}, s::Vararg{Int,N}; fourier=false, padvalue=zero(T)
) where {N,T}
    padxyzpost = Tuple(floor.(Int, (s .- size(a)) ./ 2))
    padxyzpre = padxyzpost .+ isodd.(size(a) .- s)

    # FIX: Is `fourier` correct?! <10-12-23> 
    # PERF: Should be using fftshift and ifftshift views  
    a = fourier ? fftshift(a) : a
    padded = padarray(a, Fill(padvalue, padxyzpre, padxyzpost)).parent
    return fourier ? ifftshift(padded) : padded
end
padtosize(a, size::Integer; vargs...) = padtosize(a, fill(size, ndims(a))...; vargs...)

@doc """
    roundupcenter(a::AbstractArray)
    roundupcenter((s...))

 Calculate center index of `a` rounded-up (i.e. `fft` center)

# Examples
```jldoctest; setup = :(using TransferFunctions: roundupcenter)
julia> roundupcenter(ones(15,15))
(8, 8)

julia> roundupcenter(ones(16,16))
(9, 9)

julia> roundupcenter(ones(15,15,15))
(8, 8, 8)
```
 """
# TODO: Document the tuple method <10-12-23> 
roundupcenter(s::NTuple{N,<:Int}) where {N} = round.(Int, s ./ 2, RoundUp) .+ iseven.(s)
roundupcenter(a::AbstractArray) = roundupcenter(size(a))

@doc """
     exactcenter(a::AbstractArray)

 Calculate the exact center index of `a`

# Examples
```jldoctest; setup = :(using TransferFunctions: exactcenter)
julia> exactcenter(ones(15,15))
(8.0, 8.0)

julia> exactcenter(ones(16,16))
(8.5, 8.5)

julia> exactcenter(ones(15,15,15))
(8.0, 8.0, 8.0)
```
 """
exactcenter(a::AbstractArray) = (size(a) .+ 1) ./ 2

# FIX: Type stability <13-12-23> 
# FIX: Domain error, when ∉ (-1,1) <13-12-23> 
# https://mathworld.wolfram.com/ApodizationFunction.html
abstract type Apodization end
Broadcast.broadcastable(a::Apodization) = Ref(a)
apodization(apo::Apodization, x::Real, halfwidth::Int) = apodization(apo, x / halfwidth)

@doc raw"""
 ``A(x) = 1-|x|/a``

 ``I(k) = a sinc²(π k a)``
"""
struct Bartlett <: Apodization
    a::Real
end
apodization(apo::Bartlett, x::Real) = 1 - abs(x) / apo.a
instrument(apo::Bartlett, k::Real) = apo.a * sinc(π * k * apo.a)^2

struct Blackman <: Apodization
    a::Real
end
function apodization(apo::Blackman, x::Real)
    return 21 / 50 + cospi(x / apo.a) + 2 * cospi(2 * x / apo.a) / 25
end
function instrument(apo::Blackman, k::Real)
    return apo.a * sinc(2π * k * apo.a) * (21 / 25 - 9apo.a^2 * k^2 / 25) /
           ((1 - apo.a^2 * k^2) * (1 - 4apo.a^2 * k^2))
end

struct Connes <: Apodization
    a::Real
end
apodization(apo::Connes, x::Real) = (1 - x^2 / apo.a^2)^2
function instrument(apo::Connes, k::Real)
    return 8apo.a * sqrt(2π) * SpecialFunctions.besselj(5 / 2, 2π * k * apo.a) /
           (2π * k * apo.a)^(5 / 2)
end

struct Cosine <: Apodization
    a::Real
end
apodization(apo::Cosine, x::Real) = cospi(x / 2apo.a)
instrument(apo::Cosine, k::Real) = 4apo.a * cospi(2apo.a * k) / (π * (1 - 16apo.a^2 * k^2))

struct Gaussian <: Apodization
    σ::Real
end
apodization(apo::Gaussian, x::Real) = exp(-x^2 / (2apo.σ^2))
function instrument(_::Gaussian, _::Real)
    throw(NotImplementedError("Low Priority (if you need it, please open an issue)"))
end

struct Hamming <: Apodization
    a::Real
end
apodization(apo::Hamming, x::Real) = 27 / 50 + 23cospi(x / apo.a) / 50
function instrument(apo::Hamming, k::Real)
    return apo.a(27 / 25 - 16apo.a^2 * k^2 / 25) * sinc(2π * apo.a * k) /
           (1 - 4apo.a^2 * k^2)
end

struct Hanning <: Apodization
    a::Real
end
apodization(apo::Hanning, x::Real) = cospi(x / 2apo.a)^2
instrument(apo::Hanning, k::Real) = apo.a * sinc(2π * apo.a * k) / (1 - 4apo.a^2 * k^2)

struct Welch <: Apodization
    a::Real
end
apodization(apo::Welch, x::Real) = 1 - x^2 / apo.a^2
function instrument(apo::Welch, k::Real)
    return apo.a * (sinpi(2 * k * apo.a) - 2π * apo.a * k * cospi(2apo.a * k)) /
           (2(apo.a * k * π)^3)
end

"""
    taperedges(apo, A, width::Int, [dims])
    taperedges(apo, A, (w1,...,wM), [(d1,...,dM)])
    taperedges(apo, A, ((w1_start,...,wM_start), (w1_end,...,wM_end)), [(d1,...,dM)])

Taper the edges of the input array `A` using the apodization function `apo::Apodization`, along `dims::Tuple`

- `apo::Apodization`: apodization function to apply.
- `A::AbstractArray{<:Number,N}`: The input array to taper.
- `width::NTuple{M,Int}` or `Int` or `Tuple{NTuple{M,Int},NTuple{M,Int}}`: The width of tapering at the edges.
- `dims::NTuple{M,Int}`: The dimensions along which to taper the edges.
"""
function taperedges(
    apo::Apodization,
    A::AbstractArray{<:Number,N},
    width::NTuple{M,Int},
    dims::NTuple{M,Int}=Tuple(1:N),
) where {N,M}
    return taperedges(apo, A, (width, width), dims)
end
function taperedges(
    apo::Apodization,
    A::AbstractArray{<:Number,N},
    width::Int,
    dims::NTuple{M,Int}=Tuple(1:N),
) where {N,M}
    return taperedges(apo, A, (Tuple(fill(width, M)), Tuple(fill(width, M))), dims)
end
function taperedges(
    apo::Apodization,
    A::AbstractArray{<:Number,N},
    width::Tuple{NTuple{M,Int},NTuple{M,Int}},
    dims::NTuple{M,Int}=Tuple(1:N),
) where {N,M}
    for (d, low, high) in zip(dims, width[1], width[2])
        lowedge = range(-1, 0; length=low + 1)
        for (e, i) in enumerate(firstindex(axes(A, d)):(firstindex(axes(A, d))+low-1))
            selectdim(A, d, i) .*= apodization(apo, lowedge[e])
        end
        highedge = range(0, 1; length=high + 1)
        for (e, i) in enumerate((lastindex(axes(A, d))-high+1):lastindex(axes(A, d)))
            selectdim(A, d, i) .*= apodization(apo, highedge[e])
        end
    end
    return A
end

# TODO: Add documentation <21-12-23> 
function apodize(
    apo::Apodization,
    A::AbstractArray{<:Number,N},
    cutoff::Real, # Tuple{NTuple{M,Int},NTuple{M,Int}}, # FIX: Generalize to non-symmetric cut-offs <21-12-23> 
    width::Real, # FIX: Generalize to non-symmetric widths <21-12-23> 
    # dims::NTuple{M,Int}=Tuple(1:N) # FIX: Abstract over dimensions <22-12-23> 
) where {N} # ,M}
    # PERF: Should be done with no allocation... This is the KISS solution  
    # FIX: Work for all dimensions  
    rs = [hypot(abs(x), abs(y)) for x in fftfreq(size(A, 1), size(A, 1)), y in fftfreq(size(A, 2), size(A, 2))]
    coefs = ones(eltype(A), size(rs)...)
    coefs[rs.>=cutoff] .= 0
    coefs[rs.<=(cutoff-width)] .= 1
    coefs[cutoff.>rs.>(cutoff-width)] .= map(r -> apodization(apo, r - (cutoff - width), width), rs[cutoff.>rs.>(cutoff-width)])
    A .* coefs
end
