using TransferFunctions: Length, PixelSize, Coordinate
using ImageFiltering: padarray
@reexport using ImageFiltering: Inner, Pad, Fill
using SpecialFunctions

fillsize(Δ::Length, n::Union{Int,Val})::PixelSize = ntuple(_ -> Δ, n)

# NOTE: Currently it is not used anywhere else anyway, so it is quite cheap to just throw it away. <27-08-24> 
# TODO: This should be probably left to the `ImageFiltering.jl` package. It already has a function that pads the
# needed amount <26-08-24> 
# TODO: Add examples with transfer functions. For example an model OTF, when `ifft`ed with padding of the same parity
# should be real...
"""
    TransferFunctions.padtosize(A, size...; fourier=false, padvalue=0)
    TransferFunctions.padtosize(A, size; fourier=false, padvalue=0)

Pad array `A` to a `size` with the padding value `padvalue`.

Padding is done in a way which ensures that if `a` is symmetric and decreases to zero on the edges (!!) such that
`eltype(fft(a)) <: Real`, then, if possible, output is symmetric under the DFT and `elyptype(fft(a_padded)) <: Real`).
If `a` is in the Fourier domain, then `fourier` should de set to `true`.

# Arguments:
- `size`: size of the output array. If `size isa Integer` then size along all the dimensions is `size`
- `fourier`: `a` is in the Fourier domain (default: `false`).

# Examples:
```jldoctest; setup = :(using TransferFunctions: TransferFunctions as TF; using FFTW)
julia> TF.padtosize(reshape(1:4, 2,2), 3, 3)
3×3 Matrix{Int64}:
 0  0  0
 0  1  3
 0  2  4

julia> TF.padtosize(reshape(1:4, 2,2), 4, 4)
4×4 Matrix{Int64}:
 0  0  0  0
 0  1  3  0
 0  2  4  0
 0  0  0  0

julia> TF.padtosize(reshape(1:4, 2,2), 3, 3; fourier=true)
3×3 Matrix{Int64}:
 4  2  0
 3  1  0
 0  0  0

julia> TF.padtosize(reshape(1:4, 2,2), 4, 4; fourier=true)
4×4 Matrix{Int64}:
 1  0  0  3
 0  0  0  0
 0  0  0  0
 2  0  0  4

julia> ones(1,1) |> fft
1×1 Matrix{ComplexF64}:
 1.0 + 0.0im

julia> TF.padtosize(ones(1,1), 2,2) |> fft
2×2 Matrix{ComplexF64}:
  1.0+0.0im  -1.0+0.0im
 -1.0+0.0im   1.0+0.0im

julia> ones(2,2) |> fft
2×2 Matrix{ComplexF64}:
 4.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im

julia> TF.padtosize(ones(2,2), 3, 3) |> fft # can be made symmetric
3×3 Matrix{ComplexF64}:
  4.0+0.0im  -2.0+0.0im  -2.0+0.0im
 -2.0+0.0im   1.0+0.0im   1.0+0.0im
 -2.0+0.0im   1.0+0.0im   1.0+0.0im

julia> TF.padtosize(ones(2,2), 4, 4) |> fft # cannot be symmetric
4×4 Matrix{ComplexF64}:
  4.0+0.0im  -2.0-2.0im  0.0+0.0im  -2.0+2.0im
 -2.0-2.0im   0.0+2.0im  0.0+0.0im   2.0+0.0im
  0.0+0.0im   0.0+0.0im  0.0+0.0im   0.0+0.0im
 -2.0+2.0im   2.0+0.0im  0.0+0.0im   0.0-2.0im

julia> @assert all(-10^-5 .< imag(fft(TF.padtosize(ones(2,2), 5, 5))) .< 10^-5) # .== 0 (numerical errors)

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

"""
    TransferFunctions.roundupcenter(arr::AbstractArray)

Calculate center index of `arr` rounded-up (i.e. `fft` center)

# Examples
```jldoctest; setup = :(using TransferFunctions: TransferFunctions as TF)
julia> TF.roundupcenter(ones(15,15))
CartesianIndex(8, 8)

julia> TF.roundupcenter(ones(16,16))
CartesianIndex(9, 9)

julia> TF.roundupcenter(ones(15,15,15))
CartesianIndex(8, 8, 8)
```
"""
function roundupcenter(axlims::NTuple{N,NTuple{2,Int}})::CartesianIndex{N} where {N}
    l = @. last(axlims) .- first(axlims)
    return CartesianIndex(@. first(axlims) + round(Int, l / 2, RoundUp))
end
roundupcenter(arr::AbstractArray) = roundupcenter(extrema.(axes(arr)))
roundupcenter(dims::Dims) = roundupcenter(Tuple((1, d) for d in dims))

"""
     TransferFunctions.exactcenter(a::AbstractArray)

Calculate the exact center index of `a`

# Examples
```jldoctest; setup = :(using TransferFunctions: TransferFunctions as TF)
julia> TF.exactcenter(ones(15,15))
(8.0, 8.0)

julia> TF.exactcenter(ones(16,16))
(8.5, 8.5)

julia> TF.exactcenter(ones(15,15,15))
(8.0, 8.0, 8.0)
```
"""
function exactcenter(arr::AbstractArray{N})::Coordinate where {N}
    axs = axes(arr)
    return @. minimum(axs) + (maximum(axs) - minimum(axs)) / 2
end

contained(arr::AbstractArray{T,N}, loc::Coordinate{N}) where {T,N} = all(loc .∈ axes(arr))

fftfreqs(sz::Dims{2}, Δ::PixelSize{2}) = ndgrid(fftfreq.(sz, 1 ./ Δ)...)
fftfreqs(sz::Dims{2}, Δ::Length) = fftfreqs(sz, fillsize(Δ, Val(2)))

posgrid(sz::Dims{2}, Δ::PixelSize{2}) = ndgrid(map(sz, Tuple(roundupcenter(sz)), Δ) do len, c, samp
    (range(1, len) .- c) .* samp
end...)
posgrid(sz::Dims{2}, Δ::Length) = posgrid(sz, fillsize(Δ, Val(2)))

struct OriginAt{N}
    origin::CartesianIndex{N}
end
(oat::OriginAt{N})(x::AbstractArray{<:Any,N}) where {N} = OA.Origin(CartesianIndex{N}(ntuple(_ -> 1, Val(N))) - oat.origin)(x)

# NOTE: https://github.com/JuliaLang/julia/issues/6733
"""
    @__FUNCTION__
Return current function name. For debugging and error message purposes.
"""
macro __FUNCTION__()
    return :($(esc(Expr(:isdefined, :var"#self#"))) ? $(esc(:var"#self#")) : nothing)
end
