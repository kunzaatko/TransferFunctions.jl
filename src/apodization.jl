# TODO: This should really be a separate package <26-08-24> 
# FIX: Type stability <13-12-23> 
# FIX: Domain error, when ∉ (-1,1) <13-12-23> 
# https://mathworld.wolfram.com/ApodizationFunction.html
abstract type Apodization end
Broadcast.broadcastable(a::Apodization) = Ref(a)
apodization(apo::Apodization, x::Real, halfwidth::Int) = apodization(apo, x / halfwidth)

# TODO: How to document types? <17-07-24> 

@doc raw"""
    Bartlett <: Apodization

 ``A(x) = 1-|x|/a``

``I(k) = a \mathop{sinc}²(π k a)``

Examples:
```julia
Bartlett(5) <: Apodization
```
"""
struct Bartlett <: Apodization
    a::Real
end
apodization(apo::Bartlett, x::Real) = 1 - abs(x) / apo.a
instrument(apo::Bartlett, k::Real) = apo.a * sinc(π * k * apo.a)^2


@doc raw"""
    Blackman(a::Real) <: Apodization

``A(x) = 0.42 + 0.5 \cos(\pi x/a) + 0.08 \cos(2\pi x/a)``

``I(k) = \frac{a \mathop{sinc}(2\pi k a) (0.84 - 0.36a^2k^2)}{(1 - a^2k^2)(1 - 4a^2k^2)}``

where ``a`` is the apodization parameter.
"""
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

# TODO: Documentation <17-07-24> 
@doc raw"""
    Connes(a::Real) <: Apodization
"""
struct Connes <: Apodization
    a::Real
end
apodization(apo::Connes, x::Real) = (1 - x^2 / apo.a^2)^2
function instrument(apo::Connes, k::Real)
    return 8apo.a * sqrt(2π) * SpecialFunctions.besselj(5 / 2, 2π * k * apo.a) /
           (2π * k * apo.a)^(5 / 2)
end

@doc raw"""
    Cosine(a::Real) <: Apodization

Cosine apodization function with a given amplitude factor `a`.

Formulas:
* Apodization Function: ``\cos\left(\pi \frac{x}{2a}\right)``
* Instrument Function: ``\frac{4a\cos\left(2ak\right)}{\pi\left(1 - 16a^2 k^2\right)}``
"""
struct Cosine <: Apodization
    a::Real
end
apodization(apo::Cosine, x::Real) = cospi(x / 2apo.a)
instrument(apo::Cosine, k::Real) = 4apo.a * cospi(2apo.a * k) / (π * (1 - 16apo.a^2 * k^2))

@doc raw"""
    Gaussian(σ::Real) <: Apodization

Gaussian apodization function with a given standard deviation `σ`.

Formulas:
    - Apodization Function: ``e^{-\frac{x^2}{2\sigma^2}}``
    - Instrument Function: Not implemented (Low Priority)
"""
struct Gaussian <: Apodization
    σ::Real
end
apodization(apo::Gaussian, x::Real) = exp(-x^2 / (2apo.σ^2))
function instrument(_::Gaussian, _::Real)
    throw(NotImplementedError("Low Priority (if you need it, please open an issue)"))
end

# TODO: Documentation <17-07-24> 
@doc raw"""
    Hamming(a::Real) <: Apodization
"""
struct Hamming <: Apodization
    a::Real
end
apodization(apo::Hamming, x::Real) = 27 / 50 + 23cospi(x / apo.a) / 50
function instrument(apo::Hamming, k::Real)
    return apo.a(27 / 25 - 16apo.a^2 * k^2 / 25) * sinc(2π * apo.a * k) /
           (1 - 4apo.a^2 * k^2)
end

# TODO: Documentation <17-07-24> 
@doc raw"""
    Hanning(a::Real) <: Apodization
"""
struct Hanning <: Apodization
    a::Real
end
apodization(apo::Hanning, x::Real) = cospi(x / 2apo.a)^2
instrument(apo::Hanning, k::Real) = apo.a * sinc(2π * apo.a * k) / (1 - 4apo.a^2 * k^2)

# TODO: Documentation <17-07-24> 
@doc raw"""
    Welch(a::Real) <: Apodization
"""
struct Welch <: Apodization
    a::Real
end
apodization(apo::Welch, x::Real) = 1 - x^2 / apo.a^2
function instrument(apo::Welch, k::Real)
    return apo.a * (sinpi(2 * k * apo.a) - 2π * apo.a * k * cospi(2apo.a * k)) /
           (2(apo.a * k * π)^3)
end

# TODO: There should be a mutating method and a non-mutating method that allocates the output array <26-08-24> 
# TODO: Should allow padding with some scheme from `ImageFiltering.jl` before the tapering for an unobscured data
# behaviour  <26-08-24> 
# FIX: This API should be slightly different. Right now it requires to specify the apodization. Usually the user can be
# satisfied with the default `TransferFunctions.Cosine`. The `width` parameter should be part of the `apo` instance.
# When only the width is set, the default apodization should be instantiated. The tendency should be for lower
# parameter methods be easier and just call the general methods, which are the more flexible ones. <17-07-24> 
@doc raw"""
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

# TODO: There should be a mutating method and a non-mutating method that allocates the output array <26-08-24> 
# TODO: Add documentation <21-12-23> 
@doc raw"""
    apodize(apo::Apodization, A::AbstractArray{<:Number,N}, cutoff::Real, width::Real) where {N}

Apply apodization to the input array `A`.

# Arguments
- `apo::Apodization`: The type of apodization to apply.
- `A::AbstractArray{<:Number,N}`: The input array to be apodized.
- `cutoff::Real`: The cutoff frequency for the apodization.
- `width::Real`: The width of the transition region for the apodization.
"""
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

