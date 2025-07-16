"""
`TransferFunctions.Apodization` module defines apodization functions of various formats and utility methods that are
related to apodization such as [`taperedges`](@ref) and [`apodize`](@ref).

Apodization types [`Blackman`](@ref Apodization.Blackman), [`ExactBlackman`](@ref Apodization.ExactBlackman),
[`Connes`](@ref Apodization.Connes), [`Cosine`](@ref Apodization.Cosine), [`Gaussian`](@ref Apodization.Gaussian),
[`Hamming`](@ref Apodization.Hamming), [`Welch`](@ref Apodization.Welch), [`BlackmanNuttall`](@ref
Apodization.BlackmanNuttall), [`PowerCosine`](@ref Apodization.PowerCosine), [`Triangular`](@ref
Apodization.Triangular), [`Nuttall`](@ref Apodization.Nuttall), [`SineSum`](@ref Apodization.SineSum),
[`BlackmanHarris`](@ref Apodization.BlackmanHarris), [`FlatTop`](@ref Apodization.FlatTop) and [`Hann`](@ref
Apodization.Hann)

Exports [`taperedges`](@ref), [`apodize`](@ref)
"""
module Apodization
using TransferFunctions: Size
# TODO: Use reinterpret instead of `T.` for the changes of type in the places that it is used. <09-09-24> 
using ImageFiltering: AbstractBorder, borderinstance, BorderSpecAny, Pad, Fill, Inner
using ImageFiltering: ImageFiltering as IF
using TransferFunctions: padarray

# FIX: Domain error, when ∉ (-1,1) <13-12-23> 

# https://mathworld.wolfram.com/ApodizationFunction.html
# TODO: Add documentation about what is an apodization function and how it is used. "An apodization function is ... zero-phase function ... Instrument function ... Links" <18-11-24> 
"""
    ApodizationFunction

Abstract type for apodization functions.
"""
abstract type ApodizationFunction end
Broadcast.broadcastable(a::ApodizationFunction) = Ref(a)
apodization(apo::ApodizationFunction, x::Real, halfwidth::Int) = apodization(apo, x / halfwidth)


# NOTE:  <17-05-25> 
# TODO: Same thing is done in the `slicearray.jl` `eachslice` function definition in https://github.com/JuliaLang/julia/blob/de090a92b3d564179d1fdeed7455d91c356accfb/base/slicearray.jl?plain=1#L48-L53. Use the same method. <17-05-25>

check_dims_unique(dims::Dims) = allunique(dims) || throw(ArgumentError("Dimensions in `dims` must be unique. Got `dims=$dims`."))
check_dims_bounded(A::AbstractArray{<:Any,N}, dims::Dims{M}) where {N,M} = M <= N && all(i -> i <= N, dims) || throw(DimensionMismatch("Dimensions in `dims` must be bounded `ndims(A)=$(ndims(A))`. Got `dims=$dims`."))

const Width{M} = Tuple{Size{M},Size{M}}
const SizedWidthSpec{M} = Union{Width{M},Size{M}}

# TODO: There should be a mutating method and a non-mutating method that allocates the output array <26-08-24> 
"""
    taperedges([apo=Cosine()], A, w, [border=:replicate]; dims=:)
    taperedges([apo], A, (w1,w2...,wn), [border]; dims =1:n)
    taperedges([apo], A, ((w1a,...,wna), (w1b,...,wnb)), [border]; dims=1:n)

Taper edges of width `w` of array `A` using [`apo::ApodizationFunction`](@ref ApodizationFunction). 

Widths can be same for all dimensions and directions using `w::Int`, specified separately for each of the `n` dimensions
`(w1,w2...,wn)` which means that the '_start_' and '_end_' padding will be the same. The third option is to specify the
padding fully separately for the '_start_' and '_end_' of `n` dimensions `((w1a,...,wna), (w1b,...,wnb))`. In the last
two cases, the dimensions are specified by `dims` and by default are taken as the first `n` dimensions. Border can be
`:replicate`, `:circular`, `:symmetric`, `:reflect` or [`Fill(v)`](@extref ImageFiltering
:jl:type:`ImageFiltering.Fill`).

See also [`BorderArray`](@extref ImageFiltering :std:label:`BorderArray`), [`ApodizationFunction`](@ref)
"""
taperedges(A::AbstractArray, args...; kwargs...) = taperedges(Cosine(), A, args...; kwargs...)
taperedges(apo::ApodizationFunction, A::AbstractArray, w::SizedWidthSpec{M}, border="replicate"; dims=Dims(1:M), kwargs...) where {M} = _taperedges(apo, A, w, border, dims; check_bounded=false, check_unique=false, kwargs...)
taperedges(apo::ApodizationFunction, A::AbstractArray, w, border="replicate"; dims=:, kwargs...) = _taperedges(apo, A, w, border, dims; kwargs...)

const WidthSpec = Union{<:Integer,<:Size,Tuple{Size{N},Size{N}}} where {N}

_taperedges( # STEP 1a: Fill single dim
    apo::ApodizationFunction,
    A::AbstractArray,
    w::WidthSpec,
    border::Any,
    dims::Integer; kwargs...) = _taperedges(apo, A, w, border, Dims(dims); check_unique=false, kwargs...)

_taperedges( # STEP 1b: Fill in dims for `Colon`
    apo::ApodizationFunction,
    A::AbstractArray{<:Number,N},
    w::WidthSpec,
    border::Any,
    dims::Colon; kwargs...) where {N} = _taperedges(apo, A, w, border, Dims(1:N); check_bounded=false, check_unique=false, kwargs...)

function _taperedges( # STEP 2a: Fill same `width` for all dimensions
    apo::ApodizationFunction,
    A::AbstractArray,
    w::Int,
    border::Any,
    dims::Dims{M}; kwargs...
) where {M}
    ws = ntuple(_ -> w, Val(M))
    return _taperedges(apo, A, (ws, ws), border, dims; kwargs...)
end

function _taperedges( # STEP 2b: Fill same left and right `width`s
    apo::ApodizationFunction,
    A::AbstractArray,
    w::Size{M},
    border::Any,
    dims::Dims{M}; kwargs...
) where {M}
    return _taperedges(apo, A, (w, w), border, dims; kwargs...)
end

function _taperedges( # STEP 3: Create a border instance
    apo::ApodizationFunction,
    A::AbstractArray,
    w::Width{M},
    border::AbstractString,
    dims::Dims{M}; kwargs...
) where {M}
    return _taperedges(apo, A, w, borderinstance(border), dims; kwargs...)
end

function _taperedges( # STEP 4: Set the border sizes
    apo::ApodizationFunction,
    A::AbstractArray,
    w::Width{M},
    border::BorderSpecAny,
    dims::Dims{M}; kwargs...
) where {M}
    return _taperedges(apo, A, w, specify_border(A, border, w, dims), dims; kwargs...)
end

function specify_border(::AbstractArray{<:Any,N}, border::BorderSpecAny, w::Width{M}, dims::Dims{M}) where {M,N}
    lw, rw = (ntuple(i -> i in dims ? w[lr][findfirst(j -> j == i, dims)] : 0, Val(N)) for lr in 1:2)
    return _specify_border(border, lw, rw)
end
_specify_border(border::Pad, lw::Size{N}, rw::Size{N}) where {N} = Pad(border.style, lw, rw)
_specify_border(border::Fill, lw::Size{N}, rw::Size{N}) where {N} = Fill(border.value, lw, rw)
_specify_border(::Inner, lw::Size{N}, rw::Size{N}) where {N} = Inner(lw, rw)
_specify_border(border::Union{IF.NA,IF.NoPad}, args...) = throw(ArgumentError("Border must be one of `Pad`, `Fill` and `Inner`. Got `$border`."))

# TODO: Perhaps there could be an argument to make the array odd sized for the Fourier transform. Since we do not have
# to have 0 at both edges. For the signal to be periodic, only one edge to be 0 is sufficient. An odd size is beneficial
# for a Fourier transform. <10-09-24> 
function _taperedges( # FINAL
    apo::ApodizationFunction,
    A::AbstractArray,
    ws::Width{M},
    border::AbstractBorder,
    dims::Dims{M}; check_unique=true, check_bounded=true
) where {M}
    check_unique && check_dims_unique(dims)
    check_bounded && check_dims_bounded(A, dims)

    A = padarray(A, border)

    for (dim, low, high) in zip(dims, ws[1], ws[2])
        lowedge = range(-1, 0; length=low + 1)[begin:(end-1)]
        for (e, i) in enumerate(firstindex(axes(A, dim)):(firstindex(axes(A, dim))+(low-1)))
            selectdim(A, dim, i) .*= apodization(apo, lowedge[e])
        end
        highedge = range(0, 1; length=high + 1)[(begin+1):end]
        @assert length(highedge) == length((lastindex(axes(A, dim))-(high-1)):lastindex(axes(A, dim)))
        for (e, i) in enumerate((lastindex(axes(A, dim))-(high-1)):lastindex(axes(A, dim)))
            selectdim(A, dim, i) .*= apodization(apo, highedge[e])
        end
    end
    return A
end

# TODO: There should be a mutating method and a non-mutating method that allocates the output array <26-08-24> 
# TODO: Default values for the width and cutoff <04-05-25> 
"""
    apodize(apo::ApodizationFunction, A, cutoff::Real, width::Real)

Apply apodization to the input array `A`.

# Arguments
- `apo::ApodizationFunction`: The type of apodization to apply.
- `cutoff::Real`: The cutoff frequency for the apodization.
- `width::Real`: The width of the transition region for the apodization.
"""
function apodize(
    apo::ApodizationFunction,
    A::AbstractArray{<:Number,N},
    cutoff::Real, # Tuple{NTuple{M,Int},NTuple{M,Int}}, # FIX: Generalize to non-symmetric cut-offs <21-12-23> 
    width::Real, # FIX: Generalize to non-symmetric widths <21-12-23> 
    dims::Dims=Tuple(1:N) # FIX: Abstract over dimensions <22-12-23> 
) where {N}
    # PERF: Should be done with no allocation... This is the KISS solution  
    # FIX: Work for all dimensions  
    rs = [hypot(abs(x), abs(y)) for x in fftfreq(size(A, 1), size(A, 1)), y in fftfreq(size(A, 2), size(A, 2))]
    coefs = ones(eltype(A), size(rs)...)
    coefs[rs.>=cutoff] .= 0
    coefs[rs.<=(cutoff-width)] .= 1
    coefs[cutoff.>rs.>(cutoff-width)] .= map(r -> apodization(apo, r - (cutoff - width), width), rs[cutoff.>rs.>(cutoff-width)])
    return A .* coefs
end

# https://www.gaussianwaves.com/2020/09/window-function-figure-of-merits/
# https://www.mathworks.com/help/signal/ref/bohmanwin.html
# https://www.mathworks.com/help/signal/windows.html?s_tid=CRUX_lftnav
# https://en.wikipedia.org/wiki/Window_function
# https://mathworld.wolfram.com/ApodizationFunction.html
# TODO: Warn when the apodization function is called out of the interval (-1,1) where it makes sense to sample it <18-11-24> 
# TODO: Must implement the Tuckey, Planck-taper and other windows that are usable for apodization and not edge tapering
# <10-09-24> 
# TODO: There must be a some plots in the documentation... Otherwise this is a waste <10-09-24> 
# TODO: A general way of implementation of linear combination types of apodization functions <10-09-24> 
# TODO:  A general B-Spline type window.. Triangular and Parzen windows would be part of these.
# https://en.wikipedia.org/wiki/Window_function  <10-09-24> 
# TODO: These should be implemented: https://en.wikipedia.org/wiki/Window_function instead of the ones in WolframAlpha's
# page <02-09-24> 
# TODO: How to document types? <17-07-24> 

# FIX: There are multiple Triangular types of apodization that are based on where the zero is... See
# https://en.wikipedia.org/w/index.php?title=Window_function&oldid=1237444898#Triangular_window <10-09-24>. These could
# be implemented as constants of the Triangular window similar to how PowerCosine and SineSum are defined.
"""
    Triangular <: ApodizationFunction

Formulas:
+ zero-phase function: ``w₀(r) = 1-|r|``
+ instrument function: ``I(k) = sinc²(π k)``
"""
struct Triangular <: ApodizationFunction end
apodization(::Triangular, r::Real) = 1 - abs(r)
instrument(::Triangular, k::Real) = sinc(π * k)^2

# TODO: This should be instead implemented as a polynomial with some set parameters. It could even be by using the
# B-Splines type. <10-09-24> 

"""
    Welch <: ApodizationFunction

+ zero-phase function: ``w₀(r) = 1 - r²``
"""
struct Welch <: ApodizationFunction end

apodization(::Welch, r::Real) = 1 - r^2
function instrument(::Welch, k::Real)
    return sinpi(2 * k) - 2π * k * cospi(2k) / ((k * 2π)^3)
end

# TODO: Should be implemented as polynomial <10-09-24> 
# FIX: Is this the correct definition of Connes? Wikipedia does not have a Connes window. <10-09-24> 
"""
    Connes <: ApodizationFunction
"""
struct Connes <: ApodizationFunction end
apodization(::Connes, r::Real) = (1 - r^2)^2
function instrument(::Connes, k::Real)
    return 8sqrt(2π) * SpecialFunctions.besselj(5 / 2, 2π * k) / (2π * k)^(5 / 2)
end

"""
    PowerCosine{α<:Real} <: ApodizationFunction

+ zero-phase function: ``w₀(r) = cos(πr/2)^α``

Instances: [`Cosine`](@ref) and [`Hann`](@ref)
"""
struct PowerCosine{α} <: ApodizationFunction
    function PowerCosine{α}() where {α}
        α isa Real || throw(ArgumentError("α must be a real number"))
        new{α}()
    end
end
apodization(::PowerCosine{α}, r::Real) where {α} = cospi(r / 2)^(α)
# instrument(::PowerCosine{α}, k::Real) = 4cospi(2k) / (π * (1 - 16k^2)) # TODO: <02-09-24> 

"""
    Cosine == PowerCosine{1} <: ApodizationFunction

+ zero-phase function: ``w₀(r) = cos(πr/2)``
+ instrument function: ``I(k) = 4cos(2k)/(π(1 - 16k²))``
"""
const Cosine = PowerCosine{1}

"""
    Hann == PowerCosine{2} <: ApodizationFunction

Formulas:
+ zero-phase function: ``w₀(r) = cos²(πr/2)``
"""
const Hann = PowerCosine{2}
# instrument(apo::Hann, k::Real) = sinc(2π *  k) / (1 - 4*k^2)

"""
    SineSum{N,Cs,T} <: ApodizationFunction
Sum of `N` sines apodization function with output type of `T`

+ zero-phase function: ``w₀(r) = ∑ᴺₖ₌₀ (-1)ᵏ Cs[k] cos(πk(r + 1))``

Instances: [`Hamming`](@ref), [`Nuttall`](@ref), [`BlackmanNuttall`](@ref), [`BlackmanHarris`](@ref), [`FlatTop`](@ref), [`Blackman`](@ref) and [`ExactBlackman`](@ref)
"""
struct SineSum{N,Cs,T<:Real} <: ApodizationFunction
    coefs::NTuple{N,T}
    SineSum{N,Cs}() where {N,Cs} = SineSum{N,Cs,Float64}()
    function SineSum{N,Cs,T}() where {N,T,Cs}
        N isa Int && N > 0 || throw(ArgumentError("N must be a positive integer"))
        coefs = try
            convert(NTuple{N,T}, Cs)
        catch
            throw(ArgumentError("In constructor of `SineSum{N, T, Cs}`, `Cs` must be convertible to `NTuple{N, T}`."))
        end
        return new{N,Cs,T}(coefs)
    end
end
function SineSum(Cs::Tuple{T,Vararg{T}}) where {T<:Real}
    N = length(Cs)
    return SineSum{N,Cs,T}()
end
apodization(ss::SineSum{N,Cs,T}, r::Real) where {N,T,Cs} =
    sum(zip(0:(N-1), ss.coefs)) do (k, C)
        (-1)^k * C * cospi(k * (r + 1))
    end

"""
    Hamming == SineSum{2,(25//46, 21//46)} <: ApodizationFunction
"""
const Hamming = SineSum{2,(25 // 46, 21 // 46)}

function instrument(apo::Hamming, k::Real)
    # FIX: Check whether this is correct. <10-09-24> 
    return apo.a(25 / 46 - 16 * k^2 / 25) * sinc(2π * k) / (1 - 4 * k^2)
end

"""
    Nuttall} == SineSum{4,(0.355768, 0.487396, 0.144232, 0.012604)} <: ApodizationFunction
"""
const Nuttall = SineSum{4,(0.355768, 0.487396, 0.144232, 0.012604)}

"""
    BlackmanNuttall == SineSum{4,(0.3635819, 0.4891775, 0.1365995, 0.0106411) } <: ApodizationFunction
"""
const BlackmanNuttall = SineSum{4,(0.3635819, 0.4891775, 0.1365995, 0.0106411)}

"""
    BlackmanHarris == SineSum{4,(0.35875, 0.48829, 0.14128, 0.01168)} <: ApodizationFunction
"""
const BlackmanHarris = SineSum{4,(0.35875, 0.48829, 0.14128, 0.01168)}

"""
    FlatTop == SineSum{5,(0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947365)} <: ApodizationFunction

MATLAB variant of the flat-top filter
"""
const FlatTop = SineSum{5,(0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947365)}

"""
    Blackman{α,T} == SineSum{3,((1 - α) / 2, 1 / 2, α / 2),T} <: ApodizationFunction

Default is `α = 0.16`.

See also [`ExactBlackman`](@ref)
"""
struct Blackman{α,T<:Real} <: ApodizationFunction
    parent::SineSum{3,Cs,T} where {Cs}
    Blackman{α}() where {α} = Blackman{α,Float64}()
    Blackman() = Blackman{0.16}()
    function Blackman{α,T}() where {α,T}
        αT = try
            convert(T, α)
        catch
            throw(ArgumentError("In constructor of `Blackman{T, α}`, `α` must be convertible to `T`."))
        end
        0 < αT < 1 / 2 || throw(ArgumentError("`α` must be in `(0, 1/2)"))
        return new{αT,T}(SineSum{3,((1 - αT) / 2, 1 / 2, αT / 2),T}())
    end
end
Blackman(v::T) where {T<:Real} = Blackman{v,T}()

"""
    ExactBlackman == Blackman{683 // 4652}
"""
const ExactBlackman = Blackman{683 // 4652}

apodization(apo::Blackman, r::Real) = apodization(apo.parent, r)

# FIX: Is there a way to do this properly? <10-09-24> 
# Core.isa(::Blackman{α}, ::Type{SineSum{3,Cs}}) where {α,Cs} = Cs == ((1 - α) / 2, 1 / 2, α / 2)

# function instrument(apo::Blackman, k::Real)
#     return sinc(2π * k) * (21 / 25 - 9 * k^2 / 25) / ((1 - k^2) * (1 - 4 * k^2))
# end


"""
    Gaussian{T<:Real} <: ApodizationFunction

Gaussian apodization function with a given standard deviation `σ`.

+ zero-phase function: ``w₀(r) = e^{-r²/2σ²}``
"""
struct Gaussian{T<:Real} <: ApodizationFunction
    σ::T
    function Gaussian(σ::T) where {T<:Real}
        0 < σ < 1 / 2 || throw(ArgumentError("`σ` must be in (0, 1/2)"))
        new{T}(σ)
    end
end
apodization(apo::Gaussian, r::Real) = exp(-r^2 / (2apo.σ^2))
instrument(::Gaussian, ::Real) = throw(ErrorException("Not yet implemented: Low Priority (if you need it, please open an issue)"))

public apodization, instrument
public Nuttall, BlackmanNuttall, BlackmanHarris, FlatTop, Gaussian, Blackman, SineSum, Blackman, Connes, Welch, Triangular, PowerCosine, Cosine, Hamming, Hann, ExactBlackman
export taperedges, apodize

end
