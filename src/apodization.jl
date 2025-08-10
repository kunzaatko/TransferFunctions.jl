# TODO: Rewrite docs <07-08-25> 
"""
    TransferFunctions.Apodization
Module defining apodization functions.

Apodization types [`Blackman`](@ref Apodization.Blackman), [`ExactBlackman`](@ref Apodization.ExactBlackman),
[`Connes`](@ref Apodization.Connes), [`Cosine`](@ref Apodization.Cosine), [`Gaussian`](@ref Apodization.Gaussian),
[`Hamming`](@ref Apodization.Hamming), [`Welch`](@ref Apodization.Welch), [`BlackmanNuttall`](@ref
Apodization.BlackmanNuttall), [`PowerCosine`](@ref Apodization.PowerCosine), [`Triangular`](@ref
Apodization.Triangular), [`Nuttall`](@ref Apodization.Nuttall), [`SineSum`](@ref Apodization.SineSum),
[`BlackmanHarris`](@ref Apodization.BlackmanHarris), [`FlatTop`](@ref Apodization.FlatTop) and [`Hann`](@ref
Apodization.Hann)
"""
module Apodization
using InterfaceFunctions

# https://mathworld.wolfram.com/ApodizationFunction.html
# TODO: Add documentation about what is an apodization function and how it is used. "An apodization function is ... zero-phase function ... Instrument function ... Links" <18-11-24> 
"""
    ApodizationFunction{T}
Abstract type for apodization functions with return values of type `T`.
"""
abstract type ApodizationFunction{T} end
Broadcast.broadcastable(a::ApodizationFunction) = Ref(a)
@interface apodization(apo::ApodizationFunction{T}, r::Real, halfwidth::Int) where {T} = apodization(apo, r / halfwidth)
@interface apodization(apo::ApodizationFunction{T}, r::Real) where {T}

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

# TODO: Consider adding `oftype` to the implementations of `apodization` and `instrument` to improve performance <07-08-25> 

# FIX: There are multiple Triangular types of apodization that are based on where the zero is... See
# https://en.wikipedia.org/w/index.php?title=Window_function&oldid=1237444898#Triangular_window <10-09-24>. These could
# be implemented as constants of the Triangular window similar to how PowerCosine and SineSum are defined.
"""
    Triangular{T} <: ApodizationFunction{T}
Formulas:
+ zero-phase function: ``w₀(r) = 1-|r|``
+ instrument function: ``I(k) = sinc²(π k)``
"""
struct Triangular{T} <: ApodizationFunction{T} end
apodization(::Triangular{T}, r::Real) where {T} = oneunit(T) - abs(r)
instrument(::Triangular{T}, k::Real)where {T} = sinc(π * k)^2

# TODO: This should be instead implemented as a polynomial with some set parameters. It could even be by using the
# B-Splines type. <10-09-24> 

"""
    Welch{T} <: ApodizationFunction{T}
+ zero-phase function: ``w₀(r) = 1 - r²``
"""
struct Welch{T} <: ApodizationFunction{T} end
apodization(::Welch{T}, r::Real) where {T} = convert(T, 1 - r^2)
function instrument(::Welch{T}, k::Real) where {T}
    return convert(T, sinpi(2 * k) - 2π * k * cospi(2k) / ((k * 2π)^3))
end

# TODO: Should be implemented as polynomial <10-09-24> 
# FIX: Is this the correct definition of Connes? Wikipedia does not have a Connes window. <10-09-24> 
"""
    Connes{T} <: ApodizationFunction{T}
"""
struct Connes{T} <: ApodizationFunction{T} end
apodization(::Connes{T}, r::Real) where {T} = convert(T,(1 - r^2)^2)
function instrument(::Connes{T}, k::Real) where {T}
    return convert(T, 8sqrt(2π) * SpecialFunctions.besselj(5 / 2, 2π * k) / (2π * k)^(5 / 2))
end

"""
    PowerCosine{T} <: ApodizationFunction{T}
+ zero-phase function: ``w₀(r) = cos(πr/2)^α``

Instances: [`Cosine`](@ref) and [`Hann`](@ref)
"""
struct PowerCosine{T,α,E} <: ApodizationFunction{T}
    α::E
    PowerCosine{T}(α::Real) where {T} = new{T,α,typeof(α)}(α)
    PowerCosine{T,α}() where {T,α} = new{T,α,typeof(α)}(α)
end
PowerCosine(α::Real) = PowerCosine{Any}(α)
apodization(apo::PowerCosine{T}, r::Real) where {T} = convert(T,cospi(r / 2)^(apo.α))
Base.convert(::Type{ApodizationFunction{S}}, apo::PowerCosine{T,α}) where {S,T,α} = PowerCosine{S,α}()

"""
    Cosine{T} == PowerCosine{T, 1} <: ApodizationFunction{T}
+ zero-phase function: ``w₀(r) = cos(πr/2)``
+ instrument function: ``I(k) = 4cos(2k)/(π(1 - 16k²))``
"""
const Cosine{T} = PowerCosine{T, 1}

"""
    Hann{T} == PowerCosine{T,2} <: ApodizationFunction{T}
Formulas:
+ zero-phase function: ``w₀(r) = cos²(πr/2)``
"""
const Hann{T} = PowerCosine{T,2}
# instrument(apo::Hann, k::Real) = sinc(2π *  k) / (1 - 4*k^2)

"""
    SineSum{T,Cs} <: ApodizationFunction{T}
Sum of sines with coefficients `Cs` with output type of `T`

+ zero-phase function: ``w₀(r) = ∑ᴺₖ₌₀ (-1)ᵏ Cs[k] cos(πk(r + 1))``

Instances: [`Hamming`](@ref), [`Nuttall`](@ref), [`BlackmanNuttall`](@ref), [`BlackmanHarris`](@ref), [`FlatTop`](@ref), [`Blackman`](@ref) and [`ExactBlackman`](@ref)
"""
struct SineSum{T,Cs,N} <: ApodizationFunction{T}
    coefs::NTuple{N,T}
    function SineSum{T}(coefs...) where {T}
        coefs = Tuple(convert(T, c) for c in coefs)
        new{T,coefs,length(coefs)}(coefs)
    end
    function SineSum{T, Cs}() where {T, Cs}
        coefs = Tuple(convert(T,c) for c in Cs)
        new{T,Cs,length(coefs)}(coefs)
    end
end
SineSum(coefs...) = SineSum{Any}(coefs...)
apodization(ss::SineSum{T}, r::Real) where {T} = convert(T,sum((-1)^k * C * cospi(k * (r + 1)) for (k,C) in zip(0:(length(ss.coefs)-1), ss.coefs)))
Base.convert(::Type{ApodizationFunction{S}}, ss::SineSum{T}) where {S,T} = SineSum{S}(ss.coefs...)

"""
    Hamming{T} == SineSum{T,(25//46, 21//46)} <: ApodizationFunction{T}
See [`SineSum`](@ref)
"""
const Hamming{T} = SineSum{T,(25 // 46, 21 // 46)}

# function instrument(apo::Hamming, k::Real)
#     # FIX: Check whether this is correct. <10-09-24> 
#     return apo.a(25 / 46 - 16 * k^2 / 25) * sinc(2π * k) / (1 - 4 * k^2)
# end

"""
    Nuttall{T} == SineSum{T,(0.355768, 0.487396, 0.144232, 0.012604)} <: ApodizationFunction{T}
See [`SineSum`](@ref)
"""
const Nuttall{T} = SineSum{T,(0.355768, 0.487396, 0.144232, 0.012604)}

"""
    BlackmanNuttall{T} = SineSum{T,(0.3635819, 0.4891775, 0.1365995, 0.0106411)} <: ApodizationFunction{T}
See [`SineSum`](@ref)
"""
const BlackmanNuttall{T} = SineSum{T,(0.3635819, 0.4891775, 0.1365995, 0.0106411)}

"""
    BlackmanHarris{T} == SineSum{T,(0.35875, 0.48829, 0.14128, 0.01168)} <: ApodizationFunction{T}
See [`SineSum`](@ref)
"""
const BlackmanHarris{T} = SineSum{T,(0.35875, 0.48829, 0.14128, 0.01168)}

"""
    FlatTop{T} = SineSum{T,(0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947365)} <: ApodizationFunction{T}
MATLAB variant of the flat-top filter

See [`SineSum`](@ref)
"""
const FlatTop{T} = SineSum{T,(0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947365)}

"""
    Blackman{T,α} == SineSum{T,((1 - α) / 2, 1 / 2, α / 2)} <: ApodizationFunction{T}
Blackman apodization function with the coefficient `α` (Defaults to `α = 0.16`).

See [`SineSum`](@ref)
See also [`ExactBlackman`](@ref)
"""
struct Blackman{T,α,P} <: ApodizationFunction{T}
    parent::P
    function Blackman{T,α}() where {T,α}
        zero(α) < α < oftype(α,1 / 2) || throw(ArgumentError("`α` must be in `(0, 1/2)"))
        parent = SineSum{T}((1 - α) / 2, 1 / 2, α / 2)
        return new{T,α, typeof(parent)}(parent)
    end
end
Blackman{T}(α) where {T} = Blackman{T,α}()
Blackman{T}() where {T} = Blackman{T,0.16}()
Blackman(α) = Blackman{Any}(α)
Blackman() = Blackman{Any}()
apodization(apo::Blackman{T}, r::Real) where {T} = apodization(apo.parent, r)

"""
    ExactBlackman{T} == Blackman{T,683 // 4652}
See [`Blackman`](@ref)
"""
const ExactBlackman{T} = Blackman{T,683 // 4652}

"""
    Gaussian{T} <: ApodizationFunction{T}
Gaussian apodization function with a given standard deviation `σ`.

+ zero-phase function: ``w₀(r) = e^{-r²/2σ²}``
"""
struct Gaussian{T} <: ApodizationFunction{T}
    σ::T
    function Gaussian{T}(σ) where {T}
        zero(σ) < σ < oftype(σ,1 / 2) || throw(ArgumentError("`σ` must be in (0, 1/2)"))
        new{T}(σ)
    end
end
Gaussian(σ) = Gaussian{Any}(σ)
apodization(apo::Gaussian{T}, r::Real) where {T} = convert(T,exp(-r^2 / (2apo.σ^2)))
Base.convert(::Type{ApodizationFunction{S}}, apo::Gaussian) where S = Gaussian{S}(apo.σ)

# NOTE: Default to return type `Any` but enable conversion for types that do not hold any typed fields  
const Apodization_zerofield_types = (:Triangular, :Connes, :Cosine, :Welch, :Hann, :Hamming, :Nuttall, :BlackmanHarris, :BlackmanNuttall, :FlatTop, :ExactBlackman)
for type in Apodization_zerofield_types
    @eval begin
        $type() = $type{Any}()
        Base.convert(::Type{ApodizationFunction{S}}, apo::$type) where S = $type{S}()
    end
end

public apodization, instrument
public Nuttall, BlackmanNuttall, BlackmanHarris, FlatTop, Gaussian, Blackman, SineSum, Blackman, Connes, Welch, Triangular, PowerCosine, Cosine, Hamming, Hann, ExactBlackman

end
