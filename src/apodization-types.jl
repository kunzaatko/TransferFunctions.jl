# https://www.gaussianwaves.com/2020/09/window-function-figure-of-merits/
# https://www.mathworks.com/help/signal/ref/bohmanwin.html
# https://www.mathworks.com/help/signal/windows.html?s_tid=CRUX_lftnav
# https://en.wikipedia.org/wiki/Window_function
# https://mathworld.wolfram.com/ApodizationFunction.html
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
@doc raw"""
    Triangular <: Apodization

Formulas:
+ zero-phase function: ``w₀(r) = 1-|r|``
+ instrument function: ``I(k) = \mathop{sinc}²(π k)``
"""
struct Triangular <: Apodization end
apodization(::Triangular, r::Real) = 1 - abs(r)
instrument(::Triangular, k::Real) = sinc(π * k)^2

# TODO: This should be instead implemented as a polynomial with some set parameters. It could even be by using the
# B-Splines type. <10-09-24> 
@doc raw"""
    Welch() <: Apodization

+ zero-phase function: ``w₀(r) = 1 - r²``
"""
struct Welch <: Apodization end

apodization(::Welch, r::Real) = 1 - r^2
function instrument(::Welch, k::Real)
    return sinpi(2 * k) - 2π * k * cospi(2k) / ((k * 2π)^3)
end

# TODO: Should be implemented as polynomial <10-09-24> 
# FIX: Is this the correct definition of Connes? Wikipedia does not have a Connes window. <10-09-24> 
# TODO: Documentation <17-07-24> 
@doc raw"""
    Connes(a::Real) <: Apodization
"""
struct Connes <: Apodization end
apodization(::Connes, r::Real) = (1 - r^2)^2
function instrument(::Connes, k::Real)
    return 8sqrt(2π) * SpecialFunctions.besselj(5 / 2, 2π * k) / (2π * k)^(5 / 2)
end

@doc raw"""
    PowerCosine{α<:Real} <: Apodization

+ zero-phase function: ``w₀(r) = cos(πr/2)^α``
"""
struct PowerCosine{α} <: Apodization
    function PowerCosine{α}() where {α}
        α isa Real || throw(ArgumentError("α must be a real number"))
        new{α}()
    end
end
apodization(::PowerCosine{α}, r::Real) where {α} = cospi(r / 2)^(α)
# instrument(::PowerCosine{α}, k::Real) = 4cospi(2k) / (π * (1 - 16k^2)) # TODO: <02-09-24> 

@doc raw"""
    Cosine == PowerCosine{1} <: Apodization

+ zero-phase function: ``w₀(r) = \cos(πr/2)``
+ instrument function: ``I(k) = 4\cos(2k)/(π(1 - 16k²))``
"""
const Cosine = PowerCosine{1}

@doc raw"""
Hann == PowerCosine{2} <: Apodization

Formulas:
+ zero-phase function: ``w₀(r) = \cos²(πr/2)``
"""
const Hann = PowerCosine{2}
# instrument(apo::Hann, k::Real) = sinc(2π *  k) / (1 - 4*k^2)

@doc raw"""
    SineSum{N::UInt, Cs::NTuple{N, Real}} <: Apodization

+ zero-phase function: ``w₀(r) = ∑ᴺₖ₌₀ (-1)ᵏ Cs[k] \cos(πk(r + 1))``
"""
struct SineSum{N,Cs} <: Apodization
    coefficients::NTuple{N,Real}
    function SineSum{N,Cs}() where {N,Cs}
        N isa Int && N >= 0 || throw(ArgumentError("N must be a positive integer"))
        Cs isa NTuple{N,Real} || throw(ArgumentError("Cs must be $N real coefficients"))
        new{N,Cs}()
    end
end
apodization(::SineSum{N,Cs}, r::Real) where {N,Cs} =
    sum(zip(0:(N-1), Cs)) do (k, C)
        (-1)^k * C * cospi(k * (r + 1))
    end

@doc raw"""
    Hamming == SineSum{2, (25//46, 21//46) } <: Apodization
"""
const Hamming = SineSum{2,(25 // 46, 21 // 46)}

function instrument(apo::Hamming, k::Real)
    # FIX: Check whether this is correct. <10-09-24> 
    return apo.a(25 / 46 - 16 * k^2 / 25) * sinc(2π * k) / (1 - 4 * k^2)
end

@doc raw"""
    Nuttall == SineSum{4, (0.355768, 0.487396, 0.144232, 0.012604) } <: Apodization
"""
const Nuttall = SineSum{4,(0.355768, 0.487396, 0.144232, 0.012604)}

@doc raw"""
    BlackmanNuttall == SineSum{4, (0.3635819, 0.4891775, 0.1365995, 0.0106411) } <: Apodization
"""
const BlackmanNuttall = SineSum{4,(0.3635819, 0.4891775, 0.1365995, 0.0106411)}

@doc raw"""
    BlackmanHarris == SineSum{4, (0.35875, 0.48829, 0.14128, 0.01168) } <: Apodization
"""
const BlackmanHarris = SineSum{4,(0.35875, 0.48829, 0.14128, 0.01168)}

@doc raw"""
    FlatTop == SineSum{5, (0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947365) } <: Apodization

MATLAB varaint of the flat-top filter
"""
const FlatTop = SineSum{5,(0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947365)}

@doc raw"""
    Blackman{α::Real} == SineSum{3, ((1 - α) / 2, 1 / 2, α / 2) } <: Apodization

`α = 0.16`
"""
struct Blackman{α} <: Apodization
    sine_sum::SineSum{3}
    function Blackman{α}() where {α}
        α isa Real || throw(ArgumentError("α must be a real number"))
        α ∈ OpenInterval(0, 1 / 2) || throw(ArgumentError("`α` must be in (0, 1/2)"))
        return new{α}(SineSum{3,((1 - α) / 2, 1 / 2, α / 2)}())
    end
end
Blackman() = Blackman{0.16}()
const ExactBlackman = Blackman{683 // 4652}
apodization(apo::Blackman{α}, r::Real) where {α} = apodization(apo.sine_sum, r)
# FIX: Is there a way to do this properly? <10-09-24> 
# Core.isa(::Blackman{α}, ::Type{SineSum{3,Cs}}) where {α,Cs} = Cs == ((1 - α) / 2, 1 / 2, α / 2)

# function instrument(apo::Blackman, k::Real)
#     return sinc(2π * k) * (21 / 25 - 9 * k^2 / 25) / ((1 - k^2) * (1 - 4 * k^2))
# end


@doc raw"""
    Gaussian(σ::Real) <: Apodization

Gaussian apodization function with a given standard deviation `σ`.

+ zero-phase function: ``w₀(r) = e^{-r²/2σ²}``
"""
struct Gaussian <: Apodization
    σ::Real
    function Gaussian(σ::Real)
        σ ∉ OpenInterval(0, 1 / 2) && throw(ArgumentError("`σ` must be in (0, 1/2)"))
        new(σ)
    end
end
apodization(apo::Gaussian, r::Real) = exp(-r^2 / (2apo.σ^2))
function instrument(::Gaussian, ::Real)
    throw(NotImplementedError("Low Priority (if you need it, please open an issue)"))
end
