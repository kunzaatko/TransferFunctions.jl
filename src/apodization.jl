module Apodization
using TransferFunctions: Size
# TODO: Use reinterpret instead of `T.` for the changes of type in the places that it is used. <09-09-24> 
using ImageFiltering: AbstractBorder, borderinstance, BorderSpecAny
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

# FIX: The call stack must be rewritten to construct the arguments of the function from the beginning <02-09-24> 

# FIX: docs: [`Fill(v)`](@extref Julia :jl:type:`ImageFiltering.Fill`), [`BorderType`](@extref Julia :jl:type:`ImageFiltering.BorderType`)), [`ImageFiltering.jl`](@extref).
# TODO: There should be a mutating method and a non-mutating method that allocates the output array <26-08-24> 
# TODO: Should allow padding with some scheme from `ImageFiltering.jl` before the tapering for an unobscured data
# behaviour  <26-08-24> 
# FIX: This API should be slightly different. Right now it requires to specify the apodization. Usually the user can be
# satisfied with the default `TransferFunctions.Cosine`. The `width` parameter should be part of the `apo` instance.
# When only the width is set, the default apodization should be instantiated. The tendency should be for lower
# parameter methods be easier and just call the general methods, which are the more flexible ones. <17-07-24> 
# TODO: Document: how the widths arguments work <02-09-24> 
"""
    taperedges([apo=Cosine()], A, w, [border=:replicate]; dims=:)
    taperedges([apo], A, (w1,w2...,wn), [border]; dims =1:n)
    taperedges([apo], A, ((w1a,...,wna), (w1b,...,wnb)), [border]; dims=1:n)

Taper edges of width `w` of `A` using [`apo::ApodizationFunction`](@ref ApodizationFunction).
Width can be the same `w::Int` for each dimension, specified separately for `n` dimensions `(w1,w2...,wn)`,  or
specified separately for the _start_ and _end_ of `n` dimensions `((w1a,...,wna), (w1b,...,wnb))`. In the latter two
cases, the dimensions are specified by `dims` and by default taken as the first `n` dimensions. Border can be
`:replicate`, `:circular`, `:symmetric`, `:reflect` or `Fill(v)` from `ImageFiltering.jl`.

See also `BorderType`, [`ApodizationFunction`](@ref)
"""
function taperedges( # STEP 1A: Fill the apodization type
    A::AbstractArray{<:Number,N},
    args...;
    kwargs...
) where {N}
    return taperedges(Cosine(), A, args...; kwargs...)
end
function taperedges( # STEP 2A: Fill from `width` single width
    apo::ApodizationFunction,
    A::AbstractArray{<:Number,N},
    w::Int,
    args...;
    # FIX: Instead of this, should be default `Colon()` <02-09-24> 
    dims::Dims{M}=Tuple(1:N), # TODO: Document by default taper along all the dimensions <03-09-24> 
    kwargs...
) where {N,M}
    # TODO: Use `Val` for the `fill` <30-04-25> 
    # TODO: Test whether this works for permuted order of kwargs... I.e. whether dims must be supplied as the first
    # argument or not. Otherwise, it must be done by testing if kwargs has dims in it... <02-09-24> 
    return taperedges(apo, A, (Tuple(fill(w, M)), Tuple(fill(w, M))), args...; dims, kwargs...)
end
function taperedges( # STEP 2B: Fill from `width`s for each dimension
    apo::ApodizationFunction,
    A::AbstractArray{<:Number,N},
    ws::Size{M},
    args...;
    kwargs...
) where {N,M}
    return taperedges(apo, A, (ws, ws), args...; kwargs...)
end
function taperedges( # STEP 3: Fill in the default `border`
    apo::ApodizationFunction,
    A::AbstractArray{<:Number,N},
    ws::Tuple{Size{M},Size{M}},
    args...;
    kwargs...
) where {N,M}
    return taperedges(apo, A, ws, "replicate", args...; kwargs...)
end
function taperedges( # STEP 4: Create a border instance
    apo::ApodizationFunction,
    A::AbstractArray{<:Number,N},
    ws::Tuple{Size{M},Size{M}},
    border::AbstractString,
    args...;
    kwargs...
) where {N,M}
    return taperedges(apo, A, ws, borderinstance(border), args...; kwargs...)
end
function taperedges( # STEP 5: Set the border sizes
    apo::ApodizationFunction,
    A::AbstractArray{<:Number,N},
    ws::Tuple{Size{M},Size{M}},
    border::BorderSpecAny,
    args...;
    dims=Tuple(1:M),
    kwargs...
) where {N,M}
    # FIX: Must be tested here and in the final because otherwise we would index out of bounds <09-09-24> 
    (M > N || maximum(dims) > N) && throw(ArgumentError("The number of dimensions must be less than or equal to the number of axes."))

    concrete_border = full_padding_border(border, ws, dims, N)
    return taperedges(apo, A, ws, concrete_border, args...; kwargs...)
end

# TODO: Test this <09-09-24> 
function full_padding_border(border::BorderSpecAny, ws::Tuple{Size{M},Size{M}}, dims::Dims{M}, ndims::Int) where {M}
    all_left_widths, all_right_widths = zeros(Int, ndims), zeros(Int, ndims)
    foreach(dims, ws[1], ws[2]) do dim, w_left, w_right
        all_left_widths[dim] = w_left
        all_right_widths[dim] = w_right
    end
    all_left_widths, all_right_widths = Tuple(all_left_widths), Tuple(all_right_widths)
    if border isa IF.Pad
        return IF.Pad(border.style, all_left_widths, all_right_widths)
    elseif border isa IF.Fill
        return IF.Fill(border.value, all_left_widths, all_right_widths)
    elseif border isa IF.Inner
        return IF.Inner(all_left_widths, all_right_widths)
    else
        throw(ErrorException("`NA` and `NoPad` borders should not occur here. Type is $(typeof(border))."))
    end
end

# TODO: Perhaps there could be an argument to make the array odd sized for the Fourier transform. Since we do not have
# to have 0 at both edges. For the signal to be periodic, only one edge to be 0 is sufficient. An odd size is beneficial
# for a Fourier transform. <10-09-24> 
function taperedges( # FINAL # TODO: Instead of this should be something like `_taperedges` function <02-09-24> 
    apo::ApodizationFunction,
    A::AbstractArray{<:Number,N},
    ws::Tuple{Size{M},Size{M}},
    border::AbstractBorder;
    dims::Dims{M}=Tuple(1:M) # TODO: Document that the number of dimensions is by default taken as the first M 
    # dimensions where M is the number of widths supplied <kunzaatko> 
) where {N,M}
    # TODO: How can one handle the dimensions and the various types that can define them (such as Colon())... Ask on
    # discourse and implement. <02-09-24> 
    (M > N || maximum(dims) > N) && throw(ArgumentError("The number of dimensions must be less than or equal to the number of axes."))

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
# TODO: Add documentation <21-12-23> 
"""
    apodize(apo::ApodizationFunction, A::AbstractArray{<:Number,N}, cutoff::Real, width::Real) where {N}

Apply apodization to the input array `A`.

# Arguments
- `apo::ApodizationFunction`: The type of apodization to apply.
- `A::AbstractArray{<:Number,N}`: The input array to be apodized.
- `cutoff::Real`: The cutoff frequency for the apodization.
- `width::Real`: The width of the transition region for the apodization.
"""
function apodize(
    apo::ApodizationFunction,
    A::AbstractArray{<:Number,N},
    cutoff::Real, # Tuple{NTuple{M,Int},NTuple{M,Int}}, # FIX: Generalize to non-symmetric cut-offs <21-12-23> 
    width::Real, # FIX: Generalize to non-symmetric widths <21-12-23> 
    dims::NTuple{M,Int}=Tuple(1:N) # FIX: Abstract over dimensions <22-12-23> 
) where {N,M}
    # PERF: Should be done with no allocation... This is the KISS solution  
    # FIX: Work for all dimensions  
    rs = [hypot(abs(x), abs(y)) for x in fftfreq(size(A, 1), size(A, 1)), y in fftfreq(size(A, 2), size(A, 2))]
    coefs = ones(eltype(A), size(rs)...)
    coefs[rs.>=cutoff] .= 0
    coefs[rs.<=(cutoff-width)] .= 1
    coefs[cutoff.>rs.>(cutoff-width)] .= map(r -> apodization(apo, r - (cutoff - width), width), rs[cutoff.>rs.>(cutoff-width)])
    A .* coefs
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
+ instrument function: ``I(k) = \\mathop{sinc}²(π k)``
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

+ zero-phase function: ``w₀(r) = \\cos(πr/2)^α``

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

+ zero-phase function: ``w₀(r) = \\cos(πr/2)``
+ instrument function: ``I(k) = 4\\cos(2k)/(π(1 - 16k²))``
"""
const Cosine = PowerCosine{1}

"""
    Hann == PowerCosine{2} <: ApodizationFunction

Formulas:
+ zero-phase function: ``w₀(r) = \\cos²(πr/2)``
"""
const Hann = PowerCosine{2}
# instrument(apo::Hann, k::Real) = sinc(2π *  k) / (1 - 4*k^2)

"""
    SineSum{N,T,Cs::NTuple{N,T}} <: ApodizationFunction

+ zero-phase function: ``w₀(r) = ∑ᴺₖ₌₀ (-1)ᵏ Cs[k] \\cos(πk(r + 1))``

Instances: [`Hamming`](@ref), [`Nuttall`](@ref), [`BlackmanNuttall`](@ref), [`BlackmanHarris`](@ref), [`FlatTop`](@ref), [`Blackman`](@ref) and [`ExactBlackman`](@ref)
"""
struct SineSum{N,T<:Real,Cs} <: ApodizationFunction
    function SineSum{N,T,Cs}() where {N,T,Cs}
        N isa Int && N > 0 || throw(ArgumentError("N must be a positive integer"))
        # FIX: This does not work. It does not throw the error if the type of `Cs` is not convertible to `NTuple{N, T}` <30-04-25> 
        Ts = convert(NTuple{N,T}, Cs)
        # throw(ArgumentError("In constructor of `SineSum{N, T, Cs}`, `Cs` must be convertible to `NTuple{N, T}`"))
        # rethrow(e)
        Ts isa NTuple{N,T} || throw(ArgumentError("Cs must be $N real coefficients"))
        new{N,T,Ts}()
    end
end
function SineSum(Cs::NTuple{N, T}) where {N,T<:Real}
    return SineSum{N,T,Cs}()
end
apodization(::SineSum{N,T,Cs}, r::Real) where {N,T,Cs} =
    sum(zip(0:(N-1), Cs)) do (k, C)
        (-1)^k * C * cospi(k * (r + 1))
    end

"""
    Hamming{T} == SineSum{2,T,(25//46, 21//46)} <: ApodizationFunction
"""
const Hamming{T} = SineSum{2,T,(25 // 46, 21 // 46)}
const Hamming = Hamming{Float64}

function instrument(apo::Hamming, k::Real)
    # FIX: Check whether this is correct. <10-09-24> 
    return apo.a(25 / 46 - 16 * k^2 / 25) * sinc(2π * k) / (1 - 4 * k^2)
end

"""
    Nuttall{T} == SineSum{4,T,(0.355768, 0.487396, 0.144232, 0.012604) } <: ApodizationFunction
"""
const Nuttall{T} = SineSum{4,T,(0.355768, 0.487396, 0.144232, 0.012604)}
const Nuttall = Nuttall{Float64}

"""
    BlackmanNuttall{T} == SineSum{4,T,(0.3635819, 0.4891775, 0.1365995, 0.0106411) } <: ApodizationFunction
"""
const BlackmanNuttall{T} = SineSum{4,T,(0.3635819, 0.4891775, 0.1365995, 0.0106411)}
const BlackmanNuttall = BlackmanNuttall{Float64}

"""
    BlackmanHarris{T} == SineSum{4,T,(0.35875, 0.48829, 0.14128, 0.01168) } <: ApodizationFunction
"""
const BlackmanHarris{T} = SineSum{4,T,(0.35875, 0.48829, 0.14128, 0.01168)}
const BlackmanHarris = BlackmanHarris{Float64}

"""
    FlatTop{T} == SineSum{5,T,(0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947365) } <: ApodizationFunction

MATLAB variant of the flat-top filter
"""
const FlatTop{T} = SineSum{5,T,(0.21557895, 0.41663158, 0.277263158, 0.083578947, 0.006947365)}
const FlatTop = FlatTop{Float64}

"""
    Blackman{T, α::T} == SineSum{3,T,((1 - α) / 2, 1 / 2, α / 2) } <: ApodizationFunction

Default is `α = 0.16`.

See also [`ExactBlackman`](@ref)
"""
struct Blackman{T<:Real,α} <: ApodizationFunction
    parent::SineSum{3,T}
    Blackman() = Blackman{Float64}()
    Blackman{T}() where {T} = Blackman{T,0.16}()
    function Blackman{T,α}() where {T,α}
        αT = convert(T, α)
        0 < αT < 1 / 2 || throw(ArgumentError("`α` must be in `(0, 1/2)"))
        return new{T,αT}(SineSum{3,T,((1 - αT) / 2, 1 / 2, αT / 2)}())
    end
end
Blackman(v::T) where {T<:Real} = Blackman{T,v}()

"""
    ExactBlackman{T<:Real} == Blackman{T,683 // 4652}
"""
const ExactBlackman{T} = Blackman{T,683 // 4652}
const ExactBlackman = ExactBlackman{Float64}

apodization(apo::Blackman{α}, r::Real) where {α} = apodization(apo.parent, r)

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
