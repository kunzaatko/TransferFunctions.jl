using Base: Indices
using Roots, InterfaceFunctions
using ComponentArrays
using TransferFunctions.Apodization

# TODO: Should be handled by interface functions whether the trait is implemented <20-08-25> 
"""
    intensity(psf::PointSpreadFunction, r::Length...)
The intensity of the [`PointSpreadFunction`](@ref) at the given coordinates.

For different [`PointSpreadFunctionSymmetry`](@ref) of the `psf`, it has to be defined for different arguments.
- `symmetry(psf) == ZAxisRadialSymmetry() && psf isa PointSpreadFunction{2}` then a single argument version must be defined
```julia
intensity(psf, r::Length) = # ...
```
- `symmetry(psf) == ZAxisRadialSymmetry() && psf isa PointSpreadFunction{3}` then a two argument version must be defined
```julia
intensity(psf, r::Length, z::Length) = # ...
```
- `symmetry(psf) == NoSymmetry() && psf isa PointSpreadFunction{2}` then a two argument version must be defined
```julia
intensity(psf, x::Length, y::Length) = # ...
```
- `symmetry(psf) == NoSymmetry() && psf isa PointSpreadFunction{3}` then a three argument version must be defined
```julia
intensity(psf, x::Length, y::Length, z::Length) = # ...
```
"""
@interface intensity(psf::PointSpreadFunction{3}, x::Length, y::Length, z::Length)
@interface intensity(psf::PointSpreadFunction{3}, r::Length, z::Length)
@interface intensity(psf::PointSpreadFunction{2}, x::Length, y::Length)
@interface intensity(psf::PointSpreadFunction{2}, r::Length)

@interface Base.maximum(psf::PointSpreadFunction{2}) = response(psf, 0u"nm", 0u"nm")
@interface Base.maximum(psf::PointSpreadFunction{3}) = response(psf, 0u"nm", 0u"nm", 0u"nm")

@static if VERSION >= v"1.12.0-rc1"
    using Base: Fix
    @inline function axis_HWHM_closure(psf::PointSpreadFunction{N}, dim::Int) where {N}
        psfresponse = Fix{1}(response, psf)
        halfmax = maximum(psf) / 2
        fixed = Tuple(setdiff(1:N, dim))
        axis_response = mapfoldr(F -> Fix{F}, (a, b) -> a(b, 0u"nm"), fixed; init=psfresponse)
        return x -> axis_response(x) - halfmax
    end
else
    using Base: Fix1
    @inline function axis_HWHM_closure(psf::PointSpreadFunction{2}, dim::Int)
        psfresponse = (args...) -> response(psf, args...)
        halfmax = maximum(psf) / 2
        axis_response = dim == 1 ? x -> psfresponse(x, 0u"nm") : x -> psfresponse(0u"nm", x)
        return x -> axis_response(x) - halfmax
    end
    @inline function axis_HWHM_closure(psf::PointSpreadFunction{3}, dim::Int)
        psfresponse = (args...) -> response(psf, args...)
        halfmax = maximum(psf) / 2
        axis_response = if dim == 1
            x -> psfresponse(x, 0u"nm", 0u"nm")
        elseif dim == 2
            x -> psfresponse(0u"nm", x, 0u"nm")
        else
            x -> psfresponse(0u"nm", 0u"nm", x)
        end
        return x -> axis_response(x) - halfmax
    end
end

"""
    HWHM(psf::PointSpreadFunction)
Find the HWHM of the PSF `psf` in all the axes directions.

See also [`FWHM`](@ref)

```jldoctest
julia> tf = AiryDisc(λ=488u"nm", NA=1.4);

julia> TransferFunctions.HWHM(tf)
((-89.66947452527637 nm, 89.66947452527637 nm), (-89.66947452527637 nm, 89.66947452527637 nm), (-147.04617530370933 nm, 147.04617530370933 nm))
``` 
"""
@interface function HWHM(psf::PointSpreadFunction{N}) where {N}
    ntuple(Val(N)) do n
        hwhm_closure = axis_HWHM_closure(psf, n)
        Tuple(find_zero(hwhm_closure, lr, Bisection()) for lr in ((-Inf * u"nm", 0.0u"nm"), (0.0u"nm", Inf * u"nm")))
    end
end

"""
    FWHM(psf::PointSpreadFunction)
Find the [FWHM](https://en.wikipedia.org/wiki/Full_width_at_half_maximum) of the PSF `psf` in the x and y directions.

See also [`HWHM`](@ref)

```jldoctest
julia> tf = IsotropicGaussian(488u"nm", 1.4);

julia> TF.FWHM(tf)
(179.33894905055288 nm, 179.33894905055288 nm, 93.61250264937092 nm)
``` 
"""
@interface function FWHM(psf::PointSpreadFunction)
    map(HWHM(psf)) do (axmin, axmax)
        axmax - axmin
    end
end

# TODO: Should be a better system for handling the support on which the PSF is sampled. Instead of only having the
# energy radius or the circle on which the energy is supported, it should be returned as a general window on which the
# PSF is sampled. An elipse may also be a valid support option. The window then is the bounding box of the shape. It
# should by implemented as some trait on the PSF type that determines the dispatch of the support function. The support
# function should always hold a parameter `ε` which is the energy error that is allowed in the sampling. <22-08-25>  
@doc raw"""
    energy_radius(psf::PointSpreadFunction, ε::Real)::Length
Find the radius ``R`` such that the [encircled energy](https://en.wikipedia.org/wiki/Encircled_energy) is more than 
``1 - ε``. 

I.e. find ``R`` such that ``∫_ℬ PSF ≤ 1 - ε`` where ``ℬ`` is the ball of radius ``R``.

Useful for finding the correct PSF support to sample the PSF at for a filtering operation while guaranteeing a bounded
error of the output.
"""
@interface function energy_radius(psf::PointSpreadFunction{2}, Δ, ε::Real; oversampling=8)
    # TODO: https://chatgpt.com/share/68a77e88-a55c-8013-baf9-8e87fc4000e3 <21-08-25> 
end

"""
    psf(tf::PointSpreadFunction, Δ, wh::Dims; normalize=true)
Generate a PSF array size `wh` for the model `tf` with  the pixel size `Δ`.
"""
function psf(
    tf::PointSpreadFunction{N},
    Δ::PixelSize{N},
    inds::Indices{N}; normalize=true
) where {N}
    origin = CartesianIndex(findfirst.(==(0), inds))
    data = OriginAt(origin)(response.(tf, posgrid(inds, Δ)...))
    normalize && (data ./= sum(data))
    return SpatialArray(data, Δ)
end
psf(tf::PointSpreadFunction{N}, Δ::Length, args...; kwargs...) where {N}= psf(tf, fillsize(Δ, N), args...; kwargs...)
psf(tf::PointSpreadFunction{2}, Δ::PixelSize{2}, wh::Dims{2}; kwargs...) = psf(tf, Δ, aroundorigin(wh); kwargs...)
psf(tf::PointSpreadFunction{3}, Δ::PixelSize{3}, whd::Dims{3}; kwargs...) = psf(tf, Δ, aroundorigin(whd, (0,0,whd[3] ÷ 2 - 1)); kwargs...)

radius_window(Δ::PixelSize, R::Length) = aroundorigin(map(x -> 2 * round(Int, R / x, RoundUp) + 1, Δ))

"""
    psf(tf::PointSpreadFunction{2}, Δ, ε::Real; <kwargs>)
Sample the PSF `tf` with a pixel size `Δ` over a window such that the energy error is less than `ε`.

`kwargs` are passed to the final [`psf` function](@ref psf(::PointSpreadFunction{N}, ::PixelSize{N}, ::Indices{N}) where {N}).
"""
psf(tf::PointSpreadFunction{2}, Δ::PixelSize{2}, ε::Real; kwargs...) = psf(tf, Δ, radius_window(Δ, energy_radius(tf, ε)); kwargs...)
psf(tf::PointSpreadFunction{2}, Δ::Length, ε::Real; kwargs...) = psf(tf, fillsize(Δ, 2), ε; kwargs...)


# TODO: Specify the details `ε`, border etc. <21-08-25> 
"""
    conv(img::SpatialArray{<:Real,2}, tf::PointSpreadFunction, [border=:reflect]; <kwargs>)
Convolve the image `img` with the PSF `tf`. Additional arguments are passed to `imfilter`.
"""
function conv(img::SpatialMatrix, tf::PointSpreadFunction{2}, border=:reflect; ε=0.01)
    Δ = sampling(img)
    psf_array = psf(tf, Δ, ε)
    return conv(img, psf_array, border)
end
deconv(img::SpatialArray, tf::PointSpreadFunction) = _wiener_deconv(fft(psf(tf, sampling(img), size(img))), img.parent)

# function HWHM(psf::RadialPSF{N}) where {N}
#     ntuple(Val(N)) do n
#         hwhm_closure = axis_HWHM_closure(psf, n)
#         hwhm = find_zero(hwhm_closure, (0.0u"nm", Inf * u"nm"), Bisection())
#         (-hwhm, hwhm)
#     end
# end

"""
    params(psf::PSFModel)
Get the parameters of the PSF model `psf` as a [`ComponentVector`](@extref `ComponentArrays.ComponentVector`).

This is useful for fitting the model to data.
"""
@interface params(psf::PSFModel)::ComponentVector

"""
    fit(::PSFModel, ::SpatialArray)
"""
fit(::PSFModel, ::SpatialArray) = error("TODO")
