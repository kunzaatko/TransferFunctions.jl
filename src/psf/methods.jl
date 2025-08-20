using Roots, InterfaceFunctions
using TransferFunctions.Apodization

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
@interface intensity(psf::PointSpreadFunction, args...)

@interface Base.maximum(psf::PointSpreadFunction{2}) = response(psf, 0u"nm", 0u"nm")
@interface Base.maximum(psf::PointSpreadFunction{3}) = response(psf, 0u"nm", 0u"nm", 0u"nm")

@static if VERSION >= v"1.12"
    using Base: Fix
    @inline function axis_HWHM_closure(psf::PointSpreadFunction{N}, dim::Int) where {N}
        psfresponse = Fix{1}(response, psf)
        halfmax = maximum(psf) / 2
        fixed = Tuple(setdiff(1:N, dim))
        axis_response = mapfoldr(F -> Fix{F}, (a,b) -> a(b,0u"nm"), fixed; init = psfresponse)
        return x -> axis_response(x) - halfmax
    end
else
    using Base: Fix1
    @inline function axis_HWHM_closure(psf::PointSpreadFunction{2}, dim::Int) where {N}
        psfresponse = Fix1(response, psf)
        halfmax = maximum(psf) / 2
        axis_response = dim == 1 ? x -> psfresponse(x, 0u"nm") : x -> psfresponse(0u"nm", x)
        return x -> axis_response(x) - halfmax
    end
    @inline function axis_HWHM_closure(psf::PointSpreadFunction{3}, dim::Int) where {N}
        psfresponse = Fix1(response, psf)
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
((-89.66947452527637 nm, 89.66947452527637 nm), (-89.66947452527637 nm, 89.66947452527637 nm), (-46.80625132468546 nm, 46.80625132468546 nm))
``` 
"""
@interface function HWHM(psf::PointSpreadFunction{N}) where {N}
    ntuple(Val(N)) do n
        hwhm_closure = axis_HWHM_closure(psf, n)
        Tuple(find_zero(hwhm_closure, lr, Bisection()) for lr in ((-Inf * u"nm",  0.0u"nm"), (0.0u"nm", Inf * u"nm")))
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


"""
    psf(tf::PointSpreadFunction, Δ::PixelSize{2}, wh::Dims{2}; normalize=true)
Generate a PSF array size `wh` for the model `tf` with  the pixel size `Δ`.
"""
function psf(
    tf::PointSpreadFunction,
    Δ::PixelSize{2},
    wh::Dims{2}; normalize=true
)
    data = OriginAt(roundupcenter(wh))(response.(tf, posgrid(wh, Δ)...))
    normalize && (data ./= sum(data))
    return SpatialMatrix(data, Δ)
end
psf(tf::PointSpreadFunction, Δ::Length, wh::Dims{2}; kwargs...) = psf(tf, fillsize(Δ, 2), wh; kwargs...)

"""
    psf(tf::PointSpreadFunction, Δ::PixelSize{3}, whd::Dims{3}; normalize=true)
Generate a 3D PSF array of size `whd` for the model `tf` with  the voxel size `Δ`.
"""
function psf(
    tf::PointSpreadFunction,
    Δ::PixelSize{3},
    whd::Dims{3}; normalize=true
)
    center = CartesianIndex(Tuple(roundupcenter(whd[1:2]))..., 1)
    grid = posgrid(whd, Δ; center)
    data = OriginAt(center)(response.(tf, grid...))
    normalize && (data ./= sum(data))
    return SpatialArray(data, Δ)
end
psf(tf::PointSpreadFunction, Δ::Length, wh::Dims{3}; kwargs...) = psf(tf, fillsize(Δ, 3), wh; kwargs...)

"""
    conv(tf::PointSpreadFunction, img::SpatialArray{<:Real,2}, [args...]; <kwargs>)
Convolve the image `img` with the PSF `tf`. Additional arguments are passed to `imfilter`.

# Keyword arguments
- `border=nothing` If the border is a `NamedTuple` with the keys `border` of a type compatible with
[`BorderArray`](@extref) and `apodization` of type [`Apodization.ApodizationFunction`](@ref), the border is applied to the image
with the size of the FWHM of the PSF in the corresponding directions with the `border` and `apodization` used for edge
tapering ([`taperedges`](@ref))) and the full expanded array is returned.
"""
function conv(tf::PointSpreadFunction, img::SpatialMatrix{<:Real}, args...; border=nothing)
    Δ = sampling(img)
    if border == true
        border = (border=:fill, apodization=Apodization.Cosine())
    end
    if border !== nothing
        @assert border isa NamedTuple && [:border, :apodization] ⊆ keys(border) "`border` must be a NamedTuple with keys `border` and `apodization`."
        border_widths = map(FWHM(tf), sampling(img)) do fwhm, Δ
            px_widths = fwhm ./ Δ
            abs.((floor(Int, px_widths[1]), ceil(Int, px_widths[1])))
        end
        img = taperedges(border.apodization, img, border_widths, border.border)
    end
    # FIX: I would like the SpatialArray to the be outer wrapper type <30-07-25> 
    psf_array = psf(tf, Δ, 2 .* size(img))
    psf_array ./= sum(psf_array)
    return conv(img, psf_array)
end
deconv(tf::PointSpreadFunction, img::SpatialArray) = _wiener_deconv(fft(psf(tf, sampling(img), size(img))), img.parent)

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
@interface params(psf::PSFModel)

"""
    fit(::PSFModel, ::SpatialArray)
"""
fit(::PSFModel, ::SpatialArray) = error("TODO")
