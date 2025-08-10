using InterfaceFunctions
using Roots
using TransferFunctions.Apodization

"""
    PointSpreadFunction <: LinearTransferFunction
A point spread function is a description of a transfer functions specifying the intensity transfer of a single point
source in the object plane into a region of the image plane.

# Implementation
To create a new Point spread function (PSF) `A <: PointSpreadFunction`, you must define the __intensity__ at a given
[length](@extref Unitful `Length`) coordinate `intensity(psf::A, x::Length, y::Length)`.
"""
abstract type PointSpreadFunction <: LinearTransferFunction end
Broadcast.broadcastable(tf::PointSpreadFunction) = Ref(tf)

"""
    intensity(psf::PointSpreadFunction, x::Length, y::Length)
The intensity of a `PSF` at a given location.
"""
@interface intensity(psf::PointSpreadFunction, ::Length, ::Length)
@interface Base.maximum(psf::PointSpreadFunction) = intensity(psf, 0u"nm", 0u"nm")

"""
    FWHM(psf::PointSpreadFunction)
Find the FWHM of the PSF `psf` in the x and y directions.
```jldoctest
julia> tf = TransferFunctions.BornWolf(λ=488u"nm", NA=1.4)
BornWolf{Float64}(488.0 nm, 1.4, 1.3333333333333333)

julia> TransferFunctions.FWHM(tf)
((-119.55929936703521 nm, 119.55929936703521 nm), (-119.55929936703521 nm, 119.55929936703521 nm))
``` 
"""
@interface function FWHM(psf::PointSpreadFunction)
    max = maximum(psf)
    x_fwhm_right = find_zero(x -> intensity(psf, x * 1u"nm", 0u"nm") - max/2, (0.0, Inf)) * 1u"nm"
    x_fwhm_left = find_zero(x -> intensity(psf, x * 1u"nm", 0u"nm") - max/2, (-Inf, 0.0)) * 1u"nm"
    y_fwhm_right = find_zero(y -> intensity(psf, 0u"nm", y * 1u"nm") - max/2, (0.0, Inf)) * 1u"nm"
    y_fwhm_left = find_zero(y -> intensity(psf, 0u"nm", y * 1u"nm") - max/2, (-Inf, 0.0)) * 1u"nm"
    return ((x_fwhm_left,x_fwhm_right), (y_fwhm_left, y_fwhm_right))
end

"""
    psf(tf::PointSpreadFunction, Δ::PixelSize{2}, wh::Dims{2})
Generate a PSF array size `wh` for the model `tf` with  the pixel size `Δ`.
"""
psf(
    tf::PointSpreadFunction,
    Δ::PixelSize{2},
    wh::Dims{2}
) = OriginAt(roundupcenter(wh))(intensity.(tf, posgrid(wh, Δ)...))
psf(tf::PointSpreadFunction, Δ::Length, wh::Dims{2}) = psf(tf, fillsize(Δ, 2), wh)

"""
    conv(tf::PointSpreadFunction, img::SpatialArray{<:Real,2}, [args...]; <kwargs>)
Convolve the image `img` with the PSF `tf`. Additional arguments are passed to `imfilter`.

# Keyword arguments
- `border=nothing` If the border is a `NamedTuple` with the keys `border` of a type compatible with
[`BorderArray`](@extref) and `apodization` of type [`Apodization.ApodizationFunction`](@ref), the border is applied to the image
with the size of the FWHM of the PSF in the corresponding directions with the `border` and `apodization` used for edge
tapering ([`taperedges`](@ref))) and the full expanded array is returned.
"""
function conv(tf::PointSpreadFunction, img::SpatialArray{<:Real,2}, args...; border=nothing)
    Δ = sampling(img)
    if border == true
        border = (border=:fill, apodization=Apodization.Cosine())
    end
    if border !== nothing
        @assert border isa NamedTuple && [:border, :apodization] ⊆ keys(border) "`border` must be a NamedTuple with keys `border` and `apodization`."
        border_widths = map(FWHM(tf), sampling(img)) do fwhm, Δ
            px_widths = fwhm ./ Δ
            abs.((floor(Int, px_widths[1]),ceil(Int, px_widths[1])))
        end
        img = taperedges(border.apodization, img, border_widths, border.border) 
    end
    # FIX: I would like the SpatialArray to the be outer wrapper type <30-07-25> 
    psf_array = psf(tf, Δ, 2 .* size(img))
    psf_array ./= sum(psf_array)
    imfilter!(similar(img), img, reflect(psf_array), args...)
end
deconv(tf::PointSpreadFunction, img::SpatialArray) = _wiener_deconv(fft(psf(tf, sampling(img), size(img))), img.parent)

"""
    ModelPSF <: PointSpreadFunction 
""" # TODO: Docs <24-04-25> 
abstract type ModelPSF <: PointSpreadFunction end
"""
    fit(::ModelPSF, ::SpatialArray)
""" # TODO: Docs <24-04-25> 
fit(::ModelPSF, ::SpatialArray) = error("TODO")

"""
    RadialPSF <: ModelPSF
""" # TODO: Docs <24-04-25> 
abstract type RadialPSF <: ModelPSF end
@interface intensity(psf::RadialPSF, ::Length)
intensity(psf::RadialPSF, x::Length, y::Length) = intensity(psf, hypot(x, y))

function FWHM(psf::RadialPSF)
    max = maximum(psf)
    r_fwhm = find_zero(x -> intensity(psf, x * 1u"nm") - max/2, (0.0, Inf)) * 1u"nm"
    return ((-r_fwhm, r_fwhm), (-r_fwhm, r_fwhm))
end

"""
    MeasuredPSF <: PointSpreadFunction
""" # TODO: Docs <24-04-25> 
struct MeasuredPSF{T<:PointSpreadFunction} <: PointSpreadFunction
    psf::T
    function MeasuredPSF(psf::T) where {T<:PointSpreadFunction}
        psf isa MeasuredPSF && return psf
        return new{T}(T)
    end
end

include("./psf-array.jl")

# Models
include("./gibson-lanni.jl")
include("./born-wolf.jl")

include("./estimation.jl")

export intensity, psf
export BornWolf, GibsonLanni, PSFArray
