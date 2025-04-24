"""
    PointSpreadFunction <: LinearTransferFunction

# Implementation
To create a new Point spread function (PSF) `A <: PointSpreadFunction`, you must define the __intensity__ at a given [length](@extref Unitful :std:label:`Length`) coordinate `intensity(psf::A, x::Length, y::Length)`.
"""
abstract type PointSpreadFunction <: LinearTransferFunction end
Broadcast.broadcastable(tf::PointSpreadFunction) = Ref(tf)

"""
    intensity(psf::PointSpreadFunction, ::Length, ::Length)

""" # TODO: Docs <24-04-25> 
intensity(psf::PointSpreadFunction, ::Length, ::Length) = throw_notimplemented_error(typeof(psf), :intensity)

"""
    psf(tf::PointSpreadFunction, Δ::PixelSize{2}, wh::Dims{2})
""" # TODO: Docs <24-04-25> 
psf(
    tf::PointSpreadFunction,
    Δ::PixelSize{2},
    wh::Dims{2}
) = OriginAt(roundupcenter(wh))(intensity.(tf, posgrid(wh, Δ)...))
psf(tf::PointSpreadFunction, Δ::Length, wh::Dims{2}) = psf(tf, fillsize(Δ, 2), wh)

conv(tf::PointSpreadFunction, img::SpatialArray{<:Real,2}) = imfilter(img, reflect(psf(tf, sampling(img), size(img))))
deconv(tf::PointSpreadFunction, img::SpatialArray) = _wiener_deconv(fft(psf(tf, sampling(img), size(img))), img)

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
intensity(psf::RadialPSF, ::Length) = throw_notimplemented_error(typeof(psf), :intensity)
intensity(psf::RadialPSF, x::Length, y::Length) = intensity(psf, hypot(x, y))

"""
    MeasuredPSF <: PointSpreadFunction
""" # TODO: Docs <24-04-25> 
abstract type MeasuredPSF <: PointSpreadFunction end

include("./psf-array.jl")
include("./gibson-lanni.jl")
include("./born-wolf.jl")

include("./estimation.jl")

export intensity, psf
export BornWolf, GibsonLanni, PSFArray
