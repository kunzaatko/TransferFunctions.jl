abstract type PointSpreadFunction <: LinearTransferFunction end
Broadcast.broadcastable(tf::PointSpreadFunction) = Ref(tf)
intensity(psf::PointSpreadFunction, ::Length, ::Length) = throw_notimplemented_error(typeof(psf), :intensity)
psf(
    tf::PointSpreadFunction,
    Δxy::PixelSize{2},
    wh::Dims{2}
) = OriginAt(roundupcenter(wh))(intensity.(tf, posgrid(wh, Δxy)...))
psf(tf::PointSpreadFunction, Δxy::Length, wh::Dims{2}) = psf(tf, fillsize(Δxy, 2), wh)
conv(tf::PointSpreadFunction, img::SpatialArray{<:Real,2}) = imfilter(img, reflect(psf(tf, img.Δxy, size(img))))
deconv(tf::PointSpreadFunction, img::SpatialArray) = _wiener_deconv(fft(psf(tf, img.Δxy, size(img))), img)

abstract type ModelPSF <: PointSpreadFunction end
fit(::ModelPSF, ::AbstractMatrix) = error("TODO")

abstract type RadialPSF <: ModelPSF end
intensity(psf::RadialPSF, ::Length) = throw_notimplemented_error(typeof(psf), :intensity)
intensity(psf::RadialPSF, x::Length, y::Length) = intensity(psf, hypot(x, y))

abstract type MeasuredPSF <: PointSpreadFunction end

include("./psf-array.jl")
include("./gibson-lanni.jl")
include("./born-wolf.jl")
include("./estimation.jl")
