abstract type PointSpreadFunction <: LinearTransferFunction end
Broadcast.broadcastable(tf::PointSpreadFunction) = Ref(tf)
intensity(psf::PointSpreadFunction, ::Length, ::Length) = no_implemementation_error(typeof(psf), :intensity)
function psf(
    tf::PointSpreadFunction,
    Δxy::PixelSize{2},
    wh::Dims{2}
)
    xs = isodd(wh[1]) ? (((-wh[1]-1)÷2):((wh[1]-1)÷2)) .* Δxy[1] : (((-wh[1]-2)÷2):(wh[1]÷2)) .* Δxy[1]
    ys = isodd(wh[2]) ? (((-wh[2]-1)÷2):((wh[2]-1)÷2)) .* Δxy[2] : (((-wh[2]-2)÷2):(wh[2]÷2)) .* Δxy[2]
    return centered([intensity(tf, x, y) for x in xs, y in ys])
end
psf(tf::PointSpreadFunction, Δxy::Length, wh::Dims{2}) = psf(tf, fillsize(Δxy, 2), wh)
conv(tf::PointSpreadFunction, img::SpatialImage{T,2}) where {T} = imfilter(img, psf(tf, img))
deconv(::PointSpreadFunction, ::SpatialImage) = error("TODO")

abstract type ModelPSF <: PointSpreadFunction end
fit(::ModelPSF, ::AbstractMatrix) = error("TODO")

abstract type RadialPSF <: ModelPSF end
intensity(psf::RadialPSF, ::Length) = no_implemementation_error(typeof(psf), :intensity)
intensity(psf::RadialPSF, x::Length, y::Length) = intensity(psf, hypot(x, y))

abstract type MeasuredPSF <: PointSpreadFunction end

include("./psf-array.jl")
include("./gibson-lanni.jl")
include("./born-wolf.jl")
include("./estimation.jl")
