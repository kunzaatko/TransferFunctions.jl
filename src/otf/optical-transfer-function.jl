"""
    OpticalTransferFunction <: LinearTransferFunction

# Implementation
To create a new OTF model `A <: OpticalTransferFunction`, you must define the attenuation at a given [frequency](@ref TransferFunctions.Frequency) coordinate `attenuation(model::A, kx::Frequency, ky::Frequency)`.
"""
abstract type OpticalTransferFunction <: LinearTransferFunction end
Broadcast.broadcastable(tf::OpticalTransferFunction) = Ref(tf)
conv(::OpticalTransferFunction, ::SpatialImage) = error("TODO")
deconv(::OpticalTransferFunction, ::SpatialImage) = error("TODO")
attenuation(otf::OpticalTransferFunction, ::Frequency, ::Frequency) = no_implemementation_error(typeof(otf), :attenuation)
insupport(otf::OpticalTransferFunction, fx::Frequency, fy::Frequency) = error("TODO")


"""
    RadialOTF <: OpticalTransferFunction
If the pupil function of the system is symmetric, the OTF as well as the PSF are radially symmetric which can be used to optimize the calculations.

# Implementation
- `attenuation(model::A, fr::Frequency)` 
- `cutoff(model::A, [a=0])` returning the largest frequency `f` such that `attenuation(model, f) >= a` if `a > 0` and
`attenuation(model, f) > 0` if `a = 0`.
"""
abstract type RadialOTF <: OpticalTransferFunction end
attenuation(otf::RadialOTF, ::Frequency) = no_implemementation_error(typeof(otf), :attenuation)
attenuation(otf::RadialOTF, fx::Frequency, fy::Frequency) = attenuation(otf, hypot(fx, fy))
cutoff(::RadialOTF, a=0.0) = error("TODO")
insupport(otf::RadialOTF, fx::Frequency, fy::Frequency) = hypot(fx, fy) <= cutoff(otf)


otf(
    tf::OpticalTransferFunction,
    Δxy::PixelSize{2},
    wh::Dims{2}
) = attenuation.(tf, fftfreqs(wh, Δxy)...)
otf(tf::OpticalTransferFunction, Δxy::Length, wh::Dims{2}) = otf(tf, fillsize(Δxy, 2), wh)
otf(tf::OpticalTransferFunction, img::SpatialImage{T,2}) where {T} = otf(tf, img.Δxy, size(img))

include("./otf-array.jl")
include("./circular-pupil-otf.jl")
