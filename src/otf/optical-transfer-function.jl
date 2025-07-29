using InterfaceFunctions

"""
    OpticalTransferFunction <: LinearTransferFunction

# Implementation
To create a new Optical transfer function (OTF) `A <: OpticalTransferFunction`, you must define the __attenuation__ at a given [frequency](@ref TransferFunctions.Frequency) coordinate `attenuation(otf::A, kx::Frequency, ky::Frequency)`.
"""
abstract type OpticalTransferFunction <: LinearTransferFunction end
Broadcast.broadcastable(tf::OpticalTransferFunction) = Ref(tf)
_fft_conv(otf_arr::AbstractArray, img::AbstractArray) = ifft(otf_arr .* fft(img))
conv(tf::OpticalTransferFunction, img::SpatialArray{<:Real,2}) = _fft_conv(otf(tf, img), img)
conv(tf::OpticalTransferFunction, img::SpatialArray{<:Gray,2}) = conv(tf, channelview(img))
deconv(tf::OpticalTransferFunction, img::SpatialArray{<:Real,2}) = _wiener_deconv(otf(tf, img), img)
"""
    attenuation(otf::OpticalTransferFunction, ::Frequency, ::Frequency)
""" # TODO: Docs <24-04-25> 
@interface attenuation(otf::OpticalTransferFunction, ::Frequency, ::Frequency)
"""
    insupport(otf::OpticalTransferFunction, fx::Frequency, fy::Frequency)
""" # TODO: Docs <24-04-25> 
function insupport(otf::OpticalTransferFunction, fx::Frequency, fy::Frequency)
    a = abs(attenuation(otf, fx, fy))
    return a > zero(a)
end


"""
    otf(tf::OpticalTransferFunction, Δxy::PixelSize{2}, wh::Dims{2})
""" # TODO: Docs <24-04-25> 
otf(
    tf::OpticalTransferFunction,
    Δxy::PixelSize{2},
    wh::Dims{2}
) = attenuation.(tf, fftfreqs(wh, Δxy)...)
otf(tf::OpticalTransferFunction, Δ::Length, wh::Dims{2}) = otf(tf, fillsize(Δ, 2), wh)
otf(tf::OpticalTransferFunction, img::SpatialArray{T,2}) where {T} = otf(tf, sampling(img), size(img))

"""
    RadialOTF <: OpticalTransferFunction
If the pupil function of the system is symmetric, the OTF as well as the PSF are radially symmetric which can be used to optimize the calculations.

# Implementation
- `attenuation(model::A, fr::Frequency)` 
- `cutoff(model::A, [a=0])` returning the largest frequency `f` such that `attenuation(model, f) >= a` if `a > 0` and
`attenuation(model, f) > 0` if `a = 0`.
"""
abstract type RadialOTF <: OpticalTransferFunction end
@interface attenuation(otf::RadialOTF, ::Frequency)
attenuation(otf::RadialOTF, fx::Frequency, fy::Frequency) = (@inline; attenuation(otf, hypot(fx, fy)))
"""
    cutoff(::RadialOTF, a=0.0)
""" # TODO: Docs <24-04-25> 
cutoff(::RadialOTF, a=0.0) = error("TODO")
insupport(otf::RadialOTF, fr::Frequency) = (@inline; fr <= cutoff(otf))
insupport(otf::RadialOTF, fx::Frequency, fy::Frequency) = (@inline; hypot(fx, fy) <= cutoff(otf))


include("./otf-array.jl")
include("./circular-pupil-otf.jl")

export otf, cutoff, attenuation
export CircularPupilOTF, OTFArray
