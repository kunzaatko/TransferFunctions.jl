using InterfaceFunctions

abstract type BoundaryCondition end

"""
    LinearTransferFunction <: TransferFunction
A supertype for all linear transfer functions.

A linear transfer function is that which ensures the system to have a linear response. This means that the system obeys
the superposition principle, i.e. the response to a superposition of inputs is a superposition of the corresponding
responses, and shift-invariance principle, i.e. the response to a signal is invariant to translation, which can be
restated as "The single object in the object plane produces the same output in the image plane irrespective of its
position within the object plane."

See also [`TransferFunctions`](@ref)
"""
abstract type LinearTransferFunction <: TransferFunction end
"""
    conv(ltf::LinearTransferFunction, img::SpatialArray{<:Real,2})
Transfer the image `img` using the linear transfer function `ltf`, i.e. convolve the image with the equivalent PSF.
"""
@interface conv(t::LinearTransferFunction, ::SpatialMatrix{<:Real})

"""
    deconv(ltf::LinearTransferFunction, img::SpatialArray{<:Real,2})
Deconvolve the image `img` that was transferred using the linear transfer function `ltf` using the specified algorithm.
"""
@interface deconv(t::LinearTransferFunction, ::SpatialMatrix{<:Real})
transfer(t::LinearTransferFunction, img::SpatialMatrix{<:Real}) = conv(t, img)
restore(t::LinearTransferFunction, img::SpatialMatrix{<:Real}) = deconv(t, img)

function _wiener_deconv(otf_arr::AbstractArray, img::SpatialMatrix{<:Real}; snr=100.0)
    # Convert blurred image to frequency domain
    F = fft(img)

    # Calculate Wiener filter components
    H_conj = conj(otf_arr)
    H_sq_abs = abs2.(otf_arr)
    noise_power = 1.0 / snr

    # Apply Wiener filter
    G = (H_conj ./ (H_sq_abs .+ noise_power)) .* F

    # Return real part of inverse transform
    real(ifft(G))
end

include("otf/optical-transfer-function.jl")
include("psf/point-spread-function.jl")

export conv, deconv
