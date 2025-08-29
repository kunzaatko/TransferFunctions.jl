using InterfaceFunctions

"""
    LinearTransferFunction{N} <: TransferFunction{N}
An abstract type for all linear transfer functions.

See also [`NonLinearTransferFunction`](@ref) and [`LinearShiftInvariantTransferFunction`](@ref)
"""
abstract type LinearTransferFunction{N} <: TransferFunction{N} end

# FIX: Documentation update the mismatch of Linear and LinearShiftInvariant <12-08-25> 
"""
    LinearShiftInvariantTransferFunction{N} <: LinearTransferFunction{N}
A supertype for all linear transfer functions.

A linear transfer function is that which ensures the system to have a linear response. This means that the system obeys
the superposition principle, i.e. the response to a superposition of inputs is a superposition of the corresponding
responses, and shift-invariance principle, i.e. the response to a signal is invariant to translation, which can be
restated as "The single object in the object plane produces the same output in the image plane irrespective of its
position within the object plane."

See also [`TransferFunctions`](@ref)
"""
abstract type LinearShiftInvariantTransferFunction{N} <: LinearTransferFunction{N} end
# TODO: When it is possible to mark the non-first argument as the interfacing type, it should be done with interfaces
# again <27-08-25> 
"""
    conv(img::SpatialMatrix, ltf::LinearShiftInvariantTransferFunction)
Transfer the image `img` using the linear transfer function `ltf`, i.e. convolve the image with the equivalent PSF.
"""
function conv(::SpatialMatrix, t::LinearShiftInvariantTransferFunction) end

# TODO: When it is possible to mark the non-first argument as the interfacing type, it should be done with interfaces
# again <27-08-25> 
"""
    deconv(img::SpatialArray, ltf::LinearShiftInvariantTransferFunction)
Deconvolve the image `img` that was transferred using the linear transfer function `ltf` using the specified algorithm.
"""
function deconv(::SpatialMatrix, t::LinearShiftInvariantTransferFunction) end
transfer(img::SpatialMatrix, t::LinearShiftInvariantTransferFunction) = conv(img, t)
restore(img::SpatialMatrix, t::LinearShiftInvariantTransferFunction) = deconv(img, t)

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

# TODO: Update the documentation. This is in-fact linear <15-08-25> 
# TODO: Create a structure for this type branch <01-08-25> 
"""
    ImpulseResponseMapping{2} <: LinearTransferFunction{2}
An impulse response mapping is a description of a transfer function that prescribes an impulse response to every point
in the object plane. This may be a [`point spread function`](@ref PointSpreadFunction) or some other function that
prescribes the response to the object plane point in the image plane.

An transfer defined by an impulse response mapping is a linear system, but is not translation invariant hence does not
belong to the class of linear transfer functions.
"""
struct ImpulseResponseMapping <: LinearTransferFunction{2}
end

include("otf/optical-transfer-function.jl")
include("psf/point-spread-function.jl")

export conv, deconv
