using InterfaceFunctions

"""
    LinearTransferFunction <: TransferFunction
""" # TODO: Docs <24-04-25> 
abstract type LinearTransferFunction <: TransferFunction end
"""
    conv(t::LinearTransferFunction, ::SpatialArray{<:Real,2})
""" # TODO: Docs <24-04-25> 
@interface conv(t::LinearTransferFunction, ::SpatialArray{<:Real,2})
"""
    deconv(t::LinearTransferFunction, ::SpatialArray{<:Real,2})
""" # TODO: Docs <24-04-25> 
@interface deconv(t::LinearTransferFunction, ::SpatialArray{<:Real,2})
transfer(t::LinearTransferFunction, img::SpatialArray{<:Real,2}) = conv(t, img)
restore(t::LinearTransferFunction, img::SpatialArray{<:Real,2}) = deconv(t, img)

function _wiener_deconv(otf_arr::AbstractArray, img::SpatialArray{<:Real,2}; snr=100.0)
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
