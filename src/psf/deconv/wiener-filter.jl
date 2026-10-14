# TODO: Add abstract field to allow to estimate the frequency spectral density and noise spectral density instead of
# simple SNR constant <01-09-25> 
# TODO: Documentation <01-09-25> 
"""
    WienerFilter{T<:Real}
Wiener deconvolution parameters for deconvolving images using the Wiener filter.
"""
struct WienerFilter{T<:Real}
    snr::T
end

# FIX: The filter function should only get the indices of A and generate a filter for reuse <03-09-25> 
function filter(f::WienerFilter, A, tf::PointSpreadFunction, border=:reflect; ε=0.05)
    K = psf(tf, sampling(A), ε; normalize=true)
    A_padded = border_array(A, border, kern_padding(K))
    Kv = FFTView(zeros(eltype(K), map(length, axes(A_padded))))
    for I in CartesianIndices(axes(K))
        Kv[I] = K[I]
    end

    K_fft = FFT.fft(Kv) # FIX: Should `real` be here? Does it make sense to use for PSFs that have a real output? <03-09-25> 

    flt = conj.(K_fft) ./ (abs2.(K_fft) .+ (one(eltype(K_fft)) / f.snr))
    return flt[axes(A)...]
end

# TODO: There is a need for edge-tapering here. I should add a parameter that will do that <03-09-25> 
"""
    deconv(img::SpatialArray, tf::PointSpreadFunction, wiener::WienerFilter, [border])
Deconvolve the image `img` using the point spread function `tf` and Wiener deconvolution parameters `wiener`.

The Wiener filter is applied in the frequency domain using the formula:
`H*(|H|²) / (|H|² + 1/SNR)` where `H` is the PSF transfer function and SNR is from `wiener.snr`.
"""
function deconv(A::SpatialArray, tf::PointSpreadFunction, wiener::WienerFilter, border=:reflect; ε=0.05)
    Δ = sampling(A)
    F = filter(wiener, A, tf, border; ε)
    A_fft = FFT.fft(A) # TODO: Apply the border before the filtering <03-09-25> 
    A_filtered = real(FFT.ifft(A_fft .* F)) # FIX: Is real a good idea here? <03-09-25> 
    return SpatialArray(A_filtered, Δ)
end

export WienerFilter
