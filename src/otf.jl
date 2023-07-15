@doc """
optical transfer function
""" otf

@inline @traitfn function otf(tf::TF, f_x::Frequency, f_y::Frequency) where {TF <: ClosedFormOTFModelTransferFunction; SymmetricPupilFunction{TF}}
    otf(tf, hypot(f_x, f_y))
end

function otf(tf::ClosedFormOTFModelTransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})
    fxs, fys = fftfreq(wh[1], 1 / Δxy[1]), fftfreq(wh[2], 1 / Δxy[2])
    return [otf(tf, fx, fy) for fx in fxs, fy in fys]
end

function otf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})
    tf_psf = psf(tf, wh, Δxy)
    return fft(tf_psf) ./ sum(tf_psf)
end

otf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length) = otf(tf, wh, (Δxy, Δxy))
otf(tf::TransferFunction, wh::Integer, args...) = otf(tf, (wh, wh), args...)
otf(tf::TransferFunction, img::AbstractArray, args...) = otf(tf, size(img), args...)

@doc """
 modulation transfer function
 """
mtf(tf::TransferFunction, args...; varargs...) = real.(otf(tf, args...; varargs...))

@doc """
 phase transfer function
 """
ptf(tf::TransferFunction, args...; varargs...) = imag.(otf(tf, args...; varargs...))

@traitfn function cutoff_frequency(tf::TF) where {TF <: TransferFunction; SymmetricPupilFunction{TF}}
    # TODO: Is this correct?
    all(hasfield.(TF, [:NA, :λ, :nᵢ])) ? (2 * tf.NA * tf.nᵢ) / tf.λ : nothing
end
@traitfn function otf_support(tf::TF, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length}) where {TF <: TransferFunction; SymmetricPupilFunction{TF}}
    fxs, fys = fftfreq(wh[1], 1 / Δxy[1]), fftfreq(wh[2], 1 / Δxy[2])
    return [hypot(fx, fy) < cutoff_frequency(tf) for fx in fxs, fy in fys]
end
function otf_support(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length)
    otf_support(tf, wh, (Δxy, Δxy))
end
function otf_support(tf::TransferFunction, wh::Integer, args...)
    otf_support(tf, (wh, wh), args...)
end
