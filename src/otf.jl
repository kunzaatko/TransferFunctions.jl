@doc """
optical transfer function
""" otf

@inline @traitfn function otf(tf::TF, f_x::Frequency, f_y::Frequency) where {TF <: ClosedFormOTFModelTransferFunction; SymmetricPupilFunction{TF}}
    otf(tf, hypot(f_x, f_y))
end

function otf(tf::ClosedFormOTFModelTransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})
    fxs, fys = ndgrid(fftfreq(wh[1], 1 / Δxy[1]), fftfreq(wh[2], 1 / Δxy[2]))
    return otf.(tf, fxs, fys)
end

function otf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})
    tf_psf = psf(tf, wh, Δxy)
    return fft(tf_psf) ./ sum(tf_psf)
end

otf(tf::TransferFunction, wh::Integer, args...; varargs...) = otf(tf, (wh, wh), args...; varargs...)
otf(tf::TransferFunction, wh::Tuple, Δxy::Length; varargs...) = otf(tf, wh, (Δxy, Δxy); varargs...)
otf(tf::TransferFunction, img::AbstractArray, args...; varargs...) = otf(tf, size(img), args...; varargs...)

@doc """
 modulation transfer function
 """
mtf(args...; varargs...) = real.(otf(args...; varargs...))

@doc """
 phase transfer function
 """
ptf(args...; varargs...) = imag.(otf(args...; varargs...))

@traitfn function cutoff_frequency(tf::TF) where {TF <: TransferFunction; SymmetricPupilFunction{TF}}
    if all(hasfield.(TF, [:NA, :λ, :nᵢ]))
        return (2 * tf.NA * tf.nᵢ) / tf.λ
    else
        throw(MethodError(cutoff_frequency, tf))
    end
end

@traitfn function otf_support(
    tf::TF,
    wh::Tuple{Integer,Integer},
    Δxy::Tuple{Length,Length};
    ρ::Union{Real,Tuple{Real,Real}}=(0.0, 1.0)
) where {TF <: TransferFunction; SymmetricPupilFunction{TF}}
    if ρ isa Real
        ρ = ρ > 0 ? (0, ρ) : (1 + ρ, 1)
    end
    fxs, fys = ndgrid(fftfreq(wh[1], 1 / Δxy[1]), fftfreq(wh[2], 1 / Δxy[2]))
    return ρ[1] * cutoff_frequency(tf) .<= hypot.(fxs, fys) .<= ρ[2] * cutoff_frequency(tf)
end
otf_support(tf, wh::Integer, args...; varargs...) = otf_support(tf, (wh, wh), args...; varargs...)
otf_support(tf, wh::Tuple, Δxy::Length, args...; varargs...) = otf_support(tf, wh, (Δxy, Δxy), args...; varargs...)
otf_support(tf, img::AbstractArray, args...; varargs...) = otf_support(tf, size(img), args...; varargs...)
