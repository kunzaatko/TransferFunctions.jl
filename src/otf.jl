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
@traitfn function otf_support(
    tf::TF,
    wh::Union{Tuple{Integer,Integer},Integer},
    Δxy::Union{Tuple{Length,Length},Length};
    ρ::Union{Tuple{Real,Real},Real}=1.0
) where {TF <: TransferFunction; SymmetricPupilFunction{TF}}
    wh = wh isa Integer ? (wh, wh) : wh
    Δxy = Δxy isa Length ? (Δxy, Δxy) : Δxy
    if ρ isa Real
        ρ = ρ > 0 ? (0, ρ) : (1 + ρ, 1)
    end
    fxs, fys = ndgrid(fftfreq(wh[1], 1 / Δxy[1]), fftfreq(wh[2], 1 / Δxy[2]))
    return ρ[1] * cutoff_frequency(tf) .<= hypot.(fxs, fys) .<= ρ[2] * cutoff_frequency(tf)
end
otf_support(tf::TransferFunction, img::AbstractArray, args...; varargs...) = otf_support(tf, size(img), args...; varargs...)
