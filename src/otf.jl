@doc """
optical transfer function
""" otf

@inline @traitfn function otf(tf::TF, f_x::Frequency, f_y::Frequency) where {TF <: ClosedFormOTFModel; RadiallySymmetric{TF}}
    otf(tf, hypot(f_x, f_y))
end

@doc """
    otf(tf, wh, [Δxy]; δ=(0, 0))
    otf(tf, img, [Δxy];...)


Generate an otf for the given transfer function with the size `wh` (size of `img`) and with a pixel distance of `Δxy`

!!! note
    The pixel size/distance (`Δxy`) is only required for a [model transfer function](@ref ModelTransferFunction). If it
    is  not provided for a `MeasuredTransferFunction`, then the pixel distance from the measurement is used.

# Arguments
* `tf::TransferFunction`: transfer function model/measure to generate the OTF for
* `wh::Tuple{Integer, Integer}` or `wh::Integer`: (width, height) of the generated OTF. `wh` ↦ `(wh, wh)` if `wh isa Integer`.
* `Δxy::Tuple{Length, Length}` or `Δxy::Length`: Separation of pixels in the ``x`` and ``y`` dimensions of the generated
    OTF image. `Δxy` ↦ `(Δxy, Δxy)` if `Δxy isa Length`.
* `δ::Tuple = (0,0)`: shift of the OTF in the image plane in pixels. This is useful for some algorithms, e.g. in 
    structured illumination microscopy reconstruction algorithms.
"""
# TODO: This should be shifted for a PSF that is real <15-09-23> 
function otf(
    tf::ClosedFormOTFModel,
    wh::Tuple{Integer,Integer},
    Δxy::Tuple{Length,Length};
    δ::Tuple{<:Real,<:Real}=(0, 0)
)
    fxs, fys = ndgrid(fftfreq(wh[1], 1 / Δxy[1]), fftfreq(wh[2], 1 / Δxy[2]))
    # PERF: Could be speeded up using `ShiftedArrays` with integer `δ` <20-07-23> 
    return otf.(tf,
        fxs .- δ[1] / (Δxy[1] * wh[1]),
        fys .- δ[2] / (Δxy[2] * wh[2]))
end

function otf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})
    tf_psf = psf(tf, wh, Δxy)
    return fft(tf_psf) ./ sum(tf_psf)
end

otf(tf::TransferFunction, wh::Integer, args...; varargs...) = otf(tf, (wh, wh), args...; varargs...)
otf(tf::TransferFunction, wh::Tuple, Δxy::Length; varargs...) = otf(tf, wh, (Δxy, Δxy); varargs...)
otf(tf::TransferFunction, img::AbstractMatrix, args...; varargs...) = otf(tf, size(img), args...; varargs...)

@doc """
 modulation transfer function
 """
mtf(args...; varargs...) = real.(otf(args...; varargs...))

@doc """
 phase transfer function
 """
ptf(args...; varargs...) = imag.(otf(args...; varargs...))

@traitfn function cutoff_frequency(tf::TF) where {TF <: TransferFunction; RadiallySymmetric{TF}}
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
    ρ::Union{Real,Tuple{Real,Real}}=(0.0, 1.0),
    inclusive::Union{Bool,Tuple{Bool,Bool}}=(true, true)
) where {TF <: TransferFunction; RadiallySymmetric{TF}}
    if ρ isa Real
        ρ = ρ > 0 ? (0, ρ) : (1 + ρ, 1)
    end
    inclusive = inclusive isa Bool ? (inclusive, inclusive) : inclusive

    fxs, fys = ndgrid(fftfreq(wh[1], 1 / Δxy[1]), fftfreq(wh[2], 1 / Δxy[2]))
    left = inclusive[1] ? ρ[1] * cutoff_frequency(tf) .<= hypot.(fxs, fys) : ρ[1] * cutoff_frequency(tf) .< hypot.(fxs, fys)
    right = inclusive[2] ? hypot.(fxs, fys) .< ρ[2] * cutoff_frequency(tf) : hypot.(fxs, fys) .< ρ[2] * cutoff_frequency(tf)
    return left .* right
end
otf_support(tf, wh::Integer, args...; varargs...) = otf_support(tf, (wh, wh), args...; varargs...)
otf_support(tf, wh::Tuple, Δxy::Length, args...; varargs...) = otf_support(tf, wh, (Δxy, Δxy), args...; varargs...)
otf_support(tf, img::AbstractArray, args...; varargs...) = otf_support(tf, size(img), args...; varargs...)
