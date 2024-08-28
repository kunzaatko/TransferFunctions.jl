include("models/spherical_aperture_otf.jl")

# TODO: Add docs <28-11-23> 
@doc """
optical transfer function
""" otf

# FIX: This should call the type instead of it being a function on the type?? <26-08-24> 

# TODO: This should accept as many dimensions as the transfer function allows similar to IlluminationPatterns... Right
# now this implementation does not work<28-11-23> 
# @inline @traitfn function otf(tf::TF, freqs::Vararg{Frequency,N}) where {N,TF<:ModelOTF{N};RadiallySymmetric{TF}}
@inline @traitfn function otf(OT::Type{<:Number}, tf::TF, kx::Frequency, ky::Frequency) where {N,TF<:ModelOTF{N};RadiallySymmetric{TF}}
    # FIX: The type should be passed to the method of on the model which for it to be able to optimize its run based on
    # it, and instead there should be a catch-all method that implements it if its not implemented in the implementation
    # <26-08-24> 
    OT(otf(tf, hypot(kx, ky)))
end

@inline @traitfn function otf(tf::TF, kx::Frequency, ky::Frequency) where {N,TF<:ModelOTF{N};RadiallySymmetric{TF}}
    otf(preferred_type(TF), tf, kx, ky)
end

# NOTE: Has to be defined for N-dims generally and not specific dimensions because otherwise, there could be ambiguity 
# if a concrete type implements generic N-dim `otf` method <10-12-23> 
# function otf(tf::TransferFunction{N}, wh::NTuple{N,Integer}, Δxy::NTuple{N,Length}) where {N}
#     tf_psf = ifftshift(psf(tf, wh, Δxy).parent)
#     # FIX: This is not correct for a Ndim OTF only for 2D <10-12-23> 
#     # FIX!: array must be typed because of `fft` method... This should be ensured by the `psf` method and not by this 
#     # function <10-12-23> 
#     tf_psf = promote_type(unique(typeof.(tf_psf[:]))...).(tf_psf)
#     return fft(tf_psf) ./ sum(tf_psf)
# end

# TODO: Is this correct? This should be implemented as a method for a different algorithms for determining the
# cut-off frequency <24-10-23> 
@traitfn function cutoff_frequency(tf::TF) where {TF <: TransferFunction; RadiallySymmetric{TF}}
    if all(hasfield.(TF, [:NA, :λ, :nᵢ]))
        return (2 * tf.NA * tf.nᵢ) / tf.λ
    else
        throw(MethodError(cutoff_frequency, tf))
    end
end

# TODO: Move to SIM tools <24-10-23> 
# FIX: Add support for N-dim <30-11-23> 
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
