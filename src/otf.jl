include("models/spherical-aperture-otf.jl")

# TODO: Add docs <28-11-23> 
@doc """
optical transfer function
""" otf

# TODO: This should accept as many dimensions as the transfer function allows similar to IlluminationPatterns... Right
# now this implementation does not work<28-11-23> 

# @inline @traitfn function otf(tf::TF, freqs::Vararg{Frequency,N}) where {N,TF<:ModelOTF{N};RadiallySymmetric{TF}}

# RESEARCH: How does type conversion work in other packages that let you call a function with a specified output type?
# Could probably be used to boost the performance! <29-12-24> 
@inline @traitfn function attenuation(OT::Type{<:Number}, tf::TF, freqs::Vararg{Frequency,N}) where {N,TF<:ModelOTF{N};RadiallySymmetric{TF}}
    OT(attenuation(tf, hypot(freqs)))
end

# FIX: ↑ Above does not work as a method for `attenuation(Float64, tf, 1/200u"nm", 1/200u"nm")`. Some problem with the
# specification <28-02-25> 
@inline @traitfn function attenuation(OT::Type{<:Number}, tf::TF, fx::Frequency, fy::Frequency) where {TF <: ModelOTF{2}; RadiallySymmetric{TF}}
    OT(attenuation(tf, hypot(fx, fy)))
end

@inline @traitfn function attenuation(OT::Type{<:Number}, tf::TF, kᵣ::Frequency) where {N,TF<:ModelOTF{N};RadiallySymmetric{TF}}
    # FIX: The type should be passed to the method of on the model which for it to be able to optimize its run based on
    # it, and instead there should be a catch-all method that implements it if its not implemented in the implementation
    # <26-08-24> 
    OT(attenuation(tf, kᵣ))
end

function attenuation(tf::TF, freqs::Vararg{Frequency,N}) where {N,TF<:ModelOTF{N}}
    attenuation(preferred_type(TF), tf, freqs...)
end

# FIX: ↑ Above does not work as a method for `attenuation(tf, 1/200u"nm", 1/200u"nm")`. Some problem with the
# specification <28-02-25> 
function attenuation(tf::TF, kx::Frequency, ky::Frequency) where {TF<:ModelOTF{2}}
    attenuation(preferred_type(TF), tf, kx, ky)
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
