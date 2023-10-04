@doc raw"""
Born & Wolf model of the transfer function for a circular aperture.

# Parameters
  + `λ::Unitful.Length`: (emission) wavelength 
  + `nᵢ::Number`: index of refraction of the immersion medium
  + `NA::Number`: numerical aperture

# Extended help

The Born & Wolf model is a scalar diffraction model derived for perfect systems. It assumes that the only aberration of the system is due to *defocus*. Modern microscope objectives are designed to provide optimal imaging conditions for sources located directly on the coverslip, in which case the Born & Wolf model is applicable (if the coverslip and immersion is used as designed). The model disregards spherical and higher order aberrations that are due to the source of illumination being shifted from the coverslip boundary.
"""
Base.@kwdef struct BornWolf{T<:Real} <: ClosedFormPSFModel
    λ::Length{T}
    NA::T
    nᵢ::T = 4 // 3
    function BornWolf(λ::Length{T}, NA::T, nᵢ::T) where {T<:Real}
        λ > zero(λ) || throw(DomainError(λ, "λ (wavelength) > 0"))
        NA > zero(NA) || throw(DomainError("NA (numerical aperture) > 0"))
        nᵢ > zero(nᵢ) || throw(DomainError("nᵢ(refractive index of immersion medium) > 0"))
        return new{T}(λ, NA, nᵢ)
    end
end

# Casting to promoted types
function BornWolf(λ::Length{R}, NA::Number, nᵢ::Number) where {R<:Real}
    _, NA, nᵢ = promote(ustrip(λ), NA, nᵢ)
    return BornWolf(convert(Quantity{typeof(NA)}, λ), NA, nᵢ)
end
@traitimpl RadiallySymmetric{BornWolf}
output_type(::BornWolf{T}) where {T} = T

function psf(tf::BornWolf{T}, r::Length)::T where {T}
    k = 2π / tf.λ
    k₀ = k / tf.nᵢ
    # FIX: Is this correct?! <14-07-23> 
    return r == zero(r) ? 1 : (2besselj1(k₀ * r * tf.NA) / (k₀ * r * tf.NA))^2
end
