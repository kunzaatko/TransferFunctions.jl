"""
    CircularPupilOTF(ρ₀::Real)

Ideal (aberration free) OTF of a diffraction limited imaging system with incoherent light with the cutoff-frequency `ρ₀`
"""
Base.@kwdef struct CircularPupilOTF{T<:Real} <: RadialOTF
    λ::Length{T}
    NA::T
    nᵢ::T = 4 // 3
    curvature::T
    function CircularPupilOTF(λ::Length{T}, NA::T, nᵢ::T, curvature::T) where {T<:Real}
        one(curvature) >= curvature > zero(curvature) || throw(DomainError(curvature, "Valid domain for curvature is (0,1]"))
        λ > zero(λ) || throw(DomainError(λ, "λ (wavelength) > 0"))
        NA > zero(NA) || throw(DomainError(NA, "NA (numerical aperture) > 0"))
        nᵢ > zero(nᵢ) || throw(DomainError(nᵢ, "nᵢ(refractive index of immersion medium) > 0"))
        return new{T}(λ, NA, nᵢ, curvature)
    end
end

# Casting to promoted types
function CircularPupilOTF(λ::Length{R}, NA::Real, nᵢ::Real, curvature::Real) where {R<:Real}
    _, NA, nᵢ, curvature = promote(ustrip(λ), NA, nᵢ, curvature)
    return CircularPupilOTF(convert(Quantity{typeof(NA)}, λ), NA, nᵢ, curvature)
end

attenuation_normalized(tf::CircularPupilOTF, ν::Number) = ν >= one(ν) ? 0 : (2 / π) * (acos(ν) - ν * sqrt(1 - ν * ν)) * tf.curvature^ν
# TODO: Is the immersion refractive index necessary? <14-07-23> 
attenuation(tf::CircularPupilOTF, fᵣ::Frequency) = attenuation_normalized(tf, (fᵣ * tf.λ) / (2 * tf.nᵢ * tf.NA))

function cutoff(tf::CircularPupilOTF{T}, a::Real)::Frequency where {T}
    if iszero(a)
        return (2 * tf.NA * tf.nᵢ) / tf.λ
    else
        0 <= a <= 1 || throw(DomainError(a, "OTF can only attenuate with a coefficient between 0 and 1"))
        z = find_zero(ν -> attenuation_normalized(tf, ν) - a, (0, 1), Bisection())
        return (2z * tf.nᵢ * tf.NA) / tf.λ
    end
end
