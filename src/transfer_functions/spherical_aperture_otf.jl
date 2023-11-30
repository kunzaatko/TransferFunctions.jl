@doc raw"""
    IdealCircularPupilOTF(ρ₀::Real)

Ideal (aberration free) OTF of a diffraction limited imaging system with incoherent light with the cutoff-frequency `ρ₀`

The OTF is derived from the diffraction caused by the exit pupil of the system and disregards the effect of the entrance pupil... thus assumes no reshaping of the wavefronts in the optical system. The exit pupil, being located in the optical system just before the light reaches the image plane, has a greater effect on the optical system OTF.

See also `IdealCircularPupilPSF` (TODO)

# Examples
    - TODO

# Extended help

The OTF can be written in as a function of the cutoff frequency ``ρ₀`` [^1] 

```math
    ℋ(ρ) = 
    \begin{cases}
    (2/π) \left\{
        \arccos(ρ/2ρ₀) - (ρ/2ρ₀)\sqrt{1 - (ρ/2ρ₀)²}
    \right\} & \text{ for } ρ ≤ 2ρ₀ \\
        0 & \text{ otherwise}.
    \end{cases}

```
The cutoff frequency can be written in terms of the wavelength ``λ``, distance between entrance pupil and the image plane ``f₂`` and the circular pupil radius ``w`` as

```math
    \rho_0 = w/(λ f₂).
```

[^1]:
    > Frequency Analysis of Optical Imaging Systems. In Introduction to Fourier optics; Roberts & Co: Englewood, Colo, 2005; pp. 127–172 ISBN 978-0-9747077-2-3.
"""
Base.@kwdef struct IdealOTFwithCurvature{T<:Real} <: ModelOTF{2}
    λ::Length{T}
    NA::T
    nᵢ::T = 4 // 3
    curvature::T
    function IdealOTFwithCurvature(λ::Length{T}, NA::T, nᵢ::T, curvature::T) where {T<:Real}
        one(curvature) >= curvature > zero(curvature) || throw(DomainError(curvature, "Valid domain for curvature is (0,1]"))
        λ > zero(λ) || throw(DomainError(λ, "λ (wavelength) > 0"))
        NA > zero(NA) || throw(DomainError(NA, "NA (numerical aperture) > 0"))
        nᵢ > zero(nᵢ) || throw(DomainError(nᵢ, "nᵢ(refractive index of immersion medium) > 0"))
        return new{T}(λ, NA, nᵢ, curvature)
    end
end

# Casting to promoted types
function IdealOTFwithCurvature(λ::Length{R}, NA::Real, nᵢ::Real, curvature::Real) where {R<:Real}
    _, NA, nᵢ, curvature = promote(ustrip(λ), NA, nᵢ, curvature)
    return IdealOTFwithCurvature(convert(Quantity{typeof(NA)}, λ), NA, nᵢ, curvature)
end

@traitimpl RadiallySymmetric{IdealOTFwithCurvature}
output_type(::IdealOTFwithCurvature{T}) where {T} = T

function otf(tf::IdealOTFwithCurvature, fᵣ::Frequency)
    # TODO: Is the immersion refractive index necessary? <14-07-23> 
    ν = (fᵣ * tf.λ) / (2 * tf.nᵢ * tf.NA)# normalized frequency
    return ν >= one(ν) ? 0 : (2 / π) * (acos(ν) - ν * sqrt(1 - ν * ν)) * tf.curvature^ν
end
