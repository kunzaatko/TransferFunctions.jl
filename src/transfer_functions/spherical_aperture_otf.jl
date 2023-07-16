Base.@kwdef struct IdealOTFwithCurvature{T<:Real} <: ClosedFormOTFModelTransferFunction
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

@traitimpl SymmetricPupilFunction{IdealOTFwithCurvature}
output_type(::IdealOTFwithCurvature{T}) where {T} = T

function otf(tf::IdealOTFwithCurvature, fᵣ::Frequency)
    # TODO: Is the immersion refractive index necessary? <14-07-23> 
    ν = (fᵣ * tf.λ) / (2 * tf.nᵢ * tf.NA)# normalized frequency
    return ν >= one(ν) ? 0 : (2 / π) * (acos(ν) - ν * sqrt(1 - ν * ν)) * tf.curvature^ν
end
