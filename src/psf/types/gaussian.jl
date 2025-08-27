using Roots, SpecialFunctions, ComponentArrays

@doc raw"""
    const C_AiryDisc_lateral
``c`` such that the FWHM of the Airy disc in the lateral direction is ``c * λ/NA``.
"""
const C_AiryDisc_lateral = find_zero(x -> (2besselj1(x * π) / (x * π))^2 - 1 / 2, (0.4, 0.6))

@doc raw"""
    const C_AiryDisc_axial
``c`` such that the FWHM of the Airy disc in the axial direction is ``c * λn/NA^2``.
"""
const C_AiryDisc_axial = find_zero(x -> sinc(x * π)^2 - 1 / 2, (0.1, 0.2))

"""
    IsotropicGaussian{T<:Real} <: PSFModel{3}
Airy disc scalar diffraction point spread function model for a circular aperture.

# Fields
- `λ::Length`: emission wavelength
- `NA::Real`: numerical aperture of the objective
- `n::Real`: refractive index of the immersion medium (defaults to `4//3` which is the refractive index of water)
"""
struct IsotropicGaussian{N,T<:Real} <: PSFModel{N}
    λ::Length{T}
    NA::T
    n::T # only relevant for 3D
    C_lateral::T
    C_axial::T # only relevant for 3D
    function IsotropicGaussian{N,T}(λ::Length{T}, NA::T, n::T, C_lateral::T, C_axial::T) where {T,N}
        N ∈ (2, 3) || throw(ArgumentError("`IsotropicGaussian` model can only be used for 2D and 3D PSFs. Got `N == $N`."))
        check_emission_wavelength(λ)
        check_numerical_aperture(NA)
        check_refractive_index(n)
        C_lateral > zero(C_lateral) || throw(DomainError(C_lateral, "Lateral FWHM coefficient must be a positive value. Got `C_lateral = $C_lateral`."))
        C_axial > zero(C_axial) || throw(DomainError(C_axial, "Axial FWHM coefficient must be a positive value. Got `C_axial = $C_axial`."))
        return new{N,T}(λ, NA, n, C_lateral, C_axial)
    end
end
IsotropicGaussian{N,T}(; λ::Length{T}, NA::T, n::T=T(4 // 3), C_lateral::T=T(C_AiryDisc_lateral), C_axial::T=T(C_AiryDisc_axial)) where {N,T} = IsotropicGaussian{N,T}(λ, NA, n, C_lateral, C_axial) # final 1
IsotropicGaussian{2,T}(params::ComponentVector) where {T} = IsotropicGaussian{2,T}(;λ=params.λ, NA=T(params.NA), C_lateral=T(params.C_lateral)) # -> final 1
IsotropicGaussian{3,T}(params::ComponentVector) where {T} = IsotropicGaussian{3,T}(;λ=params.λ, NA=T(params.NA), n=T(params.n), C_lateral=T(params.C_lateral), C_axial=T(params.C_axial)) # -> final 1

# TODO: Document the FWHM coefficients <18-08-25> 
"""
    IsotropicGaussian(λ::Length, NA, n=4//3)
[`IsotropicGaussian`](@ref) model point spread function with emission wavelength `λ`, numerical aperture `NA` and
immersion medium refractive index `n`.
"""
function IsotropicGaussian{N}(λ::Length{A}, NA::Real, n::Real, C_lateral::Real, C_axial::Real) where {A<:Real,N} # final 2
    T = promote_type(A, typeof(NA), typeof(n), typeof(C_lateral), typeof(C_axial))
    λ = convert(T, ustrip(λ)) * unit(λ)
    NA, n, C_lateral, C_axial = convert(T, NA), convert(T, n), convert(T, C_lateral), convert(T, C_axial)
    return IsotropicGaussian{N,T}(λ, NA, n, C_lateral, C_axial) # inner
end
IsotropicGaussian{N}(; λ, NA, n=4 // 3, C_lateral=C_AiryDisc_lateral, C_axial=C_AiryDisc_axial) where {N} = IsotropicGaussian{N}(λ, NA, n, C_lateral, C_axial) # -> final 2
IsotropicGaussian{N}(λ::Length, NA::Real; kwargs...) where {N} = IsotropicGaussian{N}(; λ, NA, kwargs...) # -> final 2
IsotropicGaussian(args...; kwargs...) = IsotropicGaussian{3}(args...; kwargs...) # -> final 2


symmetry(::IsotropicGaussian) = ZAxisRadialSymmetry()

@inline σ_xy(tf::IsotropicGaussian) = tf.C_lateral * (tf.λ / tf.NA) / (2 * √(2log(2)))
function intensity(tf::IsotropicGaussian, r::Length)
    σ = σ_xy(tf)
    w_xy = ustrip(2π * σ^2) # NOTE: dimension gets integrated out from the normalization factor
    return exp(-r^2 / (2σ^2)) / w_xy
end

@inline σ_z(tf::IsotropicGaussian) = tf.C_axial * (tf.λ * tf.n / tf.NA^2) / √(2log(2))
function axialintensity(tf::IsotropicGaussian{3}, z::Length)
    σ = σ_z(tf)
    w_z = ustrip(√(2π) * σ) # NOTE: dimension gets integrated out from the normalization factor
    return exp(-z^2 / (2σ^2)) / w_z
end

function intensity(tf::IsotropicGaussian{3}, r::Length, z::Length)
    intensity(tf, r) * axialintensity(tf, z)
end

"""
    encircled_energy(tf::IsotropicGaussian, R::Length) 
Compute the encircled energy of the Gaussian PSF `tf` for a circle of radius `R`.

An isotropic Gaussian has a closed from encircled energy of ``E(R) = 1 - exp{-R² / 2σ²}``.

See also [`energy_radius`](@ref energy_radius(::IsotropicGaussian{2}, ::Real))
"""
encircled_energy(tf::IsotropicGaussian{2}, R::Length) = 1 - exp(-R^2 / (2σ_xy(tf)^2))

"""
    energy_radius(tf::IsotropicGaussian, ε::Real)
Calculate the energy radius of the Gaussian PSF `tf` for a given error term `ε`.

An isotropic Gaussian has a closed form energy radius of ``R(ε) = σ √(2 ln(1/ε))``.

See also [`encircled_energy`](@ref encircled_energy(::IsotropicGaussian{2}, ::Length))
"""
energy_radius(tf::IsotropicGaussian{2}, ε::Real) = σ_xy(tf) * √(2log(1 / ε))

params(tf::IsotropicGaussian{2}) = ComponentVector(λ=tf.λ, NA=tf.NA, C_lateral=tf.C_lateral)
params(tf::IsotropicGaussian{3}) = ComponentVector(λ=tf.λ, NA=tf.NA, n=tf.n, C_lateral=tf.C_lateral, C_axial=tf.C_axial)

export IsotropicGaussian
