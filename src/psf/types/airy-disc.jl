"""
    AiryDisc{N, T<:Real} <: PSFModel{N}
Airy disc scalar diffraction point spread function model for a circular aperture.

# Fields
- `λ::Length`: emission wavelength
- `NA::Real`: numerical aperture of the objective
- `n::Real`: refractive index of the immersion medium (defaults to `4//3` which is the refractive index of water)
"""
struct AiryDisc{N,T<:Real} <: PSFModel{N}
    λ::Length{T}
    NA::T
    n::T # refractive index of water
    function AiryDisc{N,T}(λ::Length{T}, NA::T, n::T) where {N,T}
        N ∈ (2, 3) || throw(ArgumentError("`AiryDisc` model can only be used for 2D and 3D PSFs. Got `N == $N`."))
        check_emission_wavelength(λ)
        check_numerical_aperture(NA)
        check_refractive_index(n)
        return new{N,T}(λ, NA, n)
    end
end
AiryDisc{N}(; λ, NA, n=float(4 // 3)) where {N} = AiryDisc{N}(λ, NA, n)
AiryDisc{N}(λ::Length, NA::Real; kwargs...) where {N} = AiryDisc{N}(; λ, NA, kwargs...)
AiryDisc(args...; kwargs...) = AiryDisc{3}(args...; kwargs...) # default to the 3D full model

"""
    AiryDisc{N}(λ::Length, NA, n=4//3)
[`AiryDisc`](@ref) model point spread function with emission wavelength `λ`, numerical aperture `NA` and immersion
medium refractive index `n`.
"""
function AiryDisc{N}(λ::Length{A}, NA::Real, n::Real) where {A<:Real,N}
    T = promote_type(A, typeof(NA), typeof(n))
    λ = convert(T, ustrip(λ)) * unit(λ)
    NA, n = convert(T, NA), convert(T, n)
    return AiryDisc{N,T}(λ, NA, n) # inner
end

symmetry(::AiryDisc) = ZAxisRadialSymmetry()

"""
    intensity(psf::AiryDisc, r::Length)
Compute the intensity of the [`AiryDisc`](@ref) PSF at a distance `r` from its center.

Intensity of the Airy disc is given by ``h(r) = (2J₁(αr)/(αr))²``  where ``α = 2πNA/λ`` and ``J₁`` is the first order
Bessel function of the first kind.
"""
function intensity(tf::AiryDisc{<:Any,T}, r::Length) where {T}
    α = 2π * tf.NA / tf.λ
    return iszero(r) ? oneunit(T) : (2 * besselj1(α * r) / (α * r))^2
end

"""
    axialintensity(tf::AiryDisc{3}, z::Length)
Calculate the axial intensity of an [`AiryDisc`](@ref) PSF at the ``z``-offset `z`.

Axial intensity of the Airy disc is given by ``h(r) = sinc(ϕ_z)²``  where ``ϕ_z = π NA²z/λn``.
"""
function axialintensity(tf::AiryDisc{3,T}, z::Length) where {T}
    ϕ_z = (π * tf.NA^2 * z) / (tf.λ * tf.n) # normalized axial phase
    return iszero(z) ? oneunit(T) : sinc(ϕ_z)^2
end

"""
    intensity(tf::AiryDisc{3}, r::Length, z::Length)
Compute the intensity of the [`AiryDisc`](@ref) PSF at a distance `r` from its center and ``z``-offset `z`.

It is defined as the product of the [lateral intensity](@ref `intensity(::AiryDisc{<:Any,T}, r::Length) where {T}`) and
the [axial intensity](@ref `axialintensity(::AiryDisc{3,T}, z::Length) where {T}`).
"""
function intensity(tf::AiryDisc{3}, r::Length, z::Length)
    intensity(tf, r) * axialintensity(tf, z)
end

export AiryDisc
