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
AiryDisc{N,T}(; λ::Length{T}, NA::T, n::T=T(4 // 3)) where {N,T} = AiryDisc{N,T}(λ, NA, n) # final 1
AiryDisc{2,T}(params::ComponentVector) where {T} = AiryDisc{2,T}(;λ=params.λ, NA=T(params.NA)) # -> final 1
AiryDisc{3,T}(params::ComponentVector) where {T} = AiryDisc{3,T}(;λ=params.λ, NA=T(params.NA), n=T(params.n)) # -> final 1

"""
    AiryDisc{N}(λ::Length, NA, n=4//3)
[`AiryDisc`](@ref) model point spread function with emission wavelength `λ`, numerical aperture `NA` and immersion
medium refractive index `n`.
"""
function AiryDisc{N}(λ::Length{A}, NA::Real, n::Real) where {A<:Real,N} # final 2
    T = promote_type(A, typeof(NA), typeof(n))
    λ = convert(T, ustrip(λ)) * unit(λ)
    NA, n = convert(T, NA), convert(T, n)
    return AiryDisc{N,T}(λ, NA, n) # inner
end
AiryDisc{N}(; λ, NA, n=float(4 // 3)) where {N} = AiryDisc{N}(λ, NA, n) # -> final 2
AiryDisc{N}(λ::Length, NA::Real; kwargs...) where {N} = AiryDisc{N}(; λ, NA, kwargs...) # -> final 2
AiryDisc(args...; kwargs...) = AiryDisc{3}(args...; kwargs...) # -> final 2 (default to 3D)

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

@doc raw"""
    axialintensity(tf::AiryDisc{3}, z::Length)
Calculate the axial intensity of an [`AiryDisc`](@ref) PSF at the ``z``-offset `z`.

Axial intensity of the Airy disc is given by ``h(r) = \mathrm{sinc}(ϕ_z)²``  where ``ϕ_z = NA²z/λn``.
"""
function axialintensity(tf::AiryDisc{3,T}, z::Length) where {T}
    ϕ_z = (tf.NA^2 * z) / (tf.λ * tf.n) # normalized axial phase
    return iszero(z) ? oneunit(T) : sinc(ϕ_z)^2
end

"""
    intensity(tf::AiryDisc{3}, r::Length, z::Length)
Compute the intensity of the [`AiryDisc`](@ref) PSF at a distance `r` from its center and ``z``-offset `z`.

It is defined as the product of the [lateral intensity](@ref intensity(::AiryDisc{<:Any,T}, ::Length) where {T}) and
the [axial intensity](@ref axialintensity(::AiryDisc{3,T}, ::Length) where {T}).
"""
function intensity(tf::AiryDisc{3}, r::Length, z::Length)
    intensity(tf, r) * axialintensity(tf, z)
end

"""
    encircled_energy(tf::AiryDisc, R::Length) 
Compute the encircled energy of the Airy disc PSF `tf` for a circle of radius `R`.

An Airy disc has a closed from encircled energy of ``E(R) = 1 - J₀²(αR) - J₁²(αR)`` where ``α = 2πNA/λ``, ``J₀`` and
``J₁`` are the zeroth and first order Bessel functions of the first kind respectively.[Born - Principles of
Optics §8.5.2](@cite born2019a)

See also [`energy_radius`](@ref energy_radius(::AiryDisc{2}, ::Real))
"""
encircled_energy(tf::AiryDisc{2}, R::Length) = 1 - besselj0(2π*tf.NA*R/tf.λ)^2 - besselj1(2π*tf.NA*R/tf.λ)^2

"""
    energy_radius(tf::AiryDisc, ε::Real)
Calculate the energy radius of the Airy disc PSF `tf` for a given error term `ε`.

An Airy disc has an closed form encircled energy formula. The energy radius is computed by finding a root through
[bisection](@extref `Roots.Bisection`).
"""
energy_radius(tf::AiryDisc{2}, ε::Real) = find_zero(R -> encircled_energy(tf, R) - 1 + ε, (0.0u"nm", Inf*u"nm"), Bisection())

params(tf::AiryDisc{2}) = ComponentArray(λ=tf.λ, NA=tf.NA)
params(tf::AiryDisc{3}) = ComponentArray(λ=tf.λ, NA=tf.NA, n=tf.n)

function Base.show(io::IO, tf::AiryDisc)
    Base.showarg(io, tf, true) 
    params = [("λ", tf.λ), ("NA", tf.NA)]
    if tf isa AiryDisc{3}
        push!(params, ("n", tf.n))
    end
    print(io, "(", rounded(params...; sigdigits=3),  ")")
end

export AiryDisc
