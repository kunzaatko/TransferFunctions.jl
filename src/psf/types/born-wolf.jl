"""
    BornWolf{T<:Real} <: PSFModel{3}
Born & Wolf point spread function model.

# Fields
- `λ::Length`: emission wavelength
- `NA::Real`: numerical aperture of the objective
- `n::Real`: refractive index of the immersion medium (defaults to `4//3` which is the refractive index of water)
"""
@kwdef struct BornWolf{T<:Real} <: PSFModel{3}
    λ::Length{T}
    NA::T
    f::T
    n::T = 4 // 3 # refractive index of water
    function BornWolf{T}(λ::Length{T}, NA::T, f::T, n::T) where {T}
        λ > zero(λ) || throw(DomainError(λ, "Emission wavelength is a positive value. Got `λ = $λ`."))
        NA > zero(NA) || throw(DomainError(NA, "Numerical aperture of the objective is a positive value. Got `NA = $NA`."))
        f > zero(f) || throw(DomainError(n, "Refractive index of the immersion is a positive value. Got `n = $n`."))
        n > zero(n) || throw(DomainError(n, "Refractive index of the immersion is a positive value. Got `n = $n`."))
        return new{T}(λ, NA,f, n)
    end
end

"""
    BornWolf(λ::Length, NA, f, n=4//3)
[`BornWolf`](@ref) model point spread function with emission wavelength `λ`, numerical aperture `NA` and immersion medium
refractive index `n`.
"""
function BornWolf(λ::Length{A}, NA::Real, f::Real, n::Real) where {A<:Real}
    T = promote_type(A, typeof(NA), typeof(n))
    λ = convert(T, ustrip(λ)) * unit(λ)
    NA,f, n = convert(T, NA), convert(T,f), convert(T, n)
    return BornWolf{T}(λ, NA,f, n)
end
BornWolf(λ::Length, f::Real, NA::Real) = BornWolf(;λ,f, NA)
symmetry(::BornWolf) = ZAxisRadialSymmetry()

function intensity(tf::BornWolf{T}, r::Length)::T where {T}
    k = 2π / tf.λ
    k₀ = k / tf.n
    return r == zero(r) ? oneunit(T) : (2besselj1(k₀ * r * tf.NA) / (k₀ * r * tf.NA))^2
end

export BornWolf
