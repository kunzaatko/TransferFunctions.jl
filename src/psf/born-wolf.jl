"""
    BornWolf{T<:Real} <: RadialPSF
Born & Wolf model point spread function

# Fields
    + `λ::Length{T}`: emission wavelength 
    + `n::T`: index of refraction of the immersion medium
    + `NA::T`: numerical aperture of the objective
"""
@kwdef struct BornWolf{T<:Real} <: RadialPSF
    λ::Length{T}
    NA::T
    n_i::T = 4 // 3 # refractive index of water
    function BornWolf{T}(λ::Length{T}, NA::T, n_i::T) where {T}
        ustrip(λ) > zero(T) || throw(DomainError(λ, "Emission wavelength is a positive value. Got `λ = $λ`."))
        NA > zero(T) || throw(DomainError(NA, "Numerical aperture of the objective is a positive value. Got `NA = $NA`."))
        n_i > zero(T) || throw(DomainError(n_i, "Refractive index of the immersion is a positive value. Got `n_i = $n_i`."))
        return new{T}(λ, NA, n_i)
    end
end

"""
    BornWolf(λ::Length, NA, n_i=4//3)
[`BornWolf`](@ref) model point spread function with emission wavelength `λ`, numerical aperture `NA` and immersion medium
refractive index `n_i`.
"""
function BornWolf(λ::Length{A}, NA::B, n_i::C) where {A<:Real,B<:Real,C<:Real}
    T = promote_type(A, B, C)
    λ = convert(T, ustrip(λ)) * unit(λ)
    NA, n_i = convert(T, NA), convert(T, n_i)
    return BornWolf{T}(λ, NA, n_i)
end

function intensity(tf::BornWolf{T}, r::Length)::T where {T}
    k = 2π / tf.λ
    k₀ = k / tf.n_i
    # FIX: Is this correct?! <14-07-23> 
    # FIX: normalize <10-12-23> 
    return r == zero(r) ? 1 : (2besselj1(k₀ * r * tf.NA) / (k₀ * r * tf.NA))^2
end
