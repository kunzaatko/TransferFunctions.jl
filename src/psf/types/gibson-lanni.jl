"""
    GibsonLanni{T<:Real} <: PSFModel{2}
Gibson & Lanni model point spread function

# Fields
    + `λ::Length{T}`: emission wavelength
    + `NA::T`: numerical aperture of the objective
    + `n_i::T`: refractive index of the immersion medium
    + `n_s::T`: refractive index of the sample
    + `n_g::T`: refractive index of the coverslip
    + `t_i::Length{T}`: working distance of the objective
    + `t_g::Length{T}`: coverslip thickness
"""
@kwdef struct GibsonLanni{T<:Real} <: PSFModel{2}
    λ::Length{T}
    NA::T
    n_i::T = 1.5
    n_s::T = 4 // 3
    n_g::T = 1.5
    t_i::Length{T}
    t_g::Length{T}
    function GibsonLanni{T}(λ::Length{T}, NA::T, n_i::T, n_s::T, n_g::T, t_i::Length{T}, t_g::Length{T}) where {T}
        ustrip(λ) > zero(T) || throw(DomainError(λ, "Emission wavelength is a positive value. Got `λ = $λ`."))
        NA > zero(T) || throw(DomainError(NA, "Numerical aperture of the objective is a positive value. Got `NA = $NA`."))
        n_i > zero(T) || throw(DomainError(n_i, "Refractive index of the immersion is a positive value. Got `n_i = $n_i`."))
        return new{T}(λ, NA, n_i, n_s, n_g, t_i, t_g)
    end
end

"""
    GibsonLanni(λ::Length, NA, n_i=1.5, n_s=4//3, n_g=1.5, t_i, t_g)
[`GibsonLanni`](@ref) model point spread function with emission wavelength `λ`, numerical aperture `NA`, immersion
medium refractive index `n_i`, sample refractive index `n_s`, coverslip refractive index `n_g`, working distance `t_i`
and coverslip thickness `t_g`.
"""
function GibsonLanni(λ::Length{A}, NA::B, n_i::C, n_s::D, n_g::E, t_i::F, t_g::G) where {A<:Real,B<:Real,C<:Real,D<:Real,E<:Real,F<:Real,G<:Real}
    T = promote_type(A, B, C, D, E, F, G)
    λ = convert(T, ustrip(λ)) * unit(λ)
end
symmtery(::GibsonLanni) = ZAxisRadialSymmetry()

export GibsonLanni
