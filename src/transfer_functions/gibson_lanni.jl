@doc raw"""
Gibson & Lanni model of the transfer function for a circular aperture.

# Parameters

# Extended help

The Gibson & Lanni model is assumes that, disregarding defocus, all observed aberrations are generated by factors external to the objective (i.e. originating in the sample, coverslip and immersion medium combination). These aberrations can be characterized by the optical path difference between a ray in a perfect system (see [`BornWolf`](@ref)) and a ray under experimental conditions.
"""
Base.@kwdef struct GibsonLanni{T<:Real} <: ClosedFormPSFModel
    "numrical aperture"
    NA::T
    "wavelength"
    λ::T
    "immersion medium refractive index"
    n_i::T = 1.5
    "sample refractive index"
    n_s::T = (1 + 1 // 3)
    "coverslip refractive index"
    n_g::T = 1.5
    "working distance of the objective"
    t_i::T
    "coverslip thickness"
    t_g::T
end

@traitimpl SymmetricPupilFunction{GibsonLanni}
