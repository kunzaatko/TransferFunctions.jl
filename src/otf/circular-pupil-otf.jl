# TODO: Update the signature of the type in the docstring. <11-04-25> 
"""
    CircularPupilOTF(ρ₀::Real)

Ideal (aberration free) OTF of a diffraction limited imaging system with incoherent light with the cutoff-frequency `ρ₀`
"""
Base.@kwdef struct CircularPupilOTF{T<:Real} <: RadialOTF{2}
    λ::Length{T}
    NA::T
    n::T = float(4 // 3)
    curvature::T
    function CircularPupilOTF(λ::Length{A}, NA::B, n::C, curvature::D) where {A<:Real,B<:Real,C<:Real,D<:Real}
        oneunit(curvature) >= curvature > zero(curvature) || throw(DomainError(curvature, "Valid domain for curvature is (0,1]"))
        check_emission_wavelength(λ)
        check_numerical_aperture(NA)
        check_refractive_index(n)
        T = promote_type(A, B, C, D)
        λ = convert(T, ustrip(λ)) * unit(λ)
        return new{T}(λ, NA, n, curvature)
    end
end

attenuation_normalized(tf::CircularPupilOTF{T}, ν::Number) where {T} = ν >= oneunit(ν) ? zero(T) : (2 / π) * (acos(ν) - ν * sqrt(1 - ν * ν)) * tf.curvature^ν
# TODO: Is the immersion refractive index necessary? <14-07-23> 
attenuation(tf::CircularPupilOTF, f_r::Frequency) = attenuation_normalized(tf, (f_r * tf.λ) / (2 * tf.n * tf.NA))

function cutoff(tf::CircularPupilOTF{T}, a::Real)::Frequency where {T}
    if iszero(a)
        return (2 * tf.NA * tf.n) / tf.λ
    else
        0 <= a <= 1 || throw(DomainError(a, "OTF can only attenuate with a coefficient between 0 and 1"))
        z = find_zero(ν -> attenuation_normalized(tf, ν) - a, (0, 1), Bisection())
        return (2z * tf.n * tf.NA) / tf.λ
    end
end

export CircularPupilOTF
