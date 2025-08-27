"""
    BornWolf{N,T<:Real} <: PSFModel{N}
Born & Wolf point spread function model.

# Fields
- `λ::Length`: emission wavelength
- `NA::Real`: numerical aperture of the objective
- `n::Real`: refractive index of the immersion medium (defaults to `4//3` which is the refractive index of water)
"""
struct BornWolf{N,T<:Real} <: PSFModel{N}
    λ::Length{T}
    NA::T
    f::T
    n::T
    function BornWolf{N,T}(λ::Length{T}, NA::T, f::T, n::T) where {N,T}
        check_emission_wavelength(λ)
        check_numerical_aperture(NA)
        check_refractive_index(n)
        f > zero(f) || throw(DomainError(n, "Refractive index of the immersion is a positive value. Got `n = $n`."))

        return new{N,T}(λ, NA, f, n)
    end
end
BornWolf{N,T}(; λ::Length{T}, NA::T, f::T, n::T=T(4 // 3)) where {N,T} = BornWolf{N,T}(λ, NA, f, n) # final 1
BornWolf{2,T}(params::ComponentVector) where {T} = BornWolf{2,T}(;λ=params.λ, NA=T(params.NA), f=T(params.f)) # -> final 1
BornWolf{3,T}(params::ComponentVector) where {T} = BornWolf{3,T}(;λ=params.λ, NA=T(params.NA), f=T(params.f), n=T(params.n)) # -> final 1

"""
    BornWolf(λ::Length, NA, f, n=4//3)
[`BornWolf`](@ref) model point spread function with emission wavelength `λ`, numerical aperture `NA` and immersion medium
refractive index `n`.
"""
function BornWolf{N}(λ::Length{A}, NA::Real, f::Real, n::Real) where {N,A<:Real} # final 2
    T = promote_type(A, typeof(NA), typeof(n))
    λ = convert(T, ustrip(λ)) * unit(λ)
    NA, f, n = convert(T, NA), convert(T, f), convert(T, n)
    return BornWolf{N,T}(λ, NA, f, n)
end
BornWolf{N}(;λ, f, NA, n=float(4 // 3)) where {N} = BornWolf{N}(λ,NA,f,n) # -> final 2
BornWolf{N}(λ::Length, NA::Real, f::Real; kwargs...) where {N} = BornWolf{N}(; λ, f, NA, kwargs...) # -> final 2
BornWolf(args...; kwargs...) = BornWolf{3}(args...; kwargs...) # -> final 2 (default to the 3D full model)

symmetry(::BornWolf) = ZAxisRadialSymmetry()

function intensity(tf::BornWolf{<:Any,T}, r::Length)::T where {T}
    k = 2π / tf.λ
    k₀ = k / tf.n
    return r == zero(r) ? oneunit(T) : (2besselj1(k₀ * r * tf.NA) / (k₀ * r * tf.NA))^2
end

params(tf::BornWolf{2}) = ComponentArray(λ=tf.λ, NA=tf.NA, f=tf.f)
params(tf::BornWolf{3}) = ComponentArray(λ=tf.λ, NA=tf.NA, f=tf.f, n=tf.n)

export BornWolf
