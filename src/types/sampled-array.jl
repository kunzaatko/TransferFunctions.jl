using Base: CartesianIndices
# TODO: This should instead be defined in some package like MicroscopyCore.jl or similar <30-07-25> 
"""
    SampledArray{T,ST,N,AA<:AbstractArray} <: AbstractArray{T,N}

An `N`-dimensional sampled array with data of type `T` sample __single__ spacing of type `ST`.

The restriction of `ST` being a single type means that `SampledArray` is not able to represent arrays that have
different dimensions in each direction.

See also [`SampledMatrix`](@ref), [`SampledVector`](@ref), [`SpatialArray`](@ref)
"""
struct SampledArray{T,ST,N,AA<:AbstractArray{T}} <: AbstractArray{T,N}
    parent::AA
    sampling::NTuple{N,ST}
    function SampledArray(parent::AbstractArray{T}, sampling::Tuple{ST,Vararg{ST}}) where {T,ST}
        ndims(parent) == length(sampling) || throw(DimensionMismatch("In constructor `SampledArray(parent::AbstractArray{T},sampling::NTuple)`. The length of the sampling rate `sampling` does not match the dimensionality of the parent array `parent`, `length(sampling)==$(length(sampling))!=$(ndims(parent))==ndims(parent)`."))
        return new{T,ST,ndims(parent),typeof(parent)}(parent, sampling)
    end
    SampledArray(parent::AbstractArray{<:Any,0}, sampling::NTuple{0,<:Any}) = new{eltype(parent),eltype(sampling),0,typeof(parent)}(parent, sampling)
end

Base.size(a::SampledArray) = (@inline; size(a.parent))
Base.axes(a::SampledArray) = (@inline; axes(a.parent))
Base.parent(a::SampledArray) = a.parent
Base.similar(a::SampledArray{T,ST,N}, ::Type{S}, dims::Dims{N}) where {T,ST,N,S}  = SampledArray(similar(parent(a), S, dims), a.sampling)
Base.getindex(a::SampledArray, i) = (@inline; getindex(parent(a), i))
Base.setindex!(A::SampledArray, v, i::Int) = (@inline; setindex!(parent(A), v, i))
Base.IndexStyle(::Type{<:SampledArray{T,ST,N,AA}}) where {T,ST,N,AA} = IndexStyle(AA)
sampling(a::SampledArray) = a.sampling

"""
    SampledMatrix{T, ST} <: AbstractMatrix{T}

Two-dimensional array with elements of type `T` with a given sample spacing of type `ST`. Alias for
`SampledArray{T,ST,2}`.
"""
const SampledMatrix{T,ST} = SampledArray{T,ST,2}


"""
    SampledVector{T, ST} <: AbstractVector{T}

One-dimensional array with elements of type `T` with a given sample spacing of type `ST`. Alias for
`SampledArray{T,ST,1}`.
"""
const SampledVector{T,ST} = SampledArray{T,ST,1}

"""
    SpatialArray{T,N} <: AbstractArray{T,N}

An `N`-dimensional array representing values sampled uniformly in space. Alias for `SampledArray{T,<:Length,N}`.

See also [`SpatialMatrix`](@ref), [`SpatialVector`](@ref), [`SampledArray`](@ref)
"""
const SpatialArray{T,N} = SampledArray{T,<:Length,N}

"""
    SpatialArray(A, (Δx, Δy,...))
    SpatialArray(A, Δ::Length)

Construct a `SpatialArray` with values `A` and sampling `(Δx, Δy, ...)`. For a single `Δ` uses uniform sampling with
distance `Δ` in every direction.
"""
SpatialArray(A::AbstractArray, Δ) = SampledArray(A, Δ)
SpatialArray(A::AbstractArray{<:Any,N}, Δ::Length) where {N} = SampledArray(A, fillsize(Δ, N))

"""
    SpatialMatrix{T} <: AbstractMatrix{T}

Two-dimensional array with elements of type `T` with a given sampling distance in both directions. Alias for
`SpatialArray{T,2}`.
"""
const SpatialMatrix{T} = SpatialArray{T,2}

"""
    SpatialMatrix(M, (Δx, Δy))
    SpatialMatrix(M, Δ::Length)

Construct a `SpatialMatrix` with values `M` and sampling `(Δx, Δy)`. For a single `Δ` uses uniform sampling with the
distance `Δ` in both directions.
"""
SpatialMatrix(A::AbstractMatrix, Δ) = SpatialArray(A, Δ)

"""
    SpatialVector{T} <: AbstractVector{T}
One-dimensional array with elements of type `T` with a given sampling distance. Alias for `SpatialArray{T,1}`.
"""
const SpatialVector{T} = SpatialArray{T,1}

"""
    SpatialVector(V, Δ::Length)

Construct a `SpatialVector` with values `V` and sampling `Δ`.
"""
SpatialVector(A::AbstractVector, Δ::Length) = SpatialArray(A, (Δ,))
