using Base: CartesianIndices, @propagate_inbounds

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
    function SampledArray(parent::AbstractArray{T,N}, sampling::NTuple{N,Number}) where {T,N}
        # PERF: defined on `Number` but it promotes to a concrete type stored in the field
        sampling = promote(sampling...)
        return new{T,eltype(sampling),ndims(parent),typeof(parent)}(parent, sampling)
    end
    SampledArray(parent::AbstractArray{<:Any,0}, sampling::NTuple{0}) = new{eltype(parent),Any,0,typeof(parent)}(parent, sampling)
end
SampledArray(parent::AbstractArray{<:Any,N}, sampling::Number) where {N} = SampledArray(parent, ntuple(_ -> sampling, Val(N)))

Base.size(a::SampledArray) = (@inline; size(a.parent))
Base.axes(a::SampledArray) = (@inline; axes(a.parent))
Base.parent(a::SampledArray) = a.parent
Base.similar(a::SampledArray{T,ST,N}, ::Type{S}, dims::Dims{N}) where {T,ST,N,S} = SampledArray(similar(parent(a), S, dims), a.sampling)
@propagate_inbounds Base.getindex(a::SampledArray, i) = getindex(parent(a), i)
@propagate_inbounds Base.setindex!(A::SampledArray, v, i) = setindex!(parent(A), v, i)
Base.IndexStyle(::Type{<:SampledArray{<:Any,<:Any,<:Any,AA}}) where {AA} = IndexStyle(AA)
@inline sampling(a::SampledArray) = a.sampling

# FIX: Instead should be `location_bins` which give a vector of rectangles that are the bin location corners of the
# samples <30-07-25> 
"""
    posgrid(a::SampledArray)
Returns the grid of positions of the samples of `a`.

```jldoctest
julia> sa = SampledArray(reshape(1:16, (4,4)), 61u"nm");

julia> xs, ys = TF.posgrid(sa);

julia> xs 
4×4 Matrix{Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}:
   0.0 nm    0.0 nm    0.0 nm    0.0 nm
  61.0 nm   61.0 nm   61.0 nm   61.0 nm
 122.0 nm  122.0 nm  122.0 nm  122.0 nm
 183.0 nm  183.0 nm  183.0 nm  183.0 nm

julia> ys 
4×4 Matrix{Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}:
 0.0 nm  61.0 nm  122.0 nm  183.0 nm
 0.0 nm  61.0 nm  122.0 nm  183.0 nm
 0.0 nm  61.0 nm  122.0 nm  183.0 nm
 0.0 nm  61.0 nm  122.0 nm  183.0 nm
```
"""
@inline posgrid(a::SampledArray) = posgrid(size(a), sampling(a); center=(1, 1))

"""
    sample_vertices(a::SampledArray)
Returns the vertices of the samples of `a` as tuples.

```jldoctest
julia> sa = SampledArray(reshape(1:16, (4,4)), 61u"nm");

julia> TF.sample_vertices(sa)
4×4 Matrix{Tuple{Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}, Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}}:
 (0.0 nm, 0.0 nm)    (0.0 nm, 61.0 nm)    (0.0 nm, 122.0 nm)    (0.0 nm, 183.0 nm)
 (61.0 nm, 0.0 nm)   (61.0 nm, 61.0 nm)   (61.0 nm, 122.0 nm)   (61.0 nm, 183.0 nm)
 (122.0 nm, 0.0 nm)  (122.0 nm, 61.0 nm)  (122.0 nm, 122.0 nm)  (122.0 nm, 183.0 nm)
 (183.0 nm, 0.0 nm)  (183.0 nm, 61.0 nm)  (183.0 nm, 122.0 nm)  (183.0 nm, 183.0 nm)
```
"""
@inline sample_vertices(a::SampledArray) =
    map(posgrid(a)...) do x, y
        (x, y)
    end

"""
    SampledMatrix{T, ST} <: AbstractMatrix{T}
Two-dimensional array with elements of type `T` with a given sample spacing of type `ST`. Alias for
`SampledArray{T,ST,2}`.
"""
const SampledMatrix{T,ST} = SampledArray{T,ST,2}

"""
    SampledMatrix(M, (Δx, Δy))
    SampledMatrix(M, Δ)
Construct a `SampledMatrix` with values `M` and sampling `(Δx, Δy)`. For a single `Δ` uses uniform sampling with the
distance `Δ` in both directions.
"""
SampledMatrix(A::AbstractMatrix, Δ) = SampledArray(A, Δ)

"""
    SampledVector{T, ST} <: AbstractVector{T}
One-dimensional array with elements of type `T` with a given sample spacing of type `ST`. Alias for
`SampledArray{T,ST,1}`.
"""
const SampledVector{T,ST} = SampledArray{T,ST,1}

"""
    SampledVector(V, Δ)
Construct a `SampledVector` with values `V` and sampling `Δ`.
"""
SampledVector(A::AbstractVector, Δ::Length) = SampledArray(A, (Δ,))

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

export SpatialArray, SpatialMatrix, SpatialVector, SampledArray
