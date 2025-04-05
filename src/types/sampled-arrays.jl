"""
    SampledArray{T,ST,N,AA<:AbstractArray} <: AbstractArray{T,N}
An `N`-dimensional sampled array with data of type `T` sample spacing of type `ST`.
"""
struct SampledArray{T,ST,N,AA<:AbstractArray{<:T}} <: AbstractArray{T,N}
    parent::AA
    sampling::NTuple{N,ST}
end

Base.size(a::SampledArray) = Base.size(a.parent)
Base.axes(a::SampledArray) = Base.axes(a.parent)
Base.parent(a::SampledArray) = a.parent
Base.similar(a::SampledArray{T,ST,N}, ::Type{S}, dims::Dims{N}) where {T,ST,N,S} = typeof(a)(similar(parent(a), S, dims), a.sampling)
Base.getindex(a::SampledArray, i) = (@inline; Base.getindex(parent(a), i))
Base.IndexStyle(::Type{<:SampledArray{T,ST,N,AA}}) where {T,ST,N,AA} = Base.IndexStyle(AA)
sampling(a::SampledArray) = a.sampling

"""
    const SpatialArray{T,N} = SampledArray{T,<:Length,N}
An `N`-dimensional array representing values sampled uniformly in space.
"""
const SpatialArray{T,N} = SampledArray{T,<:Length,N}

"""
    SpatialArray(data::AbstractArray, Δ::Length)
Construct a SpatialArray with `data` values and a uniform sampling with the distance `Δ` in every direction.
"""
SpatialArray(data::AbstractArray, Δ) = SampledArray(data, Δ)
SpatialArray(data::AbstractArray{<:Any,N}, Δ::Length) where {N} = SampledArray(data, fillsize(Δ, Val(N)))
