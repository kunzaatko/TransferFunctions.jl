using Base: @propagate_inbounds

reflect(ax::AbstractUnitRange) = -last(ax):-first(ax)

"""
    ReflectedArray{T,N,AA}
An array wrapper that reflects the parent around the origin in all dimensions.
"""
struct ReflectedArray{T,N,AA<:AbstractArray{T,N}} <: AbstractArray{T,N}
    parent::AA
end

Base.parent(a::ReflectedArray) = (@inline; a.parent)
Base.size(a::ReflectedArray) = (@inline; size(parent(a)))
Base.axes(a::ReflectedArray) = (@inline; map(reflect, axes(parent(a))))
@propagate_inbounds function Base.getindex(a::ReflectedArray{<:Any,N}, I::Vararg{Int,N}) where {N}
    @boundscheck checkbounds(a, I...)
    @inbounds parent(a)[((-1) .* I)...]
end

"""
    reflect(a::AbstractArray)
Returns the [`ReflectedArray`](@ref) of `a`.

For use when convolution filtering instead of correlation filtering.
"""
reflect(a::AbstractArray) = ReflectedArray(a)

"""
    ReflectedMatrix{T} <: AbstractMatrix{T}
Two-dimensional reflected array with elements of type `T`. Alias for `ReflectedArray{T,2}`.
"""
const ReflectedMatrix{T} = ReflectedArray{T,2}

"""
    ReflectedVector{T} <: AbstractVector{T}
One-dimensional reflected array with elements of type `T`. Alias for `ReflectedArray{T,1}`.
"""
const ReflectedVector{T} = ReflectedArray{T,1}

export reflect
