using Base: OneTo
# FIX: When #58389 is merged into Base, I should reduce the number of printed type parameters <22-05-25> 
# TODO: I should fill in outer or inner dims if only one is provided. This should be done with something like this:
# `Tuple(n for n in OneTo(N) if !in(outer|inner)(n))` <18-05-25> 
"""
    Flattened{P,OM,IM,AX,IA,T,N} <: AbstractArray{T,N}

An `N`-dimensional `AbstractArray` stored as an `M`-dimensional 'outer' parent array `::P` of `K`-dimensional `T`-valued
'inner' arrays, where `N=M+K`.

It should be constructed by [`flatten`](@ref).

The constructor ensures the inner arrays have the same axes. In a sense, this can be thought as the inverse of
[`Slices`](@extref `Base.Slices`) from `Base`, which allows you to view a single array as multiple arrays. `Flattened` lets you see
multiple arrays that are nested in one container array as a single flat array.

[`parent(f::Flattened)`](@ref Flattened) will return the nested parent array.

See also [`Slices`](@extref `Base.Slices`)

# Fields
`outermap::OM` and `innermap::IM` are `M` and `K` integer long tuples respectively of the dimensions that the outer and
inner arrays represent. They are ensured to be unique and to comprise the total of `1:N` dimensions. `axes::AX` contains
the axes of the flattened array. `parent::P` holds the nested array.

# Examples
```jldoctest
julia> flatten(eachslice(reshape(1:12,3,4); dims=1)) == reshape(1:12,3,4)
true
```
"""
struct Flattened{P,OM,IM,AX,IA,T,N} <: AbstractArray{T,N}
    parent::P
    outermap::OM
    innermap::IM
    axes::AX
end

function Flattened(A::AbstractArray{IA}, outermap::OM, innermap::IM, axes::AX) where{T,OM,IM,AX,IA<:AbstractArray{T}}
    N = length(outermap) + length(innermap)
    P = typeof(A)
    Flattened{P,OM,IM,AX,IA,T,N}(A, outermap, innermap, axes) 
end

_flatten_check_dims(N) = nothing
function _flatten_check_dims(N, dim, dims...)
    1 <= dim <= N || throw(DimensionMismatch("Invalid dimension $dim"))
    dim in dims && throw(DimensionMismatch("Dimensions $dims are not unique"))
    _flatten_check_dims(N,dims...)
end
_flatten_check_axes(A) = allequal(axes, A) || throw(DimensionMismatch("Inner axes are not all the equal"))
function _flatten(A::AbstractArray{IA, NO}, outer::Dims{NO}, inner::Dims{NI}) where {T,NO,NI,IA<:AbstractArray{T,NI}}
    N = NO + NI
    _flatten_check_dims(N, inner..., outer...)
    _flatten_check_axes(A)
    ax = ntuple(Val(N)) do dim
        outerdim = findfirst(==(dim), outer)
        if !isnothing(outerdim) 
            axes(A)[outerdim]
        else # if not "outer" dim, it must be "inner", else bad arguments where supplied
            axes(first(A))[findfirst(==(dim), inner)]
        end
    end
    return Flattened(A, outer, inner, ax)
end

_default_outer(A) = Dims(1:ndims(A))
_default_inner(A) = Dims((ndims(A)+1):(ndims(A)+ndims(first(A))))

"""
    flatten(A; outer, inner)

Create a [`Flattened`](@ref) object that is a flat array of view of the nested array `A` with the dimensions of `A`
spanning `outer` dimensions of the resulting array and the dimensions of the elements of `A` spanning `inner` dimensions
of the resulting array.

See also [`eachslice`](@extref Base.eachslice)

# Examples
```jldoctest
julia> flatten([[1,2], [3,4]])
2×2 flatten(::Vector{Vector{Int64}}) with eltype Int64:
 1  2
 3  4

julia> flatten([[1,2], [3,4]]; outer=(2,), inner=(1,))
2×2 flatten(::Vector{Vector{Int64}}; outer=(2,), inner=(1,)) with eltype Int64:
 1  3
 2  4

julia> size(flatten([rand(2,3,4) for _ in 1:5, _ in 1:6])) # default is outer dimensions preceding inner dimensions
(5, 6, 2, 3, 4)
```
"""
@inline flatten(A; outer=_default_outer(A),inner=_default_inner(A)) = _flatten(A, outer, inner)

Base.axes(A::Flattened) = A.axes
Base.size(A::Flattened) = map(length, A.axes)

@inline function _inner_index(A::Flattened, c...)
    return map(l -> c[l], A.innermap)
end

@inline function _outer_index(A::Flattened, c...)
    return map(l -> c[l], A.outermap)
end

@inline function Base.getindex(A::Flattened{P,OM,IM,AX,IA,T,N}, I::Vararg{Int,N}) where {P,OM,IM,AX,IA,T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds @views A.parent[_outer_index(A,I...)...][_inner_index(A,I...)...]
end
@inline function Base.setindex!(A::Flattened{P,OM,IM,AX,IA,T,N}, val, I::Vararg{Int,N}) where {P,OM,IM,AX,IA,T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds A.parent[_outer_index(A,I...)...][_inner_index(A,I...)...] = val
end

Base.parent(s::Flattened) = s.parent

export flatten
