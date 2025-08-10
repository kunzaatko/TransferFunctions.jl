using InterfaceFunctions
using TransferFunctions.Apodization: ApodizationFunction, apodization
using Base: Indices, tail, @propagate_inbounds
# FIX: Cannot be used with RGB arrays since we need to have the apodization with output type of
# `typeof(one(eltype(parent)))` instead of `eltype(parent)` in order to be able to multiply <08-08-25> 

const OneEdge = Tuple{Int,Int}
const Edges = NTuple{N,OneEdge} where {N}

"""
    TaperedArray{T,N,AA,IA,A} <: AbstractArray{T,N}
An `N`-dimensional array with data of type `T` that has edges of the width `.edges` tapered with an `ApodizationFunction`.

It is a lazy light wrapper that only computes the tapering when indexed.

It should be constructed using [`taperedges`](@ref).
"""
struct TaperedArray{S,T,N,AA,IA,A<:ApodizationFunction{S}} <: AbstractArray{T,N}
    parent::AA
    apodization::A
    edges::Edges{N}
    inner::IA
    function TaperedArray{S}(parent::AbstractArray{<:Any,N}, apodization::ApodizationFunction{S}, edges::NTuple{N,Tuple{Int,Int}}) where {S,N}
        # FIX: Should not be a problem. There is an issue in the checking of inclusion of the index in the `inner` range <08-08-25> 
        all(splat(+).(edges) .<= length.(axes(parent))) || throw(ArgumentError("Edges that sum up to more than the axes lengths lead to undefined behaviour."))
        inner = map((a, e) -> (first(a)+e[1]):(last(a)-e[2]), axes(parent), edges)
        new{S,eltype(parent),N,typeof(parent),typeof(inner),typeof(apodization)}(parent, apodization, edges, inner)
    end
end

Base.parent(a::TaperedArray) = (@inline; a.parent)
Base.size(a::TaperedArray) = (@inline; size(a.parent))
Base.axes(a::TaperedArray) = (@inline; axes(a.parent))
@inline attenuate(apo::ApodizationFunction, s::Int, e::Int, I::Int) = apodization(apo, (I - s) / (e - s))
function attenuate(apo::ApodizationFunction{T}, ax::AbstractUnitRange, edge::OneEdge, ind::Int) where {T}
    ind in ax && return oneunit(T)::T
    ind > last(ax) && return attenuate(apo, last(ax), last(ax) + edge[2], ind)::T
    return attenuate(apo, first(ax), first(ax) - edge[1], ind)::T
end
attenuate(apo::ApodizationFunction, inds::Indices, edges::Edges, I) = (@inline; prod(attenuate_recurse(apo, inds, edges, I)))
attenuate_recurse(apo::ApodizationFunction, inds::Indices, edges::Edges, indlast::NTuple{1}) = (attenuate(apo, inds[1], edges[1], indlast[1]),)
function attenuate_recurse(apo::ApodizationFunction{T}, inds::Indices, edges::Edges, I) where {T}
    a = attenuate(apo, inds[1], edges[1], I[1])
    (a, attenuate_recurse(apo, tail(inds), tail(edges), tail(I))...)
end
@propagate_inbounds function Base.getindex(a::TaperedArray{<:Any,<:Any,N}, I::Vararg{Int,N}) where {N}
    parent(a)[I...] * attenuate(a.apodization, a.inner, a.edges, I)
end

TaperedArray(parent::AbstractArray, apo::Type{<:ApodizationFunction}, edges) = TaperedArray(parent, apo(), edges)
function TaperedArray(parent::AbstractArray, apo::ApodizationFunction, edges) 
    S = typeof(one(eltype(parent)))
    return TaperedArray{S}(parent, convert(ApodizationFunction{S}, apo), edges)
end

inferedges(parent::AbstractArray{<:Any,N}, edge::Int) where {N} = ntuple(_ -> (edge, edge), N)
inferedges(parent::AbstractArray{<:Any,N}, edges::NTuple{N,Int}) where {N} = Tuple((e, e) for e in edges)

TaperedArray{S}(parent::AbstractArray, apo::ApodizationFunction, edges) where {S} = TaperedArray{S}(parent, apo, inferedges(parent, edges))

"""
    taperedges([apo=Cosine()], A, w, [border])
    taperedges([apo=Cosine()], A, (w1,w2...,wn), [border])
    taperedges([apo=Cosine()], A, ((w1a,w1b),...,(wna,wnb), [border])

Taper edges of width `w` of array `A` using [`apo::ApodizationFunction`](@ref ApodizationFunction). 

See also [`BorderArray`](@extref), [`TaperedArray`](@ref), [`ApodizationFunction`](@ref)
"""
taperedges(A::AbstractArray, args...) = taperedges(Apodization.Cosine(), A, args...)
taperedges(apo::ApodizationFunction, A::AbstractArray, w) = TaperedArray(A, apo, w)

function taperedges(apo::ApodizationFunction, A::AbstractArray, ws, border)
    return TaperedArray(BorderArray(A, border, ws), apo, ws)
end

export taperedges
