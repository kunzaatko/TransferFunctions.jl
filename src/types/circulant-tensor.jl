using Base: @propagate_inbounds, Indices
using BlockArrays, StaticArraysCore, ImageFiltering

const SBitVector{N} = SVector{N,Bool}

# TODO: If `innerdims` are not sorted, it should permute the dims of the inner arrays. <16-04-25> 
"""
    OuterInnerArray{T,N,M,K,IA,OA} <: AbstractArray{T,N}
`N`-dimensional array that consists of an `M` dimensional 'outer' array `OA<:AbstractArray{IA,M}` of `K` dimensional `T`
valued 'inner' arrays `IA{T,K}`. The constructor ensures that `M+K==N` and that the inner arrays have the same axes.

See also [`Slices`](@extref Julia :jl:type:`Base.Slices`)

# Fields
- `outer::OA` -- Array containing 
- `isinnerdim::SVector{N,Bool}`
"""
struct OuterInnerArray{T,N,M,K,IA<:AbstractArray{<:T,K},OA<:AbstractArray{IA,M}} <: AbstractArray{T,N}
    outer::OA
    isinnerdim::SBitVector{N}
    size::NTuple{N,Int}
    axes::Indices{N}
    OuterInnerArray(OA::AbstractArray{<:IA}, innerdims::Dims) where {T,IA<:AbstractArray{T}} = OuterInnerArray{T}(OA, innerdims)
    OuterInnerArray{T}(OA::AbstractArray{<:AbstractArray{<:T,K},M}, innerdims::Dims{K}) where {T,K,M} = OuterInnerArray{T,(K + M)}(OA, innerdims)
    function OuterInnerArray{T,N}(OA::AbstractArray{<:IA,M}, innerdims::Dims{K}) where {T,N,M,K,IA}
        N == (K + M) || throw(DimensionMismatch("In constructor `OuterInnerArray{T,N}(<:AbstractArray{<:AbstractArray{<:T,K},M},<:Dims{K})`. Inner dimensionality plus outer dimensionality `K+M` must be equal to the total number of dimensions `N`, but got `K==$K`, `M==$M`, `N==$N` hence `K+M==$(K+M)!=$N==N`."))
        allequal(Base.axes.(OA)) || throw(DimensionMismatch("In constructor `OuterInnerArray{T,N}(OA::AbstractArray{IA, M},<:Dims{K})`. All of the inner arrays in must have the same axes."))
        maximum(innerdims) <= N || throw(DimensionMismatch("In constructor `OuterInnerArray{T,N}(<:AbstractArray{IA,M},d::Dims{K})`. Maximum of `d` must be less than `N`, but got `maximum(d)==$(maximum(innerdims))>$N==N`."))
        allunique(innerdims) || throw(ArgumentError("In constructor `OuterInnerArray{T,N}(<:AbstractArray, innerdims::Dims{K})`, inner dimensions `innerdims` must be unique."))

        isinnerdim = SVector(ntuple(n -> ifelse(n ∈ innerdims, true, false), Val(N)))

        array_size = Vector{Int}(undef, N)
        array_size[isinnerdim.==false] .= Base.size(OA)
        array_size[isinnerdim] .= Base.size(first(OA))
        array_size = Tuple(array_size)

        array_axes = Vector{AbstractUnitRange}(undef, N)
        array_axes[isinnerdim.==false] .= Base.axes(OA)
        array_axes[isinnerdim] .= Base.axes(first(OA))
        array_axes = Indices{N}(array_axes)

        return new{T,N,M,K,IA,typeof(OA)}(OA, isinnerdim, array_size, array_axes)
    end
end
OuterInnerArray(OA::AbstractArray{<:IA,M}) where {M,K,IA<:AbstractArray{<:Any,K}} = OuterInnerArray(OA, Dims((M+1):(M+K)))

function outerdims(A::OuterInnerArray{<:Any,N}) where {N}
    ntuple(identity, Val(N))[A.isinnerdim.==false]
end
function innerdims(A::OuterInnerArray{<:Any,N}) where {N}
    ntuple(identity, Val(N))[A.isinnerdim]
end

Base.size(A::OuterInnerArray) = A.size
Base.axes(A::OuterInnerArray) = A.axes
@inline @propagate_inbounds function Base.getindex(A::OuterInnerArray{T,N}, I::Vararg{Int,N}) where {T,N}
    outer_index = CartesianIndex(I[A.isinnerdim.==false])
    inner_index = CartesianIndex(I[A.isinnerdim])
    return A.outer[outer_index][inner_index]
end

"""
    CirculantTensor{T,N,M,AA} <: AbstractArray{T,N}
`N`-dimensional circulant tensor of the `M` dimensional array `A::AA`. The dimensionality `N` is equal to `2M`, where
the first `M` dimensions have interior filtering coordinates `interior` given by the kernel axes `kern` and the axes of
`A` and the tail `M` dimensions have the kernel coordinates `kern`.

A correlation filtering result of `A` with a kernel array `K`can be obtained by outer tensor contraction over the tail
`M` dimensions of the `CirculantTensor(A,K)`.

See also [`BlockCirculantWithCirculantBlocksMatrix`](@ref), `imfilter`
"""
struct CirculantTensor{T,N,M,AA<:AbstractArray{T,M}} <: AbstractArray{T,N}
    A::AA
    interior::Indices{M}
    kern::Indices{M}
    parent::OuterInnerArray{T,N,M,M}
    CirculantTensor(A::AbstractArray{T,M}, kern::Indices{M}) where {T,M} = CirculantTensor{T,2M}(A, kern)
    function CirculantTensor{T,N}(A::AA, kern::Indices{M}) where {T,N,M,AA<:AbstractArray{T,M}}
        N == 2M || throw(DimensionMismatch("In constructor `CirculantTensor{T,N}(<:AbstractArray{T,M}, ::Indices{M})`. `N` must be equal to `2M`, but `2M==$(2M)!=$N==N`."))
        interior_inds = interior(axes(A), shrink(axes(A), kern))
        any(iszero, length(interior_inds)) && throw(DimensionMismatch()) # TODO
        views = map(CartesianIndices(interior_inds)) do cind
            @inbounds view(A, CartesianIndices(kern) .+ cind)
        end
        views = OffsetArray(views, interior_inds...)
        parent = OuterInnerArray(views)
        return new{T,N,M,AA}(A, interior_inds, kern, parent)
    end
end

@reexport using ImageFiltering: reflect # NOTE: For convolution instead of correlation <24-04-25> 

# TODO: Same constructors as for `imfilter!` <19-04-25>
CirculantTensor(A::AbstractArray{<:Any,N}, kern::AbstractArray{<:Any,N}) where {N} = CirculantTensor(A, axes(kern))
CirculantTensor(A::AbstractArray{<:Any,N}, size::NTuple{N,Int}) where {N} = CirculantTensor(A, map(Base.OneTo, size))
Base.parent(A::CirculantTensor) = (@inline; A.parent)
Base.size(A::CirculantTensor) = (@inline; size(parent(A)))
Base.size(A::CirculantTensor, dim) = (@inline; size(parent(A), dim))
Base.getindex(A::CirculantTensor, ind...) = (@inline; getindex(parent(A), ind...))
Base.axes(A::CirculantTensor, ind...) = (@inline; axes(parent(A), ind...))
Base.similar(A::CirculantTensor, args...) = (@inline; similar(parent(A), args...))

interior(inds::Indices{N}, others::Vararg{Indices{N}}) where {N} = map(intersect, inds, others...)
"""
    shrink(inds::Indices{N}, kernel::Indices{N})
Return "valid" indices for convolution of array with axes `inds` with kernel having axes `kern`.
"""
shrink(inds::Indices{N}, kern::Indices{N}) where {N} = map(shrinkind, inds, kern)
shrinkind(ind::AbstractUnitRange, kern::AbstractUnitRange) = typeof(ind)(first(ind)-first(kern):last(ind)-last(kern))
shrinkind(ind::Base.OneTo, kern::AbstractUnitRange) = shrinkind(UnitRange(ind), kern)

# TODO: Decide what to do with the indices of the array. They will not be linearly spaced if it should make sense.
# Or it could be documented that the indices are not reflecting the actual indices of the reference array. <24-04-25> 
"""
    BlockCirculantWithCirculantBlocksMatrix{T} <: AbstractBlockMatrix{T}
""" # TODO: Docs <24-04-25> 
struct BlockCirculantWithCirculantBlocksMatrix{T} <: AbstractBlockMatrix{T}
    circulant::CirculantTensor{T,4}
    parent::SubArray{T,2,CirculantTensor{T,4}}
    # function BlockCirculantWithCirculantBlocksMatrix(circulant::CirculantTensor{T,4})
    #     parent = view(circulant, )
    # end
end

function conv_kern(A::AbstractMatrix{T}, d::Dims{2}) where {T}
    A_axs = axes(A)
    rows = Vector{SubArray}()
    for y in A_axs[2][begin:(end-d[2])]
        for x in A_axs[1][begin:(end-d[1])]
            push!(rows, view(view(A, x:(x+d[1]), y:(y+d[2])), :))
        end
    end
    return rows
end
