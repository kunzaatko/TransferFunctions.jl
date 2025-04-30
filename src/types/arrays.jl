using Base: @propagate_inbounds, Indices
using BlockArrays, ImageFiltering

"""
    SampledArray{T,ST,N,AA<:AbstractArray} <: AbstractArray{T,N}
An `N`-dimensional sampled array with data of type `T` sample spacing of type `ST`.

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
Base.similar(a::SampledArray{T,ST,N}, ::Type{S}, dims::Dims{N}) where {T,ST,N,S} = typeof(a)(similar(parent(a), S, dims), a.sampling)
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
SpatialArray(A::AbstractArray{<:Any,N}, Δ::Length) where {N} = SampledArray(A, fillsize(Δ, Val(N)))

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
    size::Size{N}
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
innersize(A::OuterInnerArray{<:Any,N}) where {N} = size(A)[A.isinnerdim]
innerlength(A::OuterInnerArray{<:Any,N}) where {N} = prod(innersize(A))
outersize(A::OuterInnerArray{<:Any,N}) where {N} = size(A)[A.isinnerdim.==false]
outerlength(A::OuterInnerArray{<:Any,N}) where {N} = prod(outersize(A))

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

See also [`FilteringMatrix`](@ref), `imfilter`
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
Base.similar(A::CirculantTensor{T}, eltype::Type{T}, dims::Dims) where {T} = (@inline; similar(parent(A), eltype, dims))

interior(inds::Indices{N}, others::Vararg{Indices{N}}) where {N} = map(intersect, inds, others...)
"""
    shrink(inds::Indices{N}, kernel::Indices{N})
Return "valid" indices for convolution of array with axes `inds` with kernel having axes `kern`.
"""
shrink(inds::Indices{N}, kern::Indices{N}) where {N} = map(shrinkind, inds, kern)
shrinkind(ind::AbstractUnitRange, kern::AbstractUnitRange) = typeof(ind)(first(ind)-first(kern):last(ind)-last(kern))
shrinkind(ind::Base.OneTo, kern::AbstractUnitRange) = shrinkind(UnitRange(ind), kern)

# NOTE: K = 2N should be satisfied by the constructor of `CirculantTensor` <26-04-25> 
# TODO: Decide what to do with the indices of the array. They will not be linearly spaced if it should make sense.
# Or it could be documented that the indices are not reflecting the actual indices of the reference array. <24-04-25> 
"""
    FilteringMatrix{T,K,P} <: AbstractMatrix{T}
A matrix that for a given kernel `K`, array size `size(A)` and a padding scheme `P` gives an array `F` such that
filtering (correlation or convolution if `reflect(K)` is used as the kernel) output of the kernel `K` and an array as
`A` can be computed as `F * A[:]`. It is a view into the [`CirculantTensor`](@ref) for the kernel array `K`.
"""
struct FilteringMatrix{T,K} <: AbstractMatrix{T}
    parent::OuterInnerArray{T,2,1,1, SubArray{}}
    array_dims::Size{K}
    function FilteringMatrix(A::CirculantTensor{<:Any,K}) where {K}
        @assert iseven(K) "`K` in a `CirculantTensor{<:Any,K}` must always be even"
        return new{eltype(A),K}(A)
    end
end

Base.parent(A::FilteringMatrix) = A.parent
function Base.size(A::FilteringMatrix)
    array_length = prod(A.array_dims)
    kern_length = length(parent(A))
    return (array_length, kern_length)
end

function Base.axes(A::FilteringMatrix)
    tensor_size = size(A.circulant)
    return (Base.OneTo(tensor_size[1] * tensor_size[3]), Base.OneTo(tensor_size[2] * tensor_size[4]))
end
BlockArrays.blockaxes(A::FilteringMatrix) = (BlockRange(axes(A.circulant, 1)), BlockRange(axes(A.circulant, 2)))
Base.getindex(A::FilteringMatrix, I::Block{1}) = error("TODO")

function conv_kern(A::AbstractMatrix{T}, s::NTuple{2,Int}) where {T}
    A_axs = axes(A)
    rows = Vector{SubArray}()
    for y in A_axs[2][begin:(end-s[2])]
        for x in A_axs[1][begin:(end-s[1])]
            push!(rows, view(view(A, x:(x+s[1]), y:(y+s[2])), :))
        end
    end
    return rows
end

public CirculantTensor, FilteringMatrix, OuterInnerArray
