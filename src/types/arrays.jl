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
`N`-dimensional array that consists of an `M` dimensional 'outer' array of type `OA<:AbstractArray{IA,M}` of `K`
dimensional `T` valued 'inner' arrays of type `IA{T,K}`. The constructor ensures that `M+K==N` and that the inner arrays
have the same axes.

See also [`Slices`](@extref Julia :jl:type:`Base.Slices`)

# Fields
- `outer::OA`
- `isinnerdim::SVector{N,Bool}`
- `size::Size{N}`
- `axes::Indices{N}`
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

        OIA_size = Vector{Int}(undef, N)
        OIA_size[isinnerdim.==false] .= Base.size(OA)
        OIA_size[isinnerdim] .= Base.size(first(OA))
        OIA_size = Tuple(OIA_size)

        OIA_axes = Vector{AbstractUnitRange}(undef, N)
        OIA_axes[isinnerdim.==false] .= Base.axes(OA)
        OIA_axes[isinnerdim] .= Base.axes(first(OA))
        OIA_axes = Indices{N}(OIA_axes)

        return new{T,N,M,K,IA,typeof(OA)}(OA, isinnerdim, OIA_size, OIA_axes)
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
    outer_index = I[A.isinnerdim.==false]
    inner_index = I[A.isinnerdim]
    return A.outer[outer_index...][inner_index...]
end

"""
    CirculantTensor{T,N,M,AA} <: AbstractArray{T,N}
`N`-dimensional circulant tensor of an `M` dimensional array `A` or type `AA<:AbstractArray{T,M}`. The dimensionality
`N` is equal to `2M`, where the first `M` dimensions have interior filtering coordinates `ct.interior` given by the
kernel axes `ct.kern` and the axes of `A` and the tail `M` dimensions have the kernel coordinates `ct.kern`.

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
    function CirculantTensor{T,N}(A::AA, KI::Indices{M}) where {T,N,M,AA<:AbstractArray{T,M}}
        N == 2M || throw(DimensionMismatch("In constructor `CirculantTensor{T,N}(<:AbstractArray{T,M}, ::Indices{M})`. `N` must be equal to `2M`, but `2M==$(2M)!=$N==N`."))
        AI = interior(axes(A), KI)
        any(iszero, length(AI)) && throw(DimensionMismatch("In constructor `CirculantTensor(<:AbstractArray, ::Indices)` the array is not large enough for the kernel. Got interior of $AI.")) # TODO: Test
        views = map(CartesianIndices(AI)) do I
            @inbounds OffsetArray(view(A, CartesianIndices(KI) .+ I), KI)
        end
        parent = OffsetArray(views, AI) |> OuterInnerArray
        return new{T,N,M,AA}(A, AI, KI, parent)
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

# NOTE: K = 2N should be satisfied by the constructor of `CirculantTensor` <26-04-25> 
# TODO: Decide what to do with the indices of the array. They will not be linearly spaced if it should make sense.
# Or it could be documented that the indices are not reflecting the actual indices of the reference array. <24-04-25> 
"""
    FilteringMatrix{T,K,P} <: AbstractMatrix{T}
A matrix that for a given kernel `K`, array size `size(A)` and a padding scheme `P` gives an array `F` such that
filtering (correlation or convolution if `reflect(K)` is used as the kernel) output of the kernel `K` and an array as
`A` can be computed as `F * A[:]`. It is a view into the [`CirculantTensor`](@ref) for the kernel array `K`.
"""
struct FilteringMatrix{T,K,CT<:CirculantTensor{T}} <: AbstractMatrix{T}
    circulant::CT
    parent::OuterInnerArray{T,2,1,1}
    function FilteringMatrix(circulant::CirculantTensor{<:Any,N}) where {N}
        @assert iseven(N) "`K` in a `CirculantTensor{<:Any,K}` must always be even" # NOTE: Never should happen if the inner constructor is used for the CirculantTensor <05-05-25> 
        K = N ÷ 2
        rows = map(eachslice(circulant, dims=Tuple((K+1):N))) do s
            view(s, :)
        end
        parent = OuterInnerArray(view(rows, :))
        return new{eltype(circulant),K,typeof(circulant)}(circulant, parent)
    end
end
FilteringMatrix(args...) = FilteringMatrix(CirculantTensor(args...))

Base.parent(A::FilteringMatrix) = A.parent
Base.size(A::FilteringMatrix) = size(parent(A))
Base.axes(A::FilteringMatrix) = axes(parent(A))
Base.getindex(A::FilteringMatrix, ind...) = (@inline; getindex(parent(A), ind...))

CirculantTensor(A::FilteringMatrix) = A.circulant

# FIX: Should actually return `Indices` <05-05-25> 
# TODO: Change these to accessor functions (interface in base) <05-05-25> 
filtered_inds(A::FilteringMatrix) = A.circulant.interior
filtered_size(A::FilteringMatrix) = length.(filtered_inds(A))
kernel_inds(A::FilteringMatrix) = A.circulant.kern
kernel_size(A::FilteringMatrix) = length.(kernel_inds(A))

# FIX: It is not always a block matrix... It has blocks only if the sizes check out. The last block is not guaranteed to
# be a Toeplitz matrix <05-05-25> 
BlockArrays.blockaxes(A::FilteringMatrix) = (BlockRange(axes(A.circulant, 1)), BlockRange(axes(A.circulant, 2)))
Base.getindex(A::FilteringMatrix, I::Block{1}) = error("TODO")

