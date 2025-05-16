using Base: @propagate_inbounds, Indices
using BlockArrays, ImageFiltering
using TensorOperations, LinearAlgebra

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

# TODO: CirculantTensor should store the border type. We could optimize for "reflect" using fft <10-05-25> 
"""
    CirculantTensor{T,N,M,AA} <: AbstractArray{T,N}
`N`-dimensional circulant tensor of an `M` dimensional array `A` or type `AA<:AbstractArray{T,M}`. The dimensionality
`N` is equal to `2M`, where the first `M` dimensions have interior filtering coordinates `ct.interior` given by the
kernel axes `ct.kern` and the axes of `A` and the tail `M` dimensions have the kernel coordinates `ct.kern`.

A correlation filtering result of `A` with a kernel array `K`can be obtained by outer tensor contraction over the tail
`M` dimensions of the `CirculantTensor(A,K)`.

See also [`FilteringMatrix`](@ref), `imfilter`
"""
struct CirculantTensor{T,M,N,AA<:AbstractArray{T,M}} <: AbstractArray{T,N}
    A::AA
    interior::Indices{M}
    kern::Indices{M}
    parent::OuterInnerArray{T,N,M,M}
    CirculantTensor(A::AbstractArray{T,M}, kern::Indices{M}) where {T,M} = CirculantTensor{T,M}(A, kern)
    function CirculantTensor{T,M}(A::AA, KI::Indices{M}) where {T,M,AA<:AbstractArray{T,M}}
        # N == 2M || throw(DimensionMismatch("In constructor `CirculantTensor{T,N}(<:AbstractArray{T,M}, ::Indices{M})`. `N` must be equal to `2M`, but `2M==$(2M)!=$N==N`."))
        AI = interior(axes(A), KI)
        any(iszero, length(AI)) && throw(DimensionMismatch("In constructor `CirculantTensor(<:AbstractArray, ::Indices)` the array is not large enough for the kernel. Got interior of $AI.")) # TODO: Test
        views = map(CartesianIndices(AI)) do I
            @inbounds OffsetArray(view(A, CartesianIndices(KI) .+ I), KI)
        end
        parent = OffsetArray(views, AI) |> OuterInnerArray
        return new{T,M,2M,AA}(A, AI, KI, parent)
    end
end

@reexport using ImageFiltering: reflect # NOTE: For convolution instead of correlation <24-04-25> 

# Step 1: Determine kernel indices
@inline CirculantTensor(A::AbstractArray{<:Any,N}, kern::AbstractArray{<:Any,N}, args...) where {N} = CirculantTensor(A, axes(kern), args...)

# Step 2: Determine border and Initialize it if it is not fully specified. (Default Inner() doesn't need action)
@inline function CirculantTensor(A::AbstractArray{<:Any,N}, kern::Indices{N}, border::AbstractString, args...) where {N}
    return CirculantTensor(A, kern, ImageFiltering.borderinstance(border), args...)
end
@inline function CirculantTensor(A::AbstractArray{<:Any,N}, kern::Indices{N}, border::Pad{0}, args...) where {N}
    border = Pad(border.style, kern)
    return CirculantTensor(A, kern, border, args...)
end
@inline function CirculantTensor(A::AbstractArray{<:Any,N}, kern::Indices{N}, border::Fill{T,0}, args...) where {N,T}
    border = Fill(border.value, kern)
    return CirculantTensor(A, kern, border, args...)
end

# Step 3: Apply the border and call the inner constructor
function CirculantTensor(A::AbstractArray{T,N}, kern::Indices{N}, border::AbstractBorder, args...) where {T,N}
    A = padarray(T, A, border)
    return CirculantTensor(A, kern) # TODO: What should happen with the args... ? <10-05-25> 
end

function conv(T::CirculantTensor, K::AbstractArray)
    @assert T.kern == axes(K)
    return _conv(T, K)
end

# PERF!: FFT should be used for a much faster implementation using the CirculantTensor structure <13-05-25> 
function _conv(T::CirculantTensor{<:Any,1}, K::AbstractVector)
    return OffsetVector(no_offset_view(T) * no_offset_view(K), T.interior)
end
function _conv(T::CirculantTensor{<:Any,2}, K::AbstractMatrix)
    @tensor A[a, b] := no_offset_view(T)[a, b, c, d] * no_offset_view(K)[c, d]
    return OffsetMatrix(A, T.interior)
end
function _conv(T::CirculantTensor{<:Any,3}, K::AbstractArray{<:Any,3})
    @tensor A[a, b, c] := no_offset_view(T)[a, b, c, d, e, f] * no_offset_view(K)[d, e, f]
    return OffsetArray(A, T.interior)
end

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
    FilteringMatrix{T,M,P} <: AbstractMatrix{T}
A matrix that for a given kernel `K`, array `A` of size `size(A)==size(K)` equal to `M` and a padding scheme `P` gives
an array `F` such that filtering (correlation or convolution if `reflect(K)` is used as the kernel) output of the kernel
`K` and an array as `A` can be computed as `F * A[:]`. It is a view into the [`CirculantTensor`](@ref) for the kernel array `K`.
"""
struct FilteringMatrix{T,K,CT<:CirculantTensor{T}} <: AbstractMatrix{T}
    circulant::CT
    parent::OuterInnerArray{T,2,1,1}
    function FilteringMatrix(circulant::CirculantTensor{<:Any,M,K}) where {M,K}
        # NOTE: Never should happen if the inner constructor is used for the CirculantTensor <05-05-25> 
        @assert iseven(K) "`K` in a `CirculantTensor{<:Any,M,K}` must always be even"
        rows = map(eachslice(circulant, dims=Tuple((M+1):K))) do s
            reshape(s, :)  # PERF: `reshape` benchmarks better than `view` for 2D matrices <11-05-25> 
        end
        parent = OuterInnerArray(view(rows, :), (1,))
        return new{eltype(circulant),M,typeof(circulant)}(circulant, parent)
    end
end
FilteringMatrix(args...) = FilteringMatrix(CirculantTensor(args...))

Base.parent(A::FilteringMatrix) = A.parent
Base.size(A::FilteringMatrix) = size(parent(A))
Base.axes(A::FilteringMatrix) = axes(parent(A))
Base.getindex(A::FilteringMatrix, ind...) = (@inline; getindex(parent(A), ind...))

CirculantTensor(A::FilteringMatrix) = A.circulant

Base.propertynames(::FilteringMatrix) = (:circulant, :parent, :Aaxes, :Kaxes)
function Base.getproperty(A::FilteringMatrix, s::Symbol)
    if s === :Aaxes
        return getfield(A, :circulant).interior
    elseif s === :Kaxes
        return getfield(A, :circulant).kern
    else
        return getfield(A, s)
    end
end

# NOTE: Unlike `stdlib` we want matmul to return the offset vectors or matrices based on the input
LinearAlgebra.matprod_dest(A::FilteringMatrix, ::AbstractVector, T::Type) = OffsetVector(Vector{T}(undef, size(A, 1)), axes(A, 2))
LinearAlgebra.matprod_dest(A::FilteringMatrix, B::AbstractMatrix, T::Type) = OffsetMatrix(Matrix{T}(undef, size(A, 1), size(B, 2)), axes(A, 1), axes(B, 2))
# NOTE: Filtering matrices other than from a vector signal source should not have an offset in the first dimension <15-05-25> 
LinearAlgebra.matprod_dest(A::Adjoint{<:Any,<:FilteringMatrix{<:Any,1}}, ::AbstractVector, T::Type) = OffsetVector(Vector{T}(undef, size(A, 1)), axes(A, 2))
LinearAlgebra.matprod_dest(A::Adjoint{<:Any,<:FilteringMatrix{<:Any,1}}, B::AbstractMatrix, T::Type) = OffsetMatrix(Matrix{T}(undef, size(A, 1), size(B, 2)), axes(A, 1), axes(B, 2))

offset_mismatch_error(ax1, ax2) = DimensionMismatch("Offsets of the axes do not match in `mul!`. Axes of the arrays must match their offsets but got `$(UnitRange(ax1)) != $(UnitRange(ax2))`.")

function check_offsets(C::AbstractVector, A::AbstractArray, B::AbstractVector)
    axes(A, 2) == axes(B, 1) || throw(offset_mismatch_error(axes(A, 2), axes(B, 1)))
    axes(A, 1) == axes(C, 1) || throw(offset_mismatch_error(axes(A, 1), axes(C, 1)))
end
function check_offsets(C::AbstractArray, A::AbstractArray, B::AbstractArray)
    axes(C, 2) == axes(B, 2) || throw(offset_mismatch_error(axes(C, 2), axes(B, 2)))
    @views check_offsets(C[:, first(axes(C, 2))], A, B[:, first(axes(B, 2))])
end

# NOTE: It is necessary to have two functions for vectors and matrices to avoid ambiguity <15-05-25> 

# STEP 1a: check and remove the offsets -> We are working further on with OffsetArray{<:Any,1 or 2,ParentArray}
function LinearAlgebra.mul!(C::AbstractVector, A::Union{Adjoint{<:Any,<:FilteringMatrix},<:FilteringMatrix}, B::AbstractVector, α::Number, β::Number)
    check_offsets(C, A, B)
    _contract!(no_offset_view(C), no_offset_view(A), no_offset_view(B), α, β)
    return OA.Origin(C)(C)
end

# STEP 1b:  check and remove the offsets -> We are working further on with OffsetMatrix{<:Any, ParentArray}
function LinearAlgebra.mul!(C::AbstractMatrix, A::Union{Adjoint{<:Any,<:FilteringMatrix},<:FilteringMatrix}, B::AbstractMatrix, α::Number, β::Number)
    check_offsets(C, A, B)
    _contract!(no_offset_view(C), no_offset_view(A), no_offset_view(B), α, β)
    return OA.Origin(C)(C)
end

@inline circulant(A::OffsetMatrix{<:Any,<:Adjoint{<:Any,<:FilteringMatrix}}) = circulant(parent(A))
@inline circulant(A::OffsetMatrix{<:Any,<:FilteringMatrix}) = circulant(parent(A))
@inline circulant(A::Adjoint{<:Any,<:FilteringMatrix}) = circulant(parent(A))
@inline circulant(A::FilteringMatrix) = A.circulant

@inline outaxes(A::OffsetMatrix{<:Any,<:Adjoint{<:Any,<:FilteringMatrix}}) = outaxes(parent(A))
@inline outaxes(A::OffsetMatrix{<:Any,<:FilteringMatrix}) = outaxes(parent(A))
@inline outaxes(A::Adjoint{<:Any,<:FilteringMatrix}) = outaxes(parent(A))
@inline outaxes(A::FilteringMatrix) = A.Aaxes

@inline kernaxes(A::OffsetMatrix{<:Any,<:Adjoint{<:Any,<:FilteringMatrix}}) = kernaxes(parent(A))
@inline kernaxes(A::OffsetMatrix{<:Any,<:FilteringMatrix}) = parent(A).Kaxes
@inline kernaxes(A::Adjoint{<:Any,<:FilteringMatrix}) = parent(A).Kaxes
@inline kernaxes(A::FilteringMatrix) = A.Kaxes

const OffsetFilteringMatrix{N} = Union{FM,OffsetMatrix{T,FM}} where {T,FM<:FilteringMatrix{T,N}}
function _contract!(C, A::OffsetFilteringMatrix{1}, B::AbstractMatrix, α, β)
    @tensor begin
        @notensor AC = no_offset_view(circulant(A))
        C[a, b] = α * AC[a, i] * B[i, b] + β * C[a, b]
    end
end
function _contract!(C, A::OffsetFilteringMatrix{1}, B::AbstractVector, α, β)
    @tensor begin
        @notensor AC = no_offset_view(circulant(A))
        C[a] = α * AC[a, i] * B[i] + β * C[a]
    end
end
function _contract!(C, A::OffsetFilteringMatrix{2}, B::AbstractMatrix, α, β)
    @tensor begin
        @notensor BT = no_offset_view(reshape(B, :, kernaxes(A)...))
        @notensor AC = no_offset_view(circulant(A))
        CM[a, b, c] := α * AC[a, b, i, j] * BT[c, i, j]
        @notensor CM = no_offset_view(reshape(CM, axes(C)))
        C[a, b] = CM[a, b] + β * C[a, b]
    end
end
function _contract!(C, A::OffsetFilteringMatrix{2}, B::AbstractVector, α, β)
    @tensor begin
        @notensor BT = no_offset_view(reshape(B, kernaxes(A)...))
        @notensor AC = no_offset_view(circulant(A))
        CM[a, b] := α * AC[a, b, i, j] * BT[i, j]
        @notensor CM = reshape(CM, :)
        C[a] = CM[a] + β * C[a]
    end
end
function _contract!(C, A::OffsetFilteringMatrix{3}, B::AbstractMatrix, α, β)
    @tensor begin
        @notensor BT = no_offset_view(reshape(B, :, A.Kaxes...))
        @notensor AC = no_offset_view(circulant(A))
        CM[a, b, c, d] := α * AC[a, b, c, i, j, k] * BT[d, i, j, k]
        @notensor CM = no_offset_view(reshape(CM, axes(C)))
        C[a, b] = CM[a, b] + β * C[a, b]
    end
end
function _contract!(C, A::OffsetFilteringMatrix{3}, B::AbstractVector, α, β)
    @tensor begin
        @notensor BT = no_offset_view(reshape(B, kernaxes(A)...))
        @notensor AC = no_offset_view(circulant(A))
        CM[a, b, c] := α * AC[a, b, c, i, j, k] * BT[i, j, k]
        @notensor CM = reshape(CM, :)
        C[a] = CM[a] + β * C[a]
    end
end

const OffsetAdjointFilteringMatrix{N} = Union{Adjoint{T,FM},OffsetMatrix{T,Adjoint{T,FM}}} where {T,FM<:FilteringMatrix{T,N}}
function _contract!(C, A::OffsetAdjointFilteringMatrix{1}, B::AbstractMatrix, α, β)
    @tensor begin
        @notensor AC = no_offset_view(circulant(A))
        C[a, b] = α * AC[i, a] * B[i, b] + β * C[a, b]
    end
end
function _contract!(C, A::OffsetAdjointFilteringMatrix{1}, B::AbstractVector, α, β)
    @tensor begin
        @notensor AC = no_offset_view(circulant(A))
        C[a] = α * AC[i, a] * B[i] + β * C[a]
    end
end
function _contract!(C, A::OffsetAdjointFilteringMatrix{2}, B::AbstractMatrix, α, β)
    @tensor begin
        @notensor BT = no_offset_view(reshape(B, :, outaxes(A)...))
        @notensor AC = no_offset_view(circulant(A))
        CM[a, b, c] := α * AC[j, i, b, a] * BT[c, i, j]
        @notensor CM = no_offset_view(reshape(CM, axes(C)))
        C[a, b] = CM[a, b] + β * C[a, b]
    end
end
function _contract!(C, A::OffsetAdjointFilteringMatrix{2}, B::AbstractVector, α, β)
    @tensor begin
        @notensor BT = no_offset_view(reshape(B, outaxes(A)...))
        @notensor AC = no_offset_view(circulant(A))
        CM[a, b] := α * AC[i, j, a, b] * BT[i, j]
        @notensor CM = reshape(CM, :)
        C[a] = CM[a] + β * C[a]
    end
end
function _contract!(C, A::OffsetAdjointFilteringMatrix{3}, B::AbstractMatrix, α, β)
    @tensor begin
        @notensor BT = no_offset_view(reshape(B, :, outaxes(A)...))
        @notensor AC = no_offset_view(circulant(A))
        CM[a, b, c, d] := α * AC[i, j, k, a, b, c] * BT[d, i, j, k]
        @notensor CM = no_offset_view(reshape(CM, axes(C)))
        C[a, b] = CM[a, b] + β * C[a, b]
    end
end
function _contract!(C, A::OffsetAdjointFilteringMatrix{3}, B::AbstractVector, α, β)
    @tensor begin
        @notensor BT = no_offset_view(reshape(B, outaxes(A)...))
        @notensor AC = no_offset_view(circulant(A))
        CM[a, b, c] := α * AC[i, j, k, a, b, c] * BT[i, j, k]
        @notensor CM = reshape(CM, :)
        C[a] = CM[a] + β * C[a]
    end
end

# FIX: It is not always a block matrix... It has blocks only if the sizes check out. The last block is not guaranteed to
# be a Toeplitz matrix <05-05-25> 
BlockArrays.blockaxes(A::FilteringMatrix) = (BlockRange(axes(A.circulant, 1)), BlockRange(axes(A.circulant, 2)))
Base.getindex(A::FilteringMatrix, I::Block{1}) = error("TODO")

# WARN: Do not commit!!! Due to LanguageServer.jl error. <05-05-25> 
public CirculantTensor, FilteringMatrix, OuterInnerArray
