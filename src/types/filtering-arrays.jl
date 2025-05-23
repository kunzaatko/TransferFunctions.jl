# TODO: Replace the CirculantTensor implementation by a simpler implementation that only calculates the array indices
# during getindex calls <23-05-25> 
# TODO: Test whether it would be faster to directly compute the indices of the array based on the kernel size and the
# array size instead of constructing the SubArrays. Benchmark these tests. <16-05-25> 
# TODO: CirculantTensor should store the border type. We could optimize for "reflect" using fft <10-05-25> 
# TODO: CirculantTensor should be parametrized by the types of the axes of the `kern` and `interior`. These should be
# passed as KAX, and IAX similarly as with the `Flattened` array in order to avoid complication <18-05-25> 
# TODO: Update the type parameter documentation <19-05-25> 
"""
    CirculantTensor{P,M,AA,T,N} <: AbstractArray{T,N}

`N`-dimensional circulant tensor of an `M`-dimensional array `A` with kernel indices `kern`. 

The dimensionality `N` is equal to `2M`, where the first `M` dimensions have interior filtering coordinates
`ct.interior` given by the kernel axes `ct.kern` and the axes of `A` and the tail `M` dimensions have the kernel
coordinates `ct.kern`.

A correlation of array `A` with the kernel `K`can be obtained by outer tensor contraction over the tail `M` dimensions
of the `circulant(A,K)` with `K`.

[`parent(ct::CirculantTensor)`](@ref) will return the [`Flattened`](@ref) array of views inside of `A`.

See also [`FilteringMatrix`](@ref), [`imfilter`](@extref), [`circulant`](@ref)
"""
struct CirculantTensor{P,M,AA,T,N} <: AbstractArray{T,N}
    parent::P 
    interior::Indices{M}
    kern::Indices{M}
end

function CirculantTensor(AA::Type{<:AbstractArray}, parent::Flattened, interior::Indices{M}, kern::Indices{M}) where {M}
    return CirculantTensor{typeof(parent),M,AA,eltype(parent),2M}(parent, interior, kern)
end

function basearray(A::CirculantTensor)
    flattened = parent(parent(A))
    return flattened isa OffsetArray ? parent(flattened) : flattened
end

@reexport using ImageFiltering: reflect # NOTE: For convolution instead of correlation <24-04-25> 

function circulant(A::AbstractArray{<:Any, M}, KI::Indices{M}) where {M}
    AI = interior(axes(A), KI)
    any(iszero, length(AI)) && throw(DimensionMismatch("In constructor `CirculantTensor(<:AbstractArray, ::Indices)` the array is not large enough for the kernel. Got interior of $AI.")) # TODO: Test
    views = map(CartesianIndices(AI)) do I
        @inbounds OffsetArray(view(A, CartesianIndices(KI) .+ I), KI)
    end
    parent = flatten(OffsetArray(views, AI); outer=Dims(1:M))
    return CirculantTensor(typeof(A), parent, AI, KI)
end

# Step 1: Determine kernel indices
@inline circulant(A::AbstractArray{<:Any,N}, kern::AbstractArray{<:Any,N}, args...) where {N} = circulant(A, axes(kern), args...)

# Step 2: Determine border and Initialize it if it is not fully specified. (Default Inner() doesn't need action)
@inline function circulant(A::AbstractArray{<:Any,N}, kern::Indices{N}, border::AbstractString, args...) where {N}
    return circulant(A, kern, ImageFiltering.borderinstance(border), args...)
end
@inline function circulant(A::AbstractArray{<:Any,N}, kern::Indices{N}, border::Pad{0}, args...) where {N}
    border = Pad(border.style, kern)
    return circulant(A, kern, border, args...)
end
@inline function circulant(A::AbstractArray{<:Any,N}, kern::Indices{N}, border::Fill{T,0}, args...) where {N,T}
    border = Fill(border.value, kern)
    return circulant(A, kern, border, args...)
end

# Step 3: Apply the border and call the inner constructor
function circulant(A::AbstractArray{T,N}, kern::Indices{N}, border::AbstractBorder, args...) where {T,N}
    A = padarray(T, A, border)
    return circulant(A, kern) # TODO: What should happen with the args... ? <10-05-25> 
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

Base.parent(A::CirculantTensor) = A.parent
Base.size(A::CirculantTensor) = size(parent(A))
Base.getindex(A::CirculantTensor, ind...) = getindex(parent(A), ind...)
Base.axes(A::CirculantTensor, ind...) = axes(parent(A), ind...)
Base.similar(A::CirculantTensor, eltype::Type, dims::Dims) = similar(parent(A), eltype, dims)

# TODO: Create a BlockMatrix interface in for the FilteringMatrix <22-05-25> 
# TODO: Decide what to do with the indices of the array. They will not be linearly spaced if it should make sense.
# Or it could be documented that the indices are not reflecting the actual indices of the reference array. <24-04-25> 
"""
    FilteringMatrix{P,M,CT,T} <: AbstractMatrix{T}

A filtering matrix for the given `M`-dimensional array `A` and kernel indices `kern`. I.e. if `K` is an `M`-dimensional
kernel with indices `kern` and `F` is the corresponding `FilteringMatrix`, then `F * K[:]` is the correlation of (or
convolution if `reflect(K)[:]` where `reflect(K)` has the indices `kern`) of the `K` with `A`.

It is implemented as a stacked view into the [`CirculantTensor`](@ref) of `A` for the kernel indices `kern`.
"""
struct FilteringMatrix{P,M,CT,T} <: AbstractMatrix{T}
    circulant::CT
    parent::P # P<:Flattened
    function FilteringMatrix(C::CirculantTensor{<:Any,M}) where {M}
        rows = map(eachslice(C, dims=Dims((M+1):2M))) do s
            reshape(s, :)  # PERF: `reshape` benchmarks better than `view` for 2D matrices <11-05-25> 
        end
        # FIX: Should be specified by only one of `inner`/`outer` when the method allows it <20-05-25> 
        parent = flatten(view(rows, :); outer=(2,), inner=(1,))
        return new{typeof(parent),M,typeof(C),eltype(C)}(C, parent)
    end
end
FilteringMatrix(args...) = FilteringMatrix(circulant(args...))

Base.parent(A::FilteringMatrix) = A.parent
Base.size(A::FilteringMatrix) = size(parent(A))
Base.axes(A::FilteringMatrix) = axes(parent(A))
Base.getindex(A::FilteringMatrix, ind...) = (@inline; getindex(parent(A), ind...))

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
# LinearAlgebra.matprod_dest(A::FilteringMatrix, ::AbstractVector, T::Type) = OffsetVector(Vector{T}(undef, size(A, 1)), axes(A, 2))
# LinearAlgebra.matprod_dest(A::FilteringMatrix, B::AbstractMatrix, T::Type) = OffsetMatrix(Matrix{T}(undef, size(A, 1), size(B, 2)), axes(A, 1), axes(B, 2))
# NOTE: Filtering matrices other than from a vector signal source should not have an offset in the first dimension <15-05-25> 
# LinearAlgebra.matprod_dest(A::Adjoint{<:Any,<:FilteringMatrix{<:Any,1}}, ::AbstractVector, T::Type) = OffsetVector(Vector{T}(undef, size(A, 1)), axes(A, 2))
# LinearAlgebra.matprod_dest(A::Adjoint{<:Any,<:FilteringMatrix{<:Any,1}}, B::AbstractMatrix, T::Type) = OffsetMatrix(Matrix{T}(undef, size(A, 1), size(B, 2)), axes(A, 1), axes(B, 2))

# TODO: Add a note to the documentation about the necessity of calling `mul!` with `no_offset_view` for a filtering matrix <22-05-25> 

_offset_mismatch_error(ax1, ax2) = DimensionMismatch("Offsets of the axes do not match in for `mul!`. Axes of the arrays must match their offsets but got `$(UnitRange(ax1)) != $(UnitRange(ax2))`.")

function mul!_check_axes(C::AbstractVector, A::AbstractArray, B::AbstractVector)
    axes(A, 1) == axes(C, 1) || throw(_offset_mismatch_error(axes(A, 1), axes(C, 1)))
    mul!_check_axes(A,B)
end

function mul!_check_axes(A::AbstractArray, B::AbstractVector)
    axes(A, 2) == axes(B, 1) || throw(_offset_mismatch_error(axes(A, 2), axes(B, 1)))
end

function mul!_check_axes(C::AbstractArray, A::AbstractArray, B::AbstractArray)
    axes(C, 2) == axes(B, 2) || throw(_offset_mismatch_error(axes(C, 2), axes(B, 2)))
    @views mul!_check_axes(C[:, first(axes(C, 2))], A, B[:, first(axes(B, 2))])
end

# FIX: It is not always a block matrix... It has blocks only if the sizes check out. The last block is not guaranteed to
# be a Toeplitz matrix <05-05-25> 
BlockArrays.blockaxes(A::FilteringMatrix) = (BlockRange(axes(A.circulant, 1)), BlockRange(axes(A.circulant, 2)))
Base.getindex(A::FilteringMatrix, I::Block{1}) = error("TODO")

export circulant
