using InterfaceFunctions
using Base: @propagate_inbounds, showarg

# FIX: The signature of the `border_array` function should be updated to `border_array(parent, padding, border)` instead
# to be consistent with other filtering methods such as `filtering_array` <22-08-25> 

"""
    AbstractBorder{T}
Super type for borders of arrays with the element type `T`

# Implementation
A new border type needs to implement 
```julia
Base.getindex(b::AbstractBorder{T}, A::AbstractArray{T,N}, I::Vararg{Int,N}) where {T,N}
```
where `A` is the parent array of the [`BorderArray`](@ref).
"""
abstract type AbstractBorder{T} end
@interface Base.getindex(b::AbstractBorder{T}, A::AbstractArray{T,N}, I::Vararg{Int,N}) where {T,N}

"""
    validextension(b::AbstractBorder, A::AbstractArray)
Returns the maximum padding extension that is valid for the border `b` and the parent array `A`.

```jldoctest
julia> TF.validextension(TF.Circular(), ones(40,40))
((40, 40), (40, 40))

julia> TF.validextension(TF.Symmetric(), ones(40,40))
((39, 39), (39, 39))

julia> TF.validextension(TF.Fill(0), ones(40,40))
((9223372036854775807, 9223372036854775807), (9223372036854775807, 9223372036854775807))
```
"""
@interface validextension(b::AbstractBorder, A::AbstractArray)
function AbstractBorder{T}(s::Symbol) where {T}
    s == :replicate && return Replicate{T}()
    s == :reflect && return Reflect{T}()
    s == :circular && return Circular{T}()
    s == :fill && return Fill{T}(0)
    s == :symmetric && return Symmetric{T}()
    throw(ArgumentError("Unknown border type $s"))
end
AbstractBorder(s::Symbol) = AbstractBorder{Any}(s)

"""
    Fill{T} <: AbstractBorder{T} 
A border that fills the values with `value::T`
"""
struct Fill{T} <: AbstractBorder{T} 
    value::T
    Fill(value) = new{typeof(value)}(value)
    Fill{T}(value::T) where {T} = new{T}(value)
end
Fill{T}() where {T} = Fill(zero(T)) 
@propagate_inbounds function Base.getindex(b::Fill{T}, A::AbstractArray{T,N}, I::Vararg{Int,N}) where {T,N} 
    if all(I .∈ axes(A)) 
        @inbounds A[I...] 
    else
        b.value
    end
end
validextension(::Fill, A::AbstractArray) = ntuple(_ -> (typemax(Int), typemax(Int)), ndims(A))
Base.convert(::Type{Fill{S}}, f::Fill{T}) where {S,T} = Fill{S}(convert(S, f.value))
Base.convert(::Type{AbstractBorder{S}}, f::Fill{T}) where {S,T} = Base.convert(Fill{S}, f)
Base.showarg(io::IO, b::Fill, toplevel) = print(io, "fill($(b.value))")

"""
    IndexMapBorder{T} <: AbstractBorder{T}
A border that maps its values the values of the parent array with a certain index mapping
"""
abstract type IndexMapBorder{T} <: AbstractBorder{T} end
@interface mapindex(b::IndexMapBorder{T}, A::AbstractArray{T,N}, I::Vararg{Int,N}) where {T,N}
"""
    mapindex(b::IndexMapBorder{T}, A::AbstractArray{T,N}, I::Vararg{Int,N}) where {T,N}
Map the [`BorderArray`](@ref) index to a parent array index with the given scheme.
"""
@propagate_inbounds function Base.getindex(b::IndexMapBorder{T}, A::AbstractArray{T,N}, I::Vararg{Int,N}) where {T,N} 
    I_mapped = mapindex(b, A, I...)
    @boundscheck checkbounds(A, I_mapped...)
    @inbounds A[I_mapped...]
end

"""
    IndependentIndexMapBorder{T} <: IndexMapBorder{T}
A border that maps the indices of the border independently for each dimension
"""
abstract type IndependentIndexMapBorder{T} <: IndexMapBorder{T} end 

"""
    Reflect{T} <: IndexMapBorder{T}
A border that reflects the values of the border around the edges of the parent array

```jldoctest
julia> TransferFunctions.BorderArray(reshape(1:16, (4,4)), TransferFunctions.Reflect, 2)
8×8 border_array(reshape(::UnitRange{Int64}, 4, 4), :Reflect) with eltype Int64 with indices -1:6×-1:6:
 6  2  2  6  10  14  14  10
 5  1  1  5   9  13  13   9
 5  1  1  5   9  13  13   9
 6  2  2  6  10  14  14  10
 7  3  3  7  11  15  15  11
 8  4  4  8  12  16  16  12
 8  4  4  8  12  16  16  12
 7  3  3  7  11  15  15  11
```
"""
struct Reflect{T} <: IndexMapBorder{T} end 
@inline mapindex(::Reflect{T}, A::AbstractArray{T,N}, I::Vararg{Int,N}) where {T,N} = map(I, axes(A)) do ind, ax
    ind ∈ ax && return ind
    ind > last(ax) && return 2last(ax) - ind + 1
    ind < first(ax) && return 2first(ax) - ind - 1
end
validextension(b::Reflect, A::AbstractArray) = map(ax -> (length(ax), length(ax)), axes(A))

"""
    Symmetric{T} <: IndexMapBorder{T}
A border that symmetrically continues the values or the parent array around the edges, i.e. mirrors around the last valid
parent index.

```jldoctest
julia> TransferFunctions.BorderArray(reshape(1:16, (4,4)), TransferFunctions.Symmetric, 2)
8×8 border_array(reshape(::UnitRange{Int64}, 4, 4), :Symmetric) with eltype Int64 with indices -1:6×-1:6:
 11  7  3  7  11  15  11  7
 10  6  2  6  10  14  10  6
  9  5  1  5   9  13   9  5
 10  6  2  6  10  14  10  6
 11  7  3  7  11  15  11  7
 12  8  4  8  12  16  12  8
 11  7  3  7  11  15  11  7
 10  6  2  6  10  14  10  6
```
"""
struct Symmetric{T} <: IndexMapBorder{T} end 
@inline mapindex(::Symmetric{T}, A::AbstractArray{T,N}, I::Vararg{Int,N}) where {T,N} = map(I, axes(A)) do ind, ax
    ind ∈ ax && return ind
    ind > last(ax) && return 2last(ax) - ind
    ind < first(ax) && return 2first(ax) - ind
end
validextension(b::Symmetric, A::AbstractArray) = map(ax -> (length(ax) - 1, length(ax) - 1), axes(A))

"""
    Replicate{T}  <: IndexMapBorder{T}
A border that replicates the closest edge value of the parent array

```jldoctest
julia> TransferFunctions.BorderArray(reshape(1:16, (4,4)), TransferFunctions.Replicate, 2)
8×8 border_array(reshape(::UnitRange{Int64}, 4, 4), :Replicate) with eltype Int64 with indices -1:6×-1:6:
 1  1  1  5   9  13  13  13
 1  1  1  5   9  13  13  13
 1  1  1  5   9  13  13  13
 2  2  2  6  10  14  14  14
 3  3  3  7  11  15  15  15
 4  4  4  8  12  16  16  16
 4  4  4  8  12  16  16  16
 4  4  4  8  12  16  16  16
```
"""
struct Replicate{T} <: IndexMapBorder{T} end 
@inline mapindex(::Replicate{T}, A::AbstractArray{T,N}, I::Vararg{Int,N}) where {T,N} = map(I, axes(A)) do ind, ax
    ind ∈ ax && return ind
    ind > last(ax) && return last(ax)
    ind < first(ax) && return first(ax)
end
validextension(::Replicate, A::AbstractArray) = ntuple(_ -> (typemax(Int), typemax(Int)), ndims(A))

"""
    Circular{T} <: IndexMapBorder{T}
A border that wraps the values around the edges of the parent array

```jldoctest
julia> TransferFunctions.BorderArray(reshape(1:16, (4,4)), TransferFunctions.Circular, 2)
8×8 border_array(reshape(::UnitRange{Int64}, 4, 4), :Circular) with eltype Int64 with indices -1:6×-1:6:
 11  15  3  7  11  15  3  7
 12  16  4  8  12  16  4  8
  9  13  1  5   9  13  1  5
 10  14  2  6  10  14  2  6
 11  15  3  7  11  15  3  7
 12  16  4  8  12  16  4  8
  9  13  1  5   9  13  1  5
 10  14  2  6  10  14  2  6
```
"""
struct Circular{T} <: IndexMapBorder{T} end
@inline mapindex(::Circular{T}, A::AbstractArray{T,N}, I::Vararg{Int,N}) where {T,N} = map(I, axes(A)) do ind, ax
    ind ∈ ax && return ind
    ind > last(ax) && return (first(ax) - 1) + (ind - last(ax))
    ind < first(ax) && return (last(ax) + 1) - (first(ax) - ind)
end
validextension(b::Circular, A::AbstractArray) = map(ax -> (length(ax), length(ax)), axes(A))

const IndexMapBorder_subtypes = (:Reflect, :Symmetric, :Circular, :Replicate)
for subtype in IndexMapBorder_subtypes
    @eval begin
        $subtype() = $subtype{Any}()
        Base.eltype(::$subtype{S}) where S = S
        Base.convert(::Type{AbstractBorder{S}}, a::$subtype) where S = $subtype{S}()
        function Base.showarg(io::IO, b::$subtype, toplevel) 
            print(io, ":", nameof(typeof(b)))
            toplevel && print(io, " with output type ", eltype(b))
        end
    end
end


"""
    InvalidBorderExtent <: Exception
Is thrown when trying to construct a `BorderArray` with padding larger than the maximum padding that would produce valid
values.
"""
struct InvalidBorderExtent <: Exception 
    attempted::Tuple
    valid::Tuple
    border::AbstractBorder
end
Base.showerror(io::IO, e::InvalidBorderExtent) = print(io, "InvalidBorderExtent: padding $(e.attempted) is invalid for border $(e.border). Extension must be within $(e.valid)")

"""
    BorderArray{T,N,AA,AB} <: AbstractArray{T,N}
A border array with a border of type `AB` with the parent array of type `AA`.

It is a light wrapper around the parent array that defines `Base.getindex` to return values with its supplied border
scheme.

```jldoctest
julia> TransferFunctions.BorderArray(reshape(1:36, (6,6)), TransferFunctions.Fill, 2)
10×10 border_array(reshape(::UnitRange{Int64}, 6, 6), fill(0)) with eltype Int64 with indices -1:8×-1:8:
 0  0  0   0   0   0   0   0  0  0
 0  0  0   0   0   0   0   0  0  0
 0  0  1   7  13  19  25  31  0  0
 0  0  2   8  14  20  26  32  0  0
 0  0  3   9  15  21  27  33  0  0
 0  0  4  10  16  22  28  34  0  0
 0  0  5  11  17  23  29  35  0  0
 0  0  6  12  18  24  30  36  0  0
 0  0  0   0   0   0   0   0  0  0
 0  0  0   0   0   0   0   0  0  0
```
"""
struct BorderArray{T,N,AA<:AbstractArray{T,N},AB<:AbstractBorder{T}} <: AbstractArray{T,N}
    parent::AA
    border::AB
    padding::NTuple{N, Tuple{Int,Int}}
    function BorderArray{T}(parent::AbstractArray{T, N}, border::AbstractBorder{T}, padding::NTuple{N, Tuple{Int, Int}}) where {T, N} 
        padding_valid = validextension(border, parent)
        all(all(v .>= p) for (v,p) in zip(padding_valid,padding)) || throw(InvalidBorderExtent(padding,padding_valid,border))
        new{T,ndims(parent),typeof(parent),typeof(border)}(parent, border, padding)
    end
end
BorderArray(parent::AbstractArray, border::Symbol, padding) = BorderArray(parent, AbstractBorder(border), padding)
BorderArray(parent::AbstractArray{T}, border::Type{Fill}, padding) where {T} = BorderArray{T}(parent, Fill{T}(), padding)
BorderArray(parent::AbstractArray, border::Type{<:AbstractBorder}, padding) = BorderArray(parent, border(), padding)
BorderArray(parent::AbstractArray{S}, border::AbstractBorder{T}, padding) where {S,T} = BorderArray{S}(parent, convert(AbstractBorder{S}, border), padding)

# TODO: Test `inferpadding` separately. It had mistakes before <22-08-25> 
inferpadding(::AbstractArray{<:Any, N}, padding::Int) where {N} = ntuple(_->(padding,padding), N)
inferpadding(::AbstractArray{<:Any, N}, padding::NTuple{N, Int}) where {N} = Tuple((p,p) for p in padding)
inferpadding(::AbstractArray{<:Any, N}, inds::Indices{N}) where {N} = Tuple(abs.(extrema(i)) for i in inds)

BorderArray{T}(parent::AbstractArray{T}, border::AbstractBorder{T}, padding) where {T}  = BorderArray{T}(parent, border, inferpadding(parent, padding))

Base.parent(a::BorderArray) = (@inline; a.parent)
Base.axes(a::BorderArray) = map(axes(parent(a)), a.padding) do ax, (padl, padr)
    UnitRange(first(ax) - padl, last(ax) + padr)
end
Base.size(a::BorderArray) = (@inline; length.(axes(a)))
@propagate_inbounds function Base.getindex(a::BorderArray{<:Any,N}, I::Vararg{Int,N}) where {N}
    @boundscheck checkbounds(a, I...)
    @inbounds getindex(a.border, parent(a), I...)
end

"""
    border_array(A, border, padding)
Construct a [`BorderArray`](@ref) of `A` with the border `border` and padding `padding`.

```jldoctest
julia> border_array(reshape(1:9, (3,3)), :circular, 2)
7×7 border_array(reshape(::UnitRange{Int64}, 3, 3), :Circular) with eltype Int64 with indices -1:5×-1:5:
 5  8  2  5  8  2  5
 6  9  3  6  9  3  6
 4  7  1  4  7  1  4
 5  8  2  5  8  2  5
 6  9  3  6  9  3  6
 4  7  1  4  7  1  4
 5  8  2  5  8  2  5
```
"""
border_array(parent, border, padding) = BorderArray(parent, border, padding)

function Base.showarg(io::IO, A::BorderArray, toplevel)
    print(io, "border_array(")
    showarg(io, parent(A), false)
    print(io, ", ")
    showarg(io, A.border, false)
    print(io, ")")
    toplevel && print(io, " with eltype ", eltype(A))
end

export border_array
