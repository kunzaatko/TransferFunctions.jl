using Base: Indices, @propagate_inbounds, split_rest, showarg

# TODO: Update the type parameter documentation <19-05-25> 

"""
    CirculantTensor{T,AA,N,KI,II} <: AbstractArray{T,N}

`N`-dimensional circulant tensor of an `M`-dimensional array `A` with kernel indices `kern`. 

The dimensionality `N` is equal to `2M`, where the first `M` dimensions have interior filtering coordinates
`ct.interior` given by the kernel axes `ct.kern` and the axes of `A` and the tail `M` dimensions have the kernel
coordinates `ct.kern`.

A correlation of array `A` with the kernel `K`can be obtained by outer tensor contraction over the tail `M` dimensions
of the `circulant(A,K)` with `K`.

See also [`FilteringMatrix`](@ref), [`BorderArray`](@ref) [`circulant`](@ref)
"""
struct CirculantTensor{T,AA,N,KI,II} <: AbstractArray{T,N}
    parent::AA
    kern::KI
    interior::II
    function CirculantTensor(parent::AbstractArray, kern::Indices)
        I_interior = interior(axes(parent), kern)
        new{eltype(parent), typeof(parent),2*length(kern),typeof(kern), typeof(I_interior)}(parent,kern,I_interior)
    end
end
CirculantTensor(A::AbstractArray, kern::AbstractArray) = CirculantTensor(A, axes(kern))

@inline Base.parent(A::CirculantTensor) = A.parent
Base.axes(A::CirculantTensor) = (A.kern..., A.interior...)
Base.size(A::CirculantTensor) = length.(axes(A))

@propagate_inbounds function Base.getindex(A::CirculantTensor{<:Any,<:Any,N}, I::Vararg{Int,N})  where {N}
    @boundscheck checkbounds(A, I...)
    KI, II = Base.split_rest(I, length(A.interior))
    @inbounds A.parent[(KI .+ II)...]
end

"""
    circulant(A, K, border)
Construct a [`CirculantTensor`](@ref) of `A` with the kernel `K` such that tensor contraction with the front indices
with `K` outputs the filtered array of `A` with `K`.

If border is specified, the array is padded with the strategy [`border`](@ref AbstractBorder) so that the full extent of
the array `A` is kept in the contraction output.

```jldoctest
julia> circulant(reshape(1:25, (5,5)), (-1:1, -1:1));

julia> circulant(reshape(1:25, (5,5)), OAs.OffsetArray(ones(3,3), -1:1, -1:1));

julia> c = circulant(reshape(1:25, (5,5)), (-1:1, -1:1), :replicate);

julia> size(c)
(3, 3, 5, 5)
```
"""
circulant(A::AbstractArray, kern) =  CirculantTensor(A, kern)
function circulant(A::AbstractArray, kern, border)
    padding = kern_padding(kern)
    PA  = BorderArray(A, border, padding)
    return CirculantTensor(PA, kern)
end

function Base.showarg(io::IO, A::CirculantTensor, toplevel)
    print(io, "circulant(")
    showarg(io, parent(A), false)
    print(io, ", ")
    print(io, A.kern)
    print(io, ")")
    toplevel && print(io, " with eltype ", eltype(A))
end

export circulant
