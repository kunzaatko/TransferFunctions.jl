using Base: DimsInteger, Indices, tail, OneTo
using Base: @propagate_inbounds

# PERF: Could hold a fast multiplicative inverse for better performance since we are using the same division every time
# in the indexing https://github.com/JuliaLang/julia/blob/master/base/multinverses.jl <10-08-25> 
# TODO: Create a BlockMatrix interface in for the FilteringMatrix <22-05-25> 
# IDEA!: An implementation of the BlockArray interface is only when the array has a certain structure size and a border.
# But there could be a type that envelopes the array of any size in the size desired for the full block structure which
# can use the numerical optimizations that lead from use of FFT in the matrix multiplication. A function could be used
# to convert the enveloped and the normal simple Array interface version.

"""
    FilteringMatrix{T,M,AA,KI,II} <: AbstractMatrix{T}
A filtering matrix for the given `M`-dimensional array `A` and kernel indices `kern`. I.e. if `K` is an `M`-dimensional
kernel with indices `kern` and `F` is the corresponding `FilteringMatrix`, then `F * K[:]` is the correlation of (or
convolution if `reflect(K)[:]` where `reflect(K)` has the indices `kern`) of the `K` with `A`.
"""
struct FilteringMatrix{T,M,AA<:AbstractArray{<:Any,M},KI<:DimsInteger{M},II<:DimsInteger{M}} <: AbstractMatrix{T}
    parent::AA
    kern::KI
    interior::II
    function FilteringMatrix(parent::AbstractArray, kern::DimsInteger)
        interior = size(parent) .- kern .+ 1
        return new{eltype(parent),length(kern),typeof(parent),typeof(kern),typeof(interior)}(parent, kern, interior)
    end
end
FilteringMatrix(parent::AbstractArray, kern::Indices) = FilteringMatrix(parent, length.(kern))
FilteringMatrix(parent::AbstractArray, kern::AbstractArray) = FilteringMatrix(parent, axes(kern))

@inline Base.parent(A::FilteringMatrix) = A.parent
Base.size(A::FilteringMatrix) = (prod(A.kern), prod(A.interior))

# NOTE: A similar conversion from the linear index to Cartesian index is done in `abstractarray.jl` with the `_ind2sub`
# function <09-08-25> 
"""
    ind2sub(inds::DimsInteger, ind::Integer)
Convert a linear index `ind` to a subscript index for the array size `inds`.
"""
ind2sub(inds::DimsInteger, ind::Integer) = (@inline; ind2sub_recurse(inds, ind - 1))
function ind2sub_recurse(indslast::NTuple{1}, ind)
    @inline
    (ind + 1,)
end
function ind2sub_recurse(inds, ind)
    @inline
    d = inds[1]
    indnext, indfirst, indlast = div(ind, d), 1, d
    (ind - indlast * indnext + indfirst, ind2sub_recurse(tail(inds), indnext)...)
end

function Base.getindex(A::FilteringMatrix, I::Vararg{Int,2})
    @boundscheck checkbounds(A, I...)
    i, j = I
    KI = ind2sub(A.kern, i) .- 1
    II = ind2sub(A.interior, j) .- 1
    PI = first.(axes(A.parent)) .+ KI .+ II
    return @inbounds A.parent[PI...]
end

# TODO: Instead of filtering_matrix, should supply `corrmatrix` and `convmatrix` functions for the corresponding
# operations <12-08-25> 

"""
    filtering_matrix(A, K, [border])
Construct a [`FilteringMatrix`](@ref) of `A` with the kernel `K` such that the `F'*K[:]` produces a filtered vector of `A`.

If border is specified, the array is padded with the strategy [`border`](@ref AbstractBorder) so that the full extent of
the array `A` is kept in the filtering output.

```jldoctest
julia> filtering_matrix(reshape(1:25, (5,5)), (-1:1, -1:1));

julia> filtering_matrix(reshape(1:25, (5,5)), OAs.OffsetArray(ones(3,3), -1:1, -1:1));

julia> filtering_matrix(reshape(1:9, (3,3)), (-1:1, -1:1), :symmetric)
9×9 filtering_matrix(border_array(reshape(::UnitRange{Int64}, 3, 3), :Symmetric), (3, 3)) with eltype Int64:
 5  4  5  2  1  2  5  4  5
 4  5  6  1  2  3  4  5  6
 5  6  5  2  3  2  5  6  5
 2  1  2  5  4  5  8  7  8
 1  2  3  4  5  6  7  8  9
 2  3  2  5  6  5  8  9  8
 5  4  5  8  7  8  5  4  5
 4  5  6  7  8  9  4  5  6
 5  6  5  8  9  8  5  6  5
```
"""
filtering_matrix(A::AbstractArray, kern) = FilteringMatrix(A, kern)
function filtering_matrix(A::AbstractArray, kern, border)
    padding = kern_padding(kern)
    PA = BorderArray(A, border, padding)
    return FilteringMatrix(PA, kern)
end

function Base.showarg(io::IO, A::FilteringMatrix, toplevel)
    print(io, "filtering_matrix(")
    showarg(io, parent(A), false)
    print(io, ", ")
    print(io, A.kern)
    print(io, ")")
    toplevel && print(io, " with eltype ", eltype(A))
end

export filtering_matrix
