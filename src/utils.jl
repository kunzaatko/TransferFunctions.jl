using TransferFunctions: Length, PixelSize, Coordinate
using Base: OneTo, Indices
using OffsetArrays, Unitful

"""
    fillsize(Δ::Length, N::Integer) => PixelSize{N}
Construct a `PixelSize` with the same sampling `Δ` in every direction.
"""
fillsize(Δ::Length, n::Integer)::PixelSize = (@inline; ntuple(_ -> Δ, Val(n)))

"""
    roundcenter(r::RoundingMode, A)
Calculate the center index of `A` (`::AbstractArray`/`::Indices`/`::Size`/`::AbstractUnitRange`) rounded with the mode `r::RoundingMode`.

See also [`roundupcenter`](@ref), [`rounddowncenter`](@ref), [`exactcenter`](@ref)
"""
roundcenter(r::RoundingMode, A::AbstractArray) = roundcenter(r, axes(A))
roundcenter(r::RoundingMode, inds::Indices) = CartesianIndex(Tuple(roundcenter.(Ref(r), inds)))
roundcenter(r::RoundingMode, d::Size) = roundcenter(r, map(OneTo, d))
roundcenter(r::RoundingMode, ind::AbstractUnitRange{T}) where {T} = first(ind) + round(T, (last(ind) - first(ind)) / 2, r)

"""
    roundupcenter(A)
Calculate center index of `A` rounded-up (i.e. [`fft`](@extref `AbstractFFTs.fft`) center).

# Examples
```jldoctest
julia> TF.roundupcenter(ones(15,15))
CartesianIndex(8, 8)

julia> TF.roundupcenter((16,16))
CartesianIndex(9, 9)
```
"""
roundupcenter(args...) = (@inline; roundcenter(RoundUp, args...))

"""
    rounddowncenter(A)
Calculate center index of `A` rounded-down (i.e. [`ifft`](@extref `AbstractFFTs.ifft`) center).

# Examples
```jldoctest
julia> TF.rounddowncenter(ones(15,15))
CartesianIndex(8, 8)

julia> TF.rounddowncenter((16,16))
CartesianIndex(8, 8)
```
"""
rounddowncenter(args...) = (@inline; roundcenter(RoundDown, args...))

"""
     exactcenter(A)

Calculate the exact center coordinate of `A`.

# Examples
```jldoctest
julia> TF.exactcenter(ones(15,15))
(8.0, 8.0)

julia> TF.exactcenter((16,16))
(8.5, 8.5)
```
"""
exactcenter(A::AbstractArray) = exactcenter(axes(A))
exactcenter(inds::Indices) = map(exactcenter, inds)
exactcenter(d::Size) = exactcenter(map(OneTo, d))
exactcenter(ind::AbstractUnitRange) = first(ind) + (last(ind) - first(ind)) / 2

"""
    contained(A, loc::Coordinate{N})
Return `true` if the coordinate `loc` is contained in the axes of array `A`. 
"""
contained(A::AbstractArray{<:Any,N}, loc::Coordinate{N,Int}) where {N} = all(loc .∈ axes(A))

"""
    interior(inds::Indices{N}, kern::Indices{N})
Return 'valid' indices for convolution of an array with indices `inds` with a kernel having the indices `kern`.
"""
interior(inds::Indices{N}, kern::Indices{N}) where {N} = map(interior, inds, kern)
interior(ind::AbstractUnitRange, kern::AbstractUnitRange) = typeof(ind)(intersect(first(ind)-first(kern):last(ind)-last(kern), ind))
interior(ind::OneTo, kern::AbstractUnitRange) = interior(UnitRange(ind), kern)

"""
    aroundorigin(s::Size, origin=(0,0,...))
Return the offset axes of an array of size `s` around the `origin`.
"""
aroundorigin(sz::Size{N}) where {N} = aroundorigin(sz, ntuple(_ -> 0, Val(N)))
aroundorigin(sz::Size{N}, o::Coordinate{N,<:Integer}) where {N} = map(aroundorigin, sz, o)
aroundorigin(inds::Indices{N}) where {N} = aroundorigin(inds, ntuple(_ -> 0, Val(N)))
aroundorigin(inds::Indices{N}, o::Coordinate{N,<:Integer}) where {N} = map(aroundorigin, inds, o)
aroundorigin(s::Integer, o::Integer=0) = aroundorigin(OneTo(s) .- rounddowncenter(OneTo(s)), o)
aroundorigin(s::AbstractUnitRange, o::Integer=0) = s .+ o
aroundorigin(s, o::CartesianIndex) = aroundorigin(s, Tuple(o))


# TODO: Generalize for n dims and define a single function for this <18-08-25> 
@inline togrid(axes::NTuple{2}) = ([x for x in axes[1], _ in axes[2]], [y for _ in axes[1], y in axes[2]])
@inline function togrid(axes::NTuple{3})
    xs = [x for x in axes[1], _ in axes[2], _ in axes[3]]
    ys = [y for _ in axes[1], y in axes[2], _ in axes[3]]
    zs = [z for _ in axes[1], _ in axes[2], z in axes[3]]
    return xs, ys, zs
end

# TODO: Use the same calling stack as in the previous methods. <05-05-25> 
"""
    freqgrid(s::Size, Δ)
Generate a frequency grid for multidimensional FFTs of the signal of size `s` with the pixel size `Δ` (i.e. sampling
rate `1/Δ`)

```jldoctest
julia> x_f, y_f = TF.freqgrid((5,5), 50u"nm");

julia> x_f
5×5 Matrix{Quantity{Float64, 𝐋^-1, Unitful.FreeUnits{(nm^-1,), 𝐋^-1, nothing}}}:
  0.0 nm^-1     0.0 nm^-1     0.0 nm^-1     0.0 nm^-1     0.0 nm^-1
  0.004 nm^-1   0.004 nm^-1   0.004 nm^-1   0.004 nm^-1   0.004 nm^-1
  0.008 nm^-1   0.008 nm^-1   0.008 nm^-1   0.008 nm^-1   0.008 nm^-1
 -0.008 nm^-1  -0.008 nm^-1  -0.008 nm^-1  -0.008 nm^-1  -0.008 nm^-1
 -0.004 nm^-1  -0.004 nm^-1  -0.004 nm^-1  -0.004 nm^-1  -0.004 nm^-1

julia> y_f' == x_f
true

julia> _,_,z_f = TF.freqgrid((4,4,3), (50u"nm", 30u"nm", 15u"nm"));

julia> z_f
4×4×3 Array{Quantity{Float64, 𝐋^-1, Unitful.FreeUnits{(nm^-1,), 𝐋^-1, nothing}}, 3}:
[:, :, 1] =
 0.0 nm^-1  0.0 nm^-1  0.0 nm^-1  0.0 nm^-1
 0.0 nm^-1  0.0 nm^-1  0.0 nm^-1  0.0 nm^-1
 0.0 nm^-1  0.0 nm^-1  0.0 nm^-1  0.0 nm^-1
 0.0 nm^-1  0.0 nm^-1  0.0 nm^-1  0.0 nm^-1

[:, :, 2] =
 0.0222222 nm^-1  0.0222222 nm^-1  0.0222222 nm^-1  0.0222222 nm^-1
 0.0222222 nm^-1  0.0222222 nm^-1  0.0222222 nm^-1  0.0222222 nm^-1
 0.0222222 nm^-1  0.0222222 nm^-1  0.0222222 nm^-1  0.0222222 nm^-1
 0.0222222 nm^-1  0.0222222 nm^-1  0.0222222 nm^-1  0.0222222 nm^-1

[:, :, 3] =
 -0.0222222 nm^-1  -0.0222222 nm^-1  -0.0222222 nm^-1  -0.0222222 nm^-1
 -0.0222222 nm^-1  -0.0222222 nm^-1  -0.0222222 nm^-1  -0.0222222 nm^-1
 -0.0222222 nm^-1  -0.0222222 nm^-1  -0.0222222 nm^-1  -0.0222222 nm^-1
 -0.0222222 nm^-1  -0.0222222 nm^-1  -0.0222222 nm^-1  -0.0222222 nm^-1
```
"""
@inline freqgrid(sz::Size{N}, Δ::PixelSize{N}) where {N} = togrid(fftfreq.(sz, 1 ./ Δ))
freqgrid(sz::Size{N}, Δ::Length) where {N} = freqgrid(sz, fillsize(Δ, N))

@enum CellPosition left mid right

@inline posaxes(axes::Indices{N}, Δ::PixelSize{N}, c::CellPosition=left) where {N} = posaxes(axes, Δ, Val(c))
@inline posaxes(axes::Indices{N}, Δ::PixelSize{N}, ::Val{left}) where {N} =
    map(axes, Δ) do ax, Δax
        ax .* Δax
    end
@inline posaxes(axes::Indices{N}, Δ::PixelSize{N}, ::Val{right}) where {N} =
    map(axes, Δ) do ax, Δax
        (ax .+ step(ax)) .* Δax
    end
@inline posaxes(axes::Indices{N}, Δ::PixelSize{N}, ::Val{mid}) where {N} =
    map(axes, Δ) do ax, Δax
        (2ax .+ step(ax)) ./ 2 .* Δax
    end

const SizeSpec{N} = Union{Size{N},Indices{N}}
@inline posaxes(sz::Size{N}, Δ::PixelSize{N}, args...; center=roundupcenter(sz)) where {N} = posaxes(Tuple(OneTo(s) .- c for (s, c) in zip(sz, Tuple(center))), Δ, args...)
@inline posaxes(sz::SizeSpec{N}, Δ::Length, args...; kwargs...) where {N} = @inline posaxes(sz, fillsize(Δ, N), args...; kwargs...)

"""
    posgrid(a::Indices, Δ)
    posgrid(s::Size, Δ; center=roundupcenter(s))
Generate a position grid for sampled images with the pixel/voxel size `Δ`.

If a size `s` is passed generate a position grid with the given size and the center in `center`. If `Δ` is a tuple of
[`Length`s](@extref Unitful `Length`) then the elements are used for the sampling in the respective dimensions. If `Δ`
is a single length, then it is used for all dimensions.


  ```jldoctest
julia> x_p, y_p = TF.posgrid((-3:5, -1:3), (50u"nm", 20u"nm"));

julia> x_p
9×5 Matrix{Quantity{Int64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}:
 -150 nm  -150 nm  -150 nm  -150 nm  -150 nm
 -100 nm  -100 nm  -100 nm  -100 nm  -100 nm
  -50 nm   -50 nm   -50 nm   -50 nm   -50 nm
    0 nm     0 nm     0 nm     0 nm     0 nm
   50 nm    50 nm    50 nm    50 nm    50 nm
  100 nm   100 nm   100 nm   100 nm   100 nm
  150 nm   150 nm   150 nm   150 nm   150 nm
  200 nm   200 nm   200 nm   200 nm   200 nm
  250 nm   250 nm   250 nm   250 nm   250 nm

julia> y_p
9×5 Matrix{Quantity{Int64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}:
 -20 nm  0 nm  20 nm  40 nm  60 nm
 -20 nm  0 nm  20 nm  40 nm  60 nm
 -20 nm  0 nm  20 nm  40 nm  60 nm
 -20 nm  0 nm  20 nm  40 nm  60 nm
 -20 nm  0 nm  20 nm  40 nm  60 nm
 -20 nm  0 nm  20 nm  40 nm  60 nm
 -20 nm  0 nm  20 nm  40 nm  60 nm
 -20 nm  0 nm  20 nm  40 nm  60 nm
 -20 nm  0 nm  20 nm  40 nm  60 nm

julia> TF.posgrid((3,3,3), 50u"nm")[3]
3×3×3 Array{Quantity{Int64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}, 3}:
[:, :, 1] =
 -50 nm  -50 nm  -50 nm
 -50 nm  -50 nm  -50 nm
 -50 nm  -50 nm  -50 nm

[:, :, 2] =
 0 nm  0 nm  0 nm
 0 nm  0 nm  0 nm
 0 nm  0 nm  0 nm

[:, :, 3] =
 50 nm  50 nm  50 nm
 50 nm  50 nm  50 nm
 50 nm  50 nm  50 nm
```
"""
posgrid(args...; kwargs...) = togrid(posaxes(args...; kwargs...))

intervalaxes(args...; kwargs...) =
    map(posaxes(args..., Val(left); kwargs...), posaxes(args..., Val(right); kwargs...)) do startax, endax
        map(startax, endax) do s, e
            s .. e
        end
    end
intervalgrid(args...; kwargs...) = togrid(intervalaxes(args...; kwargs...))

## OffsetArray helpers ##

# FIX: Add external link when `OffsetArrays` have `objects.inv` <18-09-25> 
"""
    OriginAt{N}
A helper for `OffsetArrays` which works in a similar (partly inverse) way of `OffsetArrays.Origin` and sets the origin of the argument to the index `origin`. 


```jldoctest
julia> TF.OriginAt((3,3))
TransferFunctions.OriginAt{2}(CartesianIndex(3, 3))

julia> TF.OriginAt(3, 3)
TransferFunctions.OriginAt{2}(CartesianIndex(3, 3))

julia> TF.OriginAt(4)
TransferFunctions.OriginAt{1}(CartesianIndex(4,))

julia> TF.OriginAt(CartesianIndex(-1, 1, 5))
TransferFunctions.OriginAt{3}(CartesianIndex(-1, 1, 5))
```
"""
struct OriginAt{N}
    origin::CartesianIndex{N}
end
OriginAt(ind...) = OriginAt(CartesianIndex(ind...))

# FIX: Extref when `OffsetArrays` have `objects.inv` <18-09-25> 
"""
    (OAt::OriginAt{N})(A::AbstractArray{<:Any,N})
Returns an `OffsetArray` with the parent `A` such that the origin `(0,...,0)` is at the `OAt.origin`.

```jldoctest
julia> TF.OriginAt(3, 3)(reshape(1:16, (4,4)))
4×4 OffsetArray(reshape(::UnitRange{Int64}, 4, 4), -2:1, -2:1) with eltype Int64 with indices -2:1×-2:1:
 1  5   9  13
 2  6  10  14
 3  7  11  15
 4  8  12  16
```
"""
(oat::OriginAt{N})(x::AbstractArray{<:Any,N}) where {N} = OffsetArrays.Origin(CartesianIndex{N}(ntuple(_ -> 1, Val(N))) - oat.origin)(x)

## Filtering helpers ##

"""
    kern_padding(kern)
Determine the padding necessary to keep the input array fully contained in the interior of the output when filtered with
`kern`.
"""
function kern_padding(K::Indices)
    if !all(I -> 0 ∈ I, K)
        @warn "A kernel not containing the origin may lead to unexpected filtering output sizes"
    end
    Tuple((max(0, abs(first(k))), max(0, abs(last(k)))) for k in K)
end
kern_padding(K::AbstractArray) = kern_padding(axes(K))

"""
    inner_axes(A, edges)
    inner_axes(A, K)
Determine the inner axes of the array with edges `edges` or when filtered with kernel `K`.

```jldoctest
julia> TF.inner_axes(ones(100,100), ((2,4), (1,10)))
(3:96, 2:90)

julia> TF.inner_axes(ones(100,100), OAs.OffsetArray(ones(11,11), -5:5, -3:7))
(6:95, 4:93)
```
"""
inner_axes(A::AbstractArray{<:Any,N}, edges::Edges{N}) where {N} = map((a, e) -> (first(a)+e[1]):(last(a)-e[2]), axes(A), edges)
inner_axes(A::AbstractArray, K) = inner_axes(A, kern_padding(K))

## Parameter Checking ##

check_emission_wavelength(λ) = λ > zero(λ) || throw(DomainError(λ, "Emission wavelength is a positive value. Got `λ = $λ`."))
check_numerical_aperture(NA) = NA > zero(NA) || throw(DomainError(NA, "Numerical aperture of the objective is a positive value. Got `NA = $NA`."))
check_refractive_index(n) = n > zero(n) || throw(DomainError(n, "Refractive index of the immersion is a positive value. Got `n = $n`."))

## Printing helpers ##

rounded(params::Vararg{Tuple{String,Any}}; kwargs...) = join([rounded(name, val; kwargs...) for (name, val) in params], ", ")
rounded(name::String, val; kwargs...) = name * "=" * rounded(val; kwargs...)
rounded(val::Quantity; kwargs...) = string(round(unit(val), val; kwargs...))
rounded(val; sigdigits=3, kwargs...) = string(round(val; sigdigits=sigdigits, kwargs...))
rounded(vals::AbstractVector; kwargs...) = "[" * join([rounded(val; kwargs...) for val in vals], ", ") * "]"
