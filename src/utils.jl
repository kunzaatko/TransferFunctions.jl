using TransferFunctions: Length, PixelSize, Coordinate
using Base: OneTo
using OffsetArrays

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

Calculate center index of `A` rounded-up (i.e. `fft` center).

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

Calculate center index of `A` rounded-down (i.e. `ifft` center).

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
contained(A::AbstractArray{<:Any,N}, loc::Coordinate{N, Int}) where {N} = all(loc .∈ axes(A))

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

# TODO: Use the same calling stack as in the previous methods. <05-05-25> 
"""
    fftfreqs(sz::Size{2}, Δ::PixelSize{2})
    fftfreqs(sz::Size{2}, Δ::Length)
Generate a frequency grid for 2D FFTs with the pixel size `Δ`

```jldoctest
julia> x_f, y_f = TF.fftfreqs((5,5), 50u"nm");

julia> x_f
5×5 Matrix{Quantity{Float64, 𝐋^-1, Unitful.FreeUnits{(nm^-1,), 𝐋^-1, nothing}}}:
  0.0 nm^-1     0.0 nm^-1     0.0 nm^-1     0.0 nm^-1     0.0 nm^-1
  0.004 nm^-1   0.004 nm^-1   0.004 nm^-1   0.004 nm^-1   0.004 nm^-1
  0.008 nm^-1   0.008 nm^-1   0.008 nm^-1   0.008 nm^-1   0.008 nm^-1
 -0.008 nm^-1  -0.008 nm^-1  -0.008 nm^-1  -0.008 nm^-1  -0.008 nm^-1
 -0.004 nm^-1  -0.004 nm^-1  -0.004 nm^-1  -0.004 nm^-1  -0.004 nm^-1

julia> y_f
5×5 Matrix{Quantity{Float64, 𝐋^-1, Unitful.FreeUnits{(nm^-1,), 𝐋^-1, nothing}}}:
 0.0 nm^-1  0.004 nm^-1  0.008 nm^-1  -0.008 nm^-1  -0.004 nm^-1
 0.0 nm^-1  0.004 nm^-1  0.008 nm^-1  -0.008 nm^-1  -0.004 nm^-1
 0.0 nm^-1  0.004 nm^-1  0.008 nm^-1  -0.008 nm^-1  -0.004 nm^-1
 0.0 nm^-1  0.004 nm^-1  0.008 nm^-1  -0.008 nm^-1  -0.004 nm^-1
 0.0 nm^-1  0.004 nm^-1  0.008 nm^-1  -0.008 nm^-1  -0.004 nm^-1
```
"""
@inline fftfreqs(sz::Size{2}, Δ::PixelSize{2}) = (fftfreq(sz[1], 1 / Δ[1]) * ones(sz[2])', ones(sz[1]) * fftfreq(sz[2], 1 / Δ[2])')
fftfreqs(sz::Size{2}, Δ::Length) = fftfreqs(sz, fillsize(Δ, 2))

"""
    posgrid(sz::Size{2}, Δ::PixelSize{2}; center=roundupcenter(sz))
Generate a position grid for 2D sampled images with the pixel size `Δ`

```jldoctest
julia> x_p, y_p = TF.posgrid((5, 5), 50u"nm");

julia> x_p
5×5 Matrix{Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}:
 -100.0 nm  -100.0 nm  -100.0 nm  -100.0 nm  -100.0 nm
  -50.0 nm   -50.0 nm   -50.0 nm   -50.0 nm   -50.0 nm
    0.0 nm     0.0 nm     0.0 nm     0.0 nm     0.0 nm
   50.0 nm    50.0 nm    50.0 nm    50.0 nm    50.0 nm
  100.0 nm   100.0 nm   100.0 nm   100.0 nm   100.0 nm

julia> y_p
5×5 Matrix{Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}:
 -100.0 nm  -50.0 nm  0.0 nm  50.0 nm  100.0 nm
 -100.0 nm  -50.0 nm  0.0 nm  50.0 nm  100.0 nm
 -100.0 nm  -50.0 nm  0.0 nm  50.0 nm  100.0 nm
 -100.0 nm  -50.0 nm  0.0 nm  50.0 nm  100.0 nm
 -100.0 nm  -50.0 nm  0.0 nm  50.0 nm  100.0 nm
```
"""
@inline posgrid(sz::Size{2}, Δ::PixelSize{2}; center=roundupcenter(sz)) = (((OneTo(sz[1]) .- center[1]) * Δ[1]) * ones(sz[2])', ones(sz[1]) * ((OneTo(sz[2]) .- center[2]) * Δ[2])')
posgrid(sz::Size{2}, Δ::Length) = posgrid(sz, fillsize(Δ, 2))

## OffsetArray helpers ##

struct OriginAt{N}
    origin::CartesianIndex{N}
end
(oat::OriginAt{N})(x::AbstractArray{<:Any,N}) where {N} = OffsetArrays.Origin(CartesianIndex{N}(ntuple(_ -> 1, Val(N))) - oat.origin)(x)

## Filtering helpers ##

"""
    kern_padding(kern)
Determine the padding necessary to keep the input array fully contained in the interior of the output when filtered with
`kern`.
"""
kern_padding(kern::Indices) = Tuple((max(0, abs(first(k))),max(0, abs(last(k)))) for k in kern)
kern_padding(kern::AbstractArray) = kern_padding(axes(kern))
