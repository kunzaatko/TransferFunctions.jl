using TransferFunctions: Length, PixelSize, Coordinate
using ImageFiltering: padarray
@reexport using ImageFiltering: Inner, Pad, Fill
using SpecialFunctions

fillsize(Δ::Length, n::Integer)::PixelSize = (@inline; ntuple(_ -> Δ, Val(n)))

roundcenter(r::RoundingMode, A::AbstractArray) = roundcenter(r, axes(A))
roundcenter(r::RoundingMode, inds::Indices) = CartesianIndex(Tuple(roundcenter.(Ref(r), inds)))
roundcenter(r::RoundingMode, d::Size) = roundcenter(r, map(Base.OneTo, d))
roundcenter(r::RoundingMode, ind::AbstractUnitRange{T}) where {T} = first(ind) + round(T, (last(ind) - first(ind)) / 2, r)

"""
    roundupcenter(A)

Calculate center index of `A` rounded-up (i.e. `fft` center).

# Examples
```jldoctest; setup = :(using TransferFunctions: roundupcenter)
julia> roundupcenter(ones(15,15))
CartesianIndex(8, 8)

julia> roundupcenter((16,16))
CartesianIndex(9, 9)
```
"""
roundupcenter(args...) = (@inline; roundcenter(RoundUp, args...))

"""
    rounddowncenter(A)

Calculate center index of `A` rounded-down (i.e. `ifft` center).

# Examples
```jldoctest; setup = :(using TransferFunctions: rounddowncenter)
julia> rounddowncenter(ones(15,15))
CartesianIndex(7, 7)

julia> rounddowncenter((16,16))
CartesianIndex(8, 8)
```
"""
rounddowncenter(args...) = (@inline; roundcenter(RoundDown, args...))

"""
     exactcenter(A)

Calculate the exact center coordinate of `A`.

# Examples
```jldoctest; setup = :(using TransferFunctions: exactcenter)
julia> exactcenter(ones(15,15))
(8.0, 8.0)

julia> exactcenter((16,16))
(8.5, 8.5)
```
"""
exactcenter(A::AbstractArray) = exactcenter(axes(A))
exactcenter(inds::Indices) = map(exactcenter, inds)
exactcenter(d::Size) = exactcenter(map(Base.OneTo, d))
exactcenter(ind::AbstractUnitRange) = first(ind) + (last(ind) - first(ind)) / 2

"""
    contained(A, loc::Coordinate{N})
Return `true` if the coordinate `loc` is contained in the axes of array `A`. 
"""
contained(A::AbstractArray{<:Any,N}, loc::Coordinate{N}) where {N} = all(loc .∈ axes(A))

"""
    interior(inds::Indices{N}, kern::Indices{N})
Return 'valid' indices for convolution of an array with indices `inds` with a kernel having the indices `kern`.
"""
interior(inds::Indices{N}, kern::Indices{N}) where {N} = map(interior, inds, kern)
interior(ind::AbstractUnitRange, kern::AbstractUnitRange) = typeof(ind)(intersect(first(ind)-first(kern):last(ind)-last(kern), ind))
interior(ind::Base.OneTo, kern::AbstractUnitRange) = interior(UnitRange(ind), kern)

"""
    aroundorigin(s::Size, origin=(0,0,...))
Return the offset axes of an array of size `s` around the `origin`.
"""
aroundorigin(sz::Size{N}) where {N} = aroundorigin(sz, ntuple(_ -> 0, Val(N)))
aroundorigin(sz::Size{N}, o::Coordinate{N,<:Integer}) where {N} = map(aroundorigin, sz, o)
aroundorigin(inds::Indices{N}) where {N} = aroundorigin(inds, ntuple(_ -> 0, Val(N)))
aroundorigin(inds::Indices{N}, o::Coordinate{N,<:Integer}) where {N} = map(aroundorigin, inds, o)
aroundorigin(s::Integer, o::Integer=0) = aroundorigin(Base.OneTo(s) .- rounddowncenter(Base.OneTo(s)), o)
aroundorigin(s::AbstractUnitRange, o::Integer=0) = s .+ o

# TODO: Use the same calling stack as in the previous methods. <05-05-25> 
# TODO: Instead of using `ndgrid` here, it should be `inlined` multiplication by `ones` in the size of the required <05-05-25> 
fftfreqs(sz::Size{2}, Δ::PixelSize{2}) = ndgrid(fftfreq.(sz, 1 ./ Δ)...)
fftfreqs(sz::Size{2}, Δ::Length) = fftfreqs(sz, fillsize(Δ, 2))

posgrid(sz::Size{2}, Δ::PixelSize{2}) = ndgrid(map(sz, Tuple(roundupcenter(sz)), Δ) do len, c, samp
    (range(1, len) .- c) .* samp
end...)
posgrid(sz::Size{2}, Δ::Length) = posgrid(sz, fillsize(Δ, 2))

struct OriginAt{N}
    origin::CartesianIndex{N}
end
(oat::OriginAt{N})(x::AbstractArray{<:Any,N}) where {N} = OA.Origin(CartesianIndex{N}(ntuple(_ -> 1, Val(N))) - oat.origin)(x)

# NOTE: https://github.com/JuliaLang/julia/issues/6733
"""
    @__FUNCTION__
Return current function name. For debugging and error message purposes.
"""
macro __FUNCTION__()
    return :($(esc(Expr(:isdefined, :var"#self#"))) ? $(esc(:var"#self#")) : nothing)
end
