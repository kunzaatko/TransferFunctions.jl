using Base: Indices
using Unitful: Length

@derived_dimension Frequency Unitful.𝐋^-1 true

"""
    PixelSize{N}
`N`-dimensional pixel size. Alias for `NTuple{N,Length}`
"""
const PixelSize{N} = NTuple{N,Length}

"""
    Coordinate{N,T}
`N`-dimensional coordinate. Alias for `NTuple{N,T<:Real}`
"""
const Coordinate{N,T} = NTuple{N,T} where {T<:Real}

"""
    Size{N}
`N`-dimensional size of an array. Alias for `NTuple{N,Int}`
"""
const Size{N} = NTuple{N,Int}

include("types/extension-interface.jl")
"""
    OneEdge
Alias for `Tuple{Int,Int}`.

An edge is a tuple with the interpretation of `(start, end)` that determines the padding applied in
a [`BorderArray`](@ref) or a [`TaperedArray`](@ref).

See also [`Edges`](@ref)
"""
const OneEdge = Tuple{Int,Int}

"""
    Edges{N}
Edges of an `N`-dimensional array. Alias for `NTuple{N,OneEdge}`

See also [`OneEdge`](@ref), [`TaperedArray`](@ref), [`BorderArray`](@ref)
"""
const Edges{N} = NTuple{N,OneEdge}
