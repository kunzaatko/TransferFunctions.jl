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
    SBitVector{N}
`N`-dimensional static vector with elements of type `Bool`. Alias for `SVector{N,Bool}`
"""
const SBitVector{N} = SVector{N,Bool}

"""
    Size{N}
`N`-dimensional size of an array. Alias for `NTuple{N,Int}`
"""
const Size{N} = NTuple{N,Int}

include("types/arrays.jl")
