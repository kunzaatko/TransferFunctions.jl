"""
Package for models, estimation, sampling and deconvolution with microscopy transfer functions.
"""
module TransferFunctions
using SpecialFunctions, FillArrays, FFTW, LazyGrids, Interpolations, Roots, IntervalSets, ColorTypes,
    ImageFiltering, OffsetArrays
using ImageCore
using OffsetArrays: centered, center, Origin
using Reexport
using Base: Indices

@reexport using Unitful

using Unitful: Length
@derived_dimension Frequency Unitful.𝐋^-1 true
const PixelSize{N} = NTuple{N,Length}
const Coordinate{N} = NTuple{N,Real}

include("utils.jl")
include("types.jl")
include("linear-transfer-functions.jl")

# utils
include("apodization.jl")

export psf, attenuation, intensity, restore, transfer
export otf, cutoff
export BornWolf, CircularPupilOTF
export SpatialArray

export PSFArray, OTFArray

end
