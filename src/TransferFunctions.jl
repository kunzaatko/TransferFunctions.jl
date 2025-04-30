"""
Package for models, estimation, sampling and deconvolution with microscopy transfer functions.
"""
module TransferFunctions
using SpecialFunctions, FillArrays, FFTW, LazyGrids, Interpolations, Roots, IntervalSets, ColorTypes,
    ImageFiltering, OffsetArrays
using ImageCore, StaticArraysCore
using OffsetArrays: centered, center
using OffsetArrays: OffsetArrays as OA

using Reexport
@reexport using Unitful

include("types.jl")
include("utils.jl")
include("interfaces.jl")
include("linear-transfer-functions.jl")

# utils
include("apodization.jl")

export restore, transfer
export SpatialArray, SampledArray

end
