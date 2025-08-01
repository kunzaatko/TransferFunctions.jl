"""
Package for models, estimation, sampling and deconvolution with microscopy transfer functions.
"""
module TransferFunctions
using SpecialFunctions, FFTW, Interpolations, Roots, IntervalSets, ImageFiltering, OffsetArrays
using OffsetArrays: centered, center, no_offset_view
using OffsetArrays: OffsetArrays as OA

using Reexport
@reexport using Unitful

include("types.jl")
include("utils.jl")

# utils
include("apodization.jl")

include("interfaces.jl")
include("linear-transfer-functions.jl")

include("extensions.jl")

export restore, transfer
export SpatialArray, SampledArray

end
