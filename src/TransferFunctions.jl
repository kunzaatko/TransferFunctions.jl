"""
Package for models, estimation, sampling and deconvolution with microscopy transfer functions.
"""
module TransferFunctions
using SpecialFunctions, FFTW, Roots, IntervalSets, ImageFiltering, OffsetArrays
using OffsetArrays: centered, center, no_offset_view

using Reexport
@reexport using Unitful

include("types.jl")

include("types/sampled-arrays.jl")
include("types/circulant-tensors.jl")
include("types/filtering-matrices.jl")
include("types/border-arrays.jl")
include("types/reflected-arrays.jl")

include("utils.jl")
include("apodization.jl")

include("types/tapered-arrays.jl")

include("transfer-function-interface.jl")

include("linear-transfer-functions.jl")
include("nonlinear-transfer-functions.jl")

include("package-extensions.jl")

end
