module Estimation
using InterfaceFunctions
using TransferFunctions: Size, FilteringMatrix, aroundorigin, ModelPSF, TransferFunction, PSFArray
include("estimation/beads-acquisition.jl")

abstract type EstimationAlgorithm{H<:TransferFunction} end
abstract type FromGroundTruth{H<:TransferFunction} <: EstimationAlgorithm{H} end

# TODO: Add subtypes and output type from estimation to the required interface <06-05-25> 
@interface estimate(model::Estimation.FromGroundTruth, raw::AbstractArray, gt::AbstractArray)

# TODO: Should not hold the size but instead the indices that are used for the construction of the `FilteringMatrix` <05-05-25> 
"""
   LeastSquares <: FromGroundTruth{PSFArray}

# Fields
- `size::Size{2}` -- size of the PSF to estimate
"""
struct LeastSquares <: FromGroundTruth{PSFArray}
    size::Size{2}
end

# TODO: Define for SpatialArrays and check the sampling rates of `gt` and `raw` <05-05-25> 
# TODO: Add padding options <05-05-25> 
function estimate(ls::LeastSquares, raw::AbstractMatrix{T}, gt::AbstractMatrix{T}) where {T<:Real}
    A = FilteringMatrix(gt, aroundorigin(ls.size))
    h_vec = (A' * raw[:]) \ (A * A')
    return reshape(h_vec, A.Kaxes)
end

abstract type Optimization end

struct ForwardPassInversion{M<:ModelPSF,O<:Optimization} <: FromGroundTruth{M}
end

abstract type BlindDeconvolution{H<:TransferFunction} <: EstimationAlgorithm{H} end

end
