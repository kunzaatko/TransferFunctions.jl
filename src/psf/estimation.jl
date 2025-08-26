module Estimation
using Base: Indices
using InterfaceFunctions
using TransferFunctions
using TransferFunctions: Size, aroundorigin, PSFModel, TransferFunction, PSFArray, sampling, inner_axes, PSFArray
using OffsetArrays: no_offset_view
include("estimation/beads-acquisition.jl")

abstract type EstimationAlgorithm{H<:TransferFunction} end
abstract type FromGroundTruth{H<:TransferFunction} <: EstimationAlgorithm{H} end

# TODO: Add subtypes and output type from estimation to the required interface <06-05-25> 
"""
    estimate(model::Estimation.FromGroundTruth, raw::SpatialArray, gt::SpatialArray)
Estimate the PSF from the ground truth `gt` and transferred data `raw`.
"""
@interface estimate(model::Estimation.FromGroundTruth, raw::AbstractArray, gt::AbstractArray)

# TODO: Should not hold the size but instead the indices that are used for the construction of the `FilteringMatrix` <05-05-25> 
"""
   LeastSquares <: FromGroundTruth{PSFArray}

# Fields
- `size::Size{2}` -- size of the PSF to estimate
"""
struct LeastSquares{N, I<:Indices{N}} <: FromGroundTruth{PSFArray{N}}
    indices::I
    LeastSquares(inds::Indices) = new{length(inds), typeof(inds)}(inds)
end
LeastSquares(size::Size) = LeastSquares(aroundorigin(size))

function estimate(ls::LeastSquares{N}, raw::SpatialArray{<:Any,N}, gt::SpatialArray{<:Any,N}) where {N}
    @assert sampling(raw) == sampling(gt) "sampling rates of `raw` and `gt` must match"
    F_A = filtering_matrix(reflect(gt), ls.indices)
    B = raw[inner_axes(raw, ls.indices)...]
    LhS_op = F_A * no_offset_view(reflect(F_A'))
    RhS = F_A * B[:]
    h_vec =  LhS_op \ RhS
    return PSFArray(SampledArray(reshape(h_vec, ls.indices), sampling(raw)))
end

export LeastSquares

end
