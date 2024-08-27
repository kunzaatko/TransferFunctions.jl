using Reexport
@reexport using Unitful
using Unitful: Length
@derived_dimension Frequency Unitful.𝐋^-1

# TODO: Documentation <28-11-23> 
@doc raw"""
An `N`-dimensional transfer function
"""
abstract type TransferFunction{N} end
const TF{N} = TransferFunction{N}

abstract type OpticalTransferFunction{N} <: TransferFunction{N} end
abstract type ModelOTF{N} <: OpticalTransferFunction{N} end

abstract type PointSpreadFunction{N} <: TransferFunction{N} end
abstract type ModelPSF{N} <: PointSpreadFunction{N} end

@doc raw"""
Union type for transfer function that is based on a physical model of an optical system. Contrary to
a [`MeasuredTransferFunction`](@ref), a `ModelTransferFunction` must be quantifiable at any point (in either spatial or
frequency domain).
"""
const ModelTransferFunction{N} = Union{ModelPSF{N},ModelOTF{N}}

const PixelSize{N} = NTuple{N,Length}
fillsize(Δxy::Length, n::Int)::PixelSize = fill(Δxy, n) |> Tuple

const Coordinate{N} = NTuple{N,Real}

@doc """
If the pupil function of the system is symmetric, the OTF as well as the PSF are radially symmteric which can be used to
optimize the calculations.
"""
@traitdef RadiallySymmetric{TF<:TransferFunction}
