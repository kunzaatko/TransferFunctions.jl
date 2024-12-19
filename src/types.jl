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
@doc raw"""
# Implementation
To create a new OTF model `A <: ModelOTF`, you need to implement the functions:
+ The transfer coefficient at a give frequency:
    + `attenuation(model::A, kx::Frequency, ky::Frequency)` for a non-symmetric OTF, or
    + `attenuation(model::A, kᵣ::Frequency)` for a symmetric OTF, where you also have to call `@traitimpl
        RadiallySymmetric{A}`
You should consider implementing the following optional methods:
+ `preferred_type(::Type{<:ModelOTF})` method, which returns the natural type that is returned from calling the model.
    It should be either `<:Complex` or `<:Real`.
For a symmetric OTF, you should consider implementing the optional methods:
+  `cutoff(model::A; a)` which returns the largest frequency `f::Frequency` such that `attenuation(model, f) >= a` if 
    `a > 0` and `attenuation(model, f) > 0` if `a == 0`

"""
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


# TODO: Implement a Base.show() in one method for types ModelPSF, ModelOTF using `typeof` and `nameof` <19-12-24> 
