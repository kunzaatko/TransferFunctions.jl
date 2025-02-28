# TODO: Documentation <28-11-23> 
@doc raw"""
An `N`-dimensional transfer function
"""
abstract type TransferFunction{N} end
const TF{N} = TransferFunction{N}

@doc """
If the pupil function of the system is symmetric, the OTF as well as the PSF are radially symmteric which can be used to
optimize the calculations.
"""
@traitdef RadiallySymmetric{TF<:TransferFunction}

abstract type OpticalTransferFunction{N} <: TransferFunction{N} end

@doc raw"""
# Implementation

To create a new OTF model `A <: ModelOTF{2}`, you must define the transfer coefficient at a given [frequency](@ref
TransferFunctions.Frequency):
+ `attenuation(model::A, kx::Frequency, ky::Frequency)` for a non-symmetric OTF, or
+ `attenuation(model::A, kᵣ::Frequency)` for a radially symmetric OTF and use `@traitimpl RadiallySymmetric{A}` mark the
    trait implementation.

Optionally implement the methods:
+ `preferred_type(::Type{<:ModelOTF})::T` which returns the natural type `T` that is returned from calling the model. It
    should be either `T<:Complex` or `T<:Real`. Defaults to `Float64`.
For a symmetric OTF, you should consider defining:
+  `cutoff(model::A, [a=0])` returning the largest frequency `f` such that `attenuation(model, f) >= a` if `a > 0` and
    `attenuation(model, f) > 0` if `a = 0`.
"""
abstract type ModelOTF{N} <: OpticalTransferFunction{N} end
@traitfn cutoff(tf::TF) where {TF <: ModelOTF; RadiallySymmetric{TF}} = cutoff(tf, 0.0)

abstract type PointSpreadFunction{N} <: TransferFunction{N} end
abstract type ModelPSF{N} <: PointSpreadFunction{N} end

@doc raw"""
`Union` type for a transfer function that is based on a physical model of an optical system. Contrary to
a [`MeasuredTransferFunction`](@ref), a `ModelTransferFunction` must be quantifiable at any point (in either the spatial
or the frequency domain).
"""
const ModelTransferFunction{N} = Union{ModelPSF{N},ModelOTF{N}}

const PixelSize{N} = NTuple{N,Length}
fillsize(Δxy::Length, n::Int)::PixelSize = fill(Δxy, n) |> Tuple

const Coordinate{N} = NTuple{N,Real}

# TODO: Implement a Base.show() in one method for types ModelPSF, ModelOTF using `typeof` and `nameof` <19-12-24> 
