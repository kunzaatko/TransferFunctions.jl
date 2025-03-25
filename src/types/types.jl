## Unitful ##
@reexport using Unitful
using Unitful: Length
@derived_dimension Frequency Unitful.𝐋^-1 true

## Abstract Types ##
"""
    TransferFunction{N}
Super-type for any `N`-dimensional transfer function

See also [`OpticalTransferFunction`](@ref), [`PointSpreadFunction`](@ref)
"""
abstract type TransferFunction{N} end
const TF{N} = TransferFunction{N} # NOTE: Internal abbreviation not meant for public use  
"""
    OpticalTransferFunction{N} <: TransferFunction{N}
Super-type for any `N`-dimensional optical transfer function

See also [`PointSpreadFunction`](@ref), [`TransferFunction`](@ref)
"""
abstract type OpticalTransferFunction{N} <: TransferFunction{N} end
"""
    PointSpreadFunction{N} <: TransferFunction{N}
Super-type for any `N`-dimensional point spread function

See also [`OpticalTransferFunction`](@ref), [`TransferFunction`](@ref)
"""
abstract type PointSpreadFunction{N} <: TransferFunction{N} end

# IDEA: Maybe this should instead be a parametric trait that holds the coefficients that the transfer functions is
# symmetric in. For instance, the radially symmetric means two arguments with the same ̃x₁² + ̃x₂². Elliptic symmetry
# could also be work if we could hold the coefficients a and b such that a ̃x₁² + b ̃x₂² same imply the same result of
# attenuation... How can one express this programmatically? It would also be useful to express axis symmetry, i.e. if x₁
# ↔ x₂ gives the same answer. <01-03-25> 
@doc """
If the pupil function of the system is symmetric, the OTF as well as the PSF are radially symmteric which can be used to
optimize the calculations.
"""
@traitdef RadiallySymmetric{TF<:TransferFunction}

# TODO: For a non-radially-symmetric OTF, the attenuation should instead be `attenuation(model::A, ks::Vararg{Frequency,
# N})` and for a larger than 2D OTF it should hold the frequencies in the larger dimensions `attenuation(model::A,
# kr::Frequency, kz::Frequency)` <01-03-25> 
# TODO: ↑ A similar argument should be made for the `cutoff` function <01-03-25> 
# TODO: Add examples of use of a ModelOTF <01-03-25> 
"""
    ModelOTF{N} <: OpticalTransferFunction{N}
An `N`-dimensional optical transfer function based on a physical model (in contrast to a [measurement](@ref
MeasuredOTF)).

See also [`CircularPupilOTF`](@ref)

# Implementation

To create a new OTF model `A <: ModelOTF{2}`, you must define the attenuation (transfer coefficient) at a given
[frequency](@ref TransferFunctions.Frequency):
+ `attenuation(model::A, kx::Frequency, ky::Frequency)` for a non-symmetric OTF, or
+ `attenuation(model::A, kᵣ::Frequency)` for a radially symmetric OTF and use `@traitimpl RadiallySymmetric{A}` to mark
    the symmetry trait implementation.

Optionally implement the methods:
+ `preferred_type(::Type{<:ModelOTF})::T` which returns the natural type `T` that is returned from calling the model. It
    should be either `T<:Complex` or `T<:Real`. Defaults to `Float64`.
For a symmetric OTF, you should consider defining:
+  `cutoff(model::A, [a=0])` returning the largest frequency `f` such that `attenuation(model, f) >= a` if `a > 0` and
    `attenuation(model, f) > 0` if `a = 0`.
"""
abstract type ModelOTF{N} <: OpticalTransferFunction{N} end
@traitfn cutoff(tf::TF) where {TF <: ModelOTF; RadiallySymmetric{TF}} = cutoff(tf, 0.0)

# TODO: Documentation. Similarly to the ModelOTF <01-03-25> 
abstract type ModelPSF{N} <: PointSpreadFunction{N} end

@doc raw"""
`Union` type for a transfer function that is based on a physical model of an optical system. Contrary to
a [`MeasuredTransferFunction`](@ref), a `ModelTransferFunction` must be quantifiable at any point (in either the spatial
or the frequency domain).
"""
const ModelTransferFunction{N} = Union{ModelPSF{N},ModelOTF{N}}

const PixelSize{N} = NTuple{N,Length} # NOTE: Internal shortcut
fillsize(Δxy::Length, n::Int)::PixelSize = fill(Δxy, n) |> Tuple

const Coordinate{N} = NTuple{N,Real} # NOTE: Internal shortcut

struct Sampled{}
end

# TODO: Implement a Base.show() in one method for types ModelPSF, ModelOTF using `typeof` and `nameof` <19-12-24> 


# Measured transfer functions
include("measured-otf.jl")
include("measured-psf.jl")

@doc raw"""
A `Union` type for the measurement of a transfer function of an optical system. It can be either a PSF measurement 
[`MeasuredPSF`](@ref) or a [`MeasuredOTF`](@ref).
"""
const MeasuredTransferFunction{N} = Union{MeasuredOTF{<:Number,N},MeasuredPSF{<:Number,N}}

# Model transfer functions
include("otf.jl")
include("psf.jl")
include("pupil.jl")

# Sampled transfer function API
include("sampled-otf.jl")
include("sampled-psf.jl")

"""
    SampledTransferFunction{N}
A `Union` type of `N`-dimensional sampled transfer functions.
"""
const SampledTransferFunction{N} = Union{SampledOTF{N},SampledPSF{N},MeasuredTransferFunction{N}}

include("function-overloads.jl")
