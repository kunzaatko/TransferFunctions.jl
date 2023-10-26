module TransferFunctions
using SimpleTraits
using Unitful
using Unitful: Length
using SpecialFunctions
using OffsetArrays: centered
using FillArrays
using FFTW
using Reexport
using LazyGrids

# IDEA: Add defocus and other aberration modifiers. Look into https://github.com/RainerHeintzmann/PointSpreadFunctions.jl
# which implements these "simulations" <24-10-23> 

abstract type TransferFunction end
output_type(::TransferFunction) = Complex{Float64}

# TODO: Update docs <24-10-23> 
@doc raw"""
An abstract type for any transfer function that is based on a physical model of an optical system. Contrary to a 
[`MeasuredTransferFunction`](@ref), a `ModelTransferFunction` must be quantifiable at any point (in either spatial or 
frequency domain). The model can be of an OTF (`ClosedFormOTFModel <: ModelTransferFunction`) or a PSF 
(`ClosedFormPSFModel <: ModelTransferFunction`) or a PupilFunction.
"""
abstract type ModelTransferFunction <: TransferFunction end

# NOTE: Allows broadcasting `func.(tf::ModelTransferFunction, a:b)`where `func` can be any of `psf`,`otf`,`mtf`, etc.
Broadcast.broadcastable(tf::ModelTransferFunction) = Ref(tf)

# TODO: Add docs for implementation <02-10-23> 
abstract type ModelPSF <: ModelTransferFunction end
abstract type ModelOTF <: ModelTransferFunction end

@doc raw"""
An abstract type for the measurement of the transfer function of an optical system. It can be either a PSF measurement 
[`MeasuredPSF`](@ref) or a [`MeasuredOTF`](@ref).
"""
abstract type MeasuredTransferFunction <: TransferFunction end

@doc """
If the pupil function of the system is symmetric, the OTF as well as the PSF are radially symmteric which can be used to optimize the calculations
"""
@traitdef RadiallySymmetric{TF<:TransferFunction}

@derived_dimension Frequency Unitful.𝐋^-1

include("otf.jl")
include("psf.jl")
include("pupil.jl")

include("SIM_transfer_function_utils.jl")

include("transfer_functions/gibson_lanni.jl")
include("transfer_functions/born_wolf.jl")
include("transfer_functions/spherical_aperture_otf.jl")
include("transfer_functions/measured_psf.jl")
include("transfer_functions/measured_otf.jl")

export psf, otf, mtf, ptf, apsf, ipsf, pupil
export cutoff_frequency, resolution_limit
export BornWolf, IdealOTFwithCurvature
export MeasuredPSF
export shift
@reexport using Unitful

end
