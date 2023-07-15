@doc raw"""
This module implements diffraction transfer function models that are widely used in microscopy. The base type for any transfer function is `TransferFunction` and further disambiguates to `MeasuredTransferFunction` and `ModelTransferFunction`.

!!! warning
    All of the models assume *incoherent illumination* sources (i.e. where the radiation from the source is inocoherent). This is an approximation for most optical setups, but one that is used almost universally (e.g. in fluorescence microscopy).
"""
module TransferFunctions
using SimpleTraits
using Unitful
using Unitful: Length
using SpecialFunctions
using OffsetArrays: centered
using FillArrays
using FFTW # TODO: should I use AbstractFFTs?

abstract type TransferFunction end
output_type(::TransferFunction) = Complex{Float64}

# NOTE: Allows broadcasting `func.(tf, a:b)`where `func` can be any of `psf`,`otf`,`mtf`, etc.
Broadcast.broadcastable(tf::TransferFunction) = Ref(tf)

abstract type ModelTransferFunction <: TransferFunction end
abstract type ClosedFormPSFModelTransferFunction <: ModelTransferFunction end
abstract type ClosedFormOTFModelTransferFunction <: ModelTransferFunction end
abstract type MeasuredTransferFunction <: TransferFunction end

"""
If the pupil function of the system is symmetric, the OTF as well as the PSF are radially symmteric which can be used to optimize the calculations
"""
@traitdef SymmetricPupilFunction{TF<:TransferFunction}

@derived_dimension Frequency Unitful.𝐋^-1

include("otf.jl")
include("psf.jl")
include("pupil.jl")

include("transfer_functions/gibson_lanni.jl")
include("transfer_functions/born_wolf.jl")
include("transfer_functions/spherical_aperture_otf.jl")

export psf, otf, mtf, ptf, apsf, pupil
export cutoff_frequency, resolution_limit
export BornWolf, IdealOTFwithCurvature

end
