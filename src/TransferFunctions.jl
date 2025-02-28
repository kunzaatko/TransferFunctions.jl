# TODO: Deconvolution techniques here too. So RL-deconv and Weiner deconvolution <26-08-24> 
# FIX: Adapt to dimensionality <28-11-23> 

# FIX: Generic methods without argument types should be made concrete, because if not, the error when supplying a wrong
# type of argument is not a MethodError for that function but another or worse, it can be a different error type
# entirely <12-12-23> 

module TransferFunctions

using SimpleTraits, SpecialFunctions, FillArrays, FFTW, LazyGrids, Interpolations, Roots, IntervalSets, ColorTypes, ImageFiltering, OffsetArrays
using OffsetArrays: centered, center, Origin

using Reexport
@reexport using Unitful
using Unitful: Length
@derived_dimension Frequency Unitful.𝐋^-1 true

using Base: Indices

# RESEARCH: In the `ImageFiltering.jl` package in `src/utils.jl` there is a function `freqkernel` which is similar to
# `psf2otf`. This should be investigated <12-08-24> 

# IDEA: Add defocus and other aberration modifiers. Look into
# https://github.com/RainerHeintzmann/PointSpreadFunctions.jl which implements these "simulations" <24-10-23> 

### source files

# Common type system and `Base` function overloads
include("types.jl")
include("common.jl")
include("utils.jl")

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

include("sampled-tf-base-overloads.jl")


export psf, otf, mtf, ptf, apsf, ipsf, pupil, attenuation, support
export cutoff, resolution_limit
export BornWolf, IdealOTFwithCurvature

export MeasuredPSF, MeasuredOTF
export SampledPSF, SampledOTF
export apply

end
