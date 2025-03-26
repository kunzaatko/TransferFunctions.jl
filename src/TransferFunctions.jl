# TODO: DimensionalData should be added as an extension and a "Sync" for sampling. By default, we should use only
# internal representation... <11-03-25> 
# TODO: Deconvolution techniques here too. So RL-deconv and Weiner deconvolution <26-08-24> 
# TODO: Consider the PointSpreadFunction and the OpticalTransferFunction interface with `Interfaces.jl` for easier
# testing and documentation <04-03-25> 

# FIX: Adapt to dimensionality <28-11-23> 

# FIX: Generic methods without argument types should be made concrete, because if not, the error when supplying a wrong
# type of argument is not a MethodError for that function but another or worse, it can be a different error type
# entirely <12-12-23> 

# TODO: Add exported types programmatically. How do they do this in the generation of default documentation? `names(::Module)`? Does it work? <11-03-25> 
"""
Package for models, estimation, sampling and deconvolution with microscopy transfer functions.
"""
module TransferFunctions

using SimpleTraits, SpecialFunctions, FillArrays, FFTW, LazyGrids, Interpolations, Roots, IntervalSets, ColorTypes,
    ImageFiltering, OffsetArrays
using OffsetArrays: centered, center, Origin

using Reexport

using Base: Indices

# RESEARCH: In the `ImageFiltering.jl` package in `src/utils.jl` there is a function `freqkernel` which is similar to
# `psf2otf`. This should be investigated <12-08-24> 

# IDEA: Add defocus and other aberration modifiers. Look into
# https://github.com/RainerHeintzmann/PointSpreadFunctions.jl which implements these "simulations". Simulations via
# semi-groups? <24-10-23> 

### source files

# Types for representing a transfer function and abstract functions on these types
include("types/types.jl")

# Utility functions for working with transfer functions and internal functions
include("utils/utils.jl")

# transfer function models
include("models/spherical-aperture-otf.jl")
include("models/gibson-lanni.jl")
include("models/born-wolf.jl")

# Model of a Beads 
include("estimation/estimation.jl")

export psf, otf, mtf, ptf, apsf, ipsf, pupil, attenuation, support
export cutoff, resolution_limit
export BornWolf, CircularPupilOTF

export MeasuredPSF, MeasuredOTF
export SampledPSF, SampledOTF
export apply

end
