# TODO: Deconvolution techniques here too. So RL-deconv and Weiner deconvolution <26-08-24> 

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
