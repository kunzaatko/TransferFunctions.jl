# FIX: Adapt to dimensionality <28-11-23> 
module TransferFunctions

using Reexport
@reexport using Unitful
using Unitful: Length
@derived_dimension Frequency Unitful.𝐋^-1

using SimpleTraits
using SpecialFunctions
using OffsetArrays: centered
using FillArrays
using FFTW
using LazyGrids

# IDEA: Add defocus and other aberration modifiers. Look into https://github.com/RainerHeintzmann/PointSpreadFunctions.jl
# which implements these "simulations" <24-10-23> 

### source files

# type system
include("common.jl")

# Functions

# POLICY: Any function that requires the image dimensions `Δxy` should have is as its last argument!
include("otf.jl")
include("psf.jl")
include("pupil.jl")

# Helper API
include("transfer_function_api.jl")

# utils
include("SIM_transfer_function_utils.jl")
include("utils.jl")

# models
include("transfer_functions/gibson_lanni.jl")
include("transfer_functions/born_wolf.jl")
include("transfer_functions/spherical_aperture_otf.jl")
include("transfer_functions/measured_psf.jl")
include("transfer_functions/measured_otf.jl")

export psf, otf, mtf, ptf, apsf, ipsf, pupil
export cutoff_frequency, resolution_limit
export BornWolf, IdealOTFwithCurvature
export MeasuredPSF, MeasuredOTF
export shift

end
