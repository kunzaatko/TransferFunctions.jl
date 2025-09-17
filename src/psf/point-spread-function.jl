# TODO: Make ModelPSF a trait instead of a subtype. This would allow the augmentations to be models if the parents are
# a model <28-08-25> 
# TODO: Try to define with [DensityInterface](https://juliamath.github.io/DensityInterface.jl/dev/) <21-08-25> 
using InterfaceFunctions

"""
    PointSpreadFunction{N} <: LinearShiftInvariantTransferFunction{N}
A point spread function is a description of a transfer functions specifying the intensity transfer of a single point
source in the object plane into a region of the image plane.

# Implementation
To create a new Point spread function (PSF) `A <: PointSpreadFunction`, you must define the __intensity__ at a given
[length](@extref Unitful `Length`) coordinate `intensity(psf::A, x::Length, y::Length)`.
"""
abstract type PointSpreadFunction{N} <: LinearShiftInvariantTransferFunction{N} end
Broadcast.broadcastable(tf::PointSpreadFunction) = Ref(tf)

"""
    PointSpreadFunctionSymmetry
A marker type that indicates the symmetries of a `PointSpreadFunction` which may lead to optimized methods.

Subtypes thus far include [`ZAxisRadialSymmetry`](@ref) and [`NoSymmetry`](@ref).
"""
abstract type PointSpreadFunctionSymmetry end
Broadcast.broadcastable(s::PointSpreadFunctionSymmetry) = Ref(s)

"""
    ZAxisRadialSymmetry <: PointSpreadFunctionSymmetry
``z`` axis radial symmetry of a `PointSpreadFunction` indicates that the PSF is symmetric in the radially around the
z-axis.
"""
struct ZAxisRadialSymmetry <: PointSpreadFunctionSymmetry end

"""
    NoSymmetry <: PointSpreadFunctionSymmetry
`NoSymmetry` marks a `PointSpreadFunction` that is not symmetric in any way. You need to compute every ``x``, ``y``,
``z`` coordinate separately. 
"""
struct NoSymmetry <: PointSpreadFunctionSymmetry end

"""
    symmtery(psf::PointSpreadFunction)
Marks the symmetries of a `PointSpreadFunction`.
"""
@interface symmetry(::PointSpreadFunction) = NoSymmetry()

"""
    intensity(psf::PointSpreadFunction, x::Length, y::Length)
The intensity of a `PSF` at a given location.
"""
@interface response(psf::PointSpreadFunction, args...) = response(symmetry(psf), psf, args...)
@inline response(::ZAxisRadialSymmetry, psf::PointSpreadFunction{2}, x::Length, y::Length) = intensity(psf, hypot(x, y))
@inline response(::ZAxisRadialSymmetry, psf::PointSpreadFunction{3}, x::Length, y::Length, z::Length) = intensity(psf, hypot(x, y), z)
@inline response(::NoSymmetry, psf::PointSpreadFunction{2}, x::Length, y::Length) = intensity(psf, x, y)
@inline response(::NoSymmetry, psf::PointSpreadFunction{3}, x::Length, y::Length, z::Length) = intensity(psf, x, y, z)

"""
    PSFModel <: PointSpreadFunction 
`PointSpreadFunction` subtype for PSF models. A `PSF` model is a `PointSpreadFunction` that has parameters accessible
through [`params`](@ref) and can be fitted to data.
"""
abstract type PSFModel{N} <: PointSpreadFunction{N} end

include("methods.jl")

# Augmentations
include("./types/rotated-psf.jl")
include("./types/scaled-psf.jl")

# Models
include("./types/airy-disc.jl")
include("./types/gaussian.jl")

# TODO: Fix these models <20-08-25> 
include("./types/gibson-lanni.jl")
include("./types/born-wolf.jl")

include("./types/psf-array.jl")

include("./estimation.jl")

export response, psf, otf
