"""
    module SyntheticData
A module that holds the methods for generating synthetic data as imitations of the objects that are relevant in
microscopy.

Every model is a subtype of the [`SyntheticData.SyntheticModel`](@ref) abstract type. Every instance of this type has
a [`SyntheticData.groundtruth`](@ref) interface which returns a [`SpatialArray`](@ref) with the sampled model.
"""
module SyntheticData # methods in TestExt
using InterfaceFunctions

"""
    SyntheticModel{N}
`SyntheticModel{N}` is an abstract type that represents a model which allows generating a ground truth image of
dimensionality `N`.

It itself may be randomized but the generation of the ground truth image from the instance must be deterministic.

Subtypes include [`Beads`](@ref)
"""
abstract type SyntheticModel{N} end

"""
    groundtruth(s::SyntheticModel)
Generate the ground truth image for the `SyntheticModel` `s`.
"""
@interface function groundtruth(::SyntheticModel) end

include("synthetic-data/beads.jl")

end

export SyntheticData
