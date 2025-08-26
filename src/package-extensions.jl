function scalebar end
function scalebar! end
function scalebarformat end

export scalebar, scalebar!

module SyntheticData # methods in TestExt
using InterfaceFunctions
    
"""
    SyntheticModel{N}
`SyntheticModel{N}` is an abstract type that represents a model which allows generating a ground truth image of
dimensionality `N`.

It itself may be randomized but the generation of the ground truth image itself from the instance must be deterministic.
"""
abstract type SyntheticModel{N} end

"""
    beads(N::Int, d::Length, Δxy::PixelSize, wh::Size; <kwargs>)
Generate synthetic the `Beads` `SyntheticModel` with `N` beads of the diameter `d` 

# Keyword Arguments
- `spacing::Length=2.5d`: spacing between the beads
- `maxiters::Int=30`: maximum number of attempts to generate the beads locations with the given spacing
- `xdist=Uniform(0,wh[1])`: distribution for sampling the x-coordinate
- `ydist=Uniform(0,wh[2])`: distribution for sampling the y-coordinate
- `α=0u"nm^-1`: evanescent wave attenuation constant
"""
function beads end

"""
    groundtruth(s::SyntheticModel)
Generate the ground truth image for the `SyntheticModel` `s`.
"""
@interface function groundtruth(::SyntheticModel) end
end

export SyntheticData
