```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

```@setup border-arrays
using TransferFunctions, Statistics,  TestImages
using MakieMaestro.Recipes
using OffsetArrays: no_offset_view
A = testimage("mandril_color")
```

# [Border Arrays](@id border-arrays-manual)

A `BorderArray` is a light wrapper on a `parent::AbstractArray` that expands the domain by a border with the prescribed
extent given in the constructor by `edges`. The border is added by a rule that prescribes a value to a given index
that is not contained in the domain of the `parent`.

```@docs
TransferFunctions.BorderArray
border_array
```

A `BorderArray` is useful for expanding the domain of an array during filtering (i.e. convolution or correlation) to
avoid boundary effects while preserving the original size of the array in the output. 

!!! tip 
    While this is a common practice for making calculations in the Fourier domain and then returning to the natural
    domain of the array, when you want to make some parameter estimation or want to determine a value directly from the
    Fourier images of arrays, using [`TaperedArray`](@ref)s or a combination of a [`BorderArray`](@ref) and
    a [`TaperedArray`](@ref) can be more appropriate.

# Border Types

The procedure of determining the value of the border of the `BorderArray` is determined by the `border` field which has
a subtype of `AbstractBorder`.

```@docs
TransferFunctions.AbstractBorder
```

For different applications, different border types are optimal.

```@docs
TransferFunctions.Replicate
TransferFunctions.Symmetric
TransferFunctions.Reflect
TransferFunctions.Circular
TransferFunctions.Fill
```

If you want to simply remove the boundary effects of filtering and do not care about the values of the edges of the
output (or the filtering kernel is small enough compared to the array that the edges are not important) you can fill the
border with zeros or some other value 

```@example border-arrays
# For a filtering kernel with axes (-30:30, -30:30)
A_zero_padded = border_array(A, :fill, 30) # single padding for all the dimensions and left/right edges
A_mean_padded = border_array(A, TransferFunctions.Fill(mean(A)), (30, 30)) # (left, right) padding 
A_unequally_padded = border_array(A, TransferFunctions.Fill(mean(A)), ((25, 20), (40, 35))) # In the order of axes i.e. ((top, bottom), (left, right))
nothing # hide
```

```@makie border-arrays; basename="fill", formats=:png
f,_,_ = Recipes.mosaic(Recipes.image!, collect(no_offset_view(A_zero_padded))', collect(no_offset_view(A_mean_padded))'; axis=(;title=["Zero Padded", "Mean Padded"], yreversed=true))
f
```

If the border values are important to you, it is better to use some border scheme that relates to the values of the
parent array in the locality of the border or the opposing edge of the array. The types of borders that are implemented
are [`:replicate`](@ref TransferFunctions.Replicate), [`:reflect`](@ref TransferFunctions.Reflect), [`:circular`](@ref
TransferFunctions.Circular) and [`:symmetric`](@ref TransferFunctions.Symmetric).

```@makie border-arrays; basename="border_arrays_mosaic", formats=:png, size=(1, 1)
borders = [:replicate, :reflect, :circular, :symmetric]
f,_,_ = Recipes.mosaic(borders...; axis=(;yreversed=true, title=map(string, borders)), nrows=2) do ax,border
    ba = border_array(A, border, 30)
    Recipes.image!(ax,collect(no_offset_view(ba))')
end
f
```

If the border type does not support the extent of the padding that is given in the constructor, the
[`TransferFunctions.InvalidBorderExtent`](@ref) is thrown.

```@docs
TransferFunctions.InvalidBorderExtent
```

# Functions 
Border arrays are useful for extending the domain of an array for [filtering](@ref filtering-manual).

This motivates some helper functions for sampling an extended array at given indices.

To determine the padding necessary for a given kernel in filtering, you can use `kern_padding`

```@docs
TransferFunctions.kern_padding
```

For determining the indices that are necessary in the parent to facilitate a filtering, you can use `outer_axes`
```@docs
TransferFunctions.outer_axes
```

For extending the domain of an array to some given indices with a border strategy, you can use `padtoaxes`. This may
also be useful if you do not want to have all the indices of the parent in the output of a filtering.
```@docs
padtoaxes
```
