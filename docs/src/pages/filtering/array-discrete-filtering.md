```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

```@setup filter
using TransferFunctions,  TestImages
using MakieMaestro.Recipes
using OffsetArrays: no_offset_view
A = testimage("mandril_color")
```

# [Filtering](@id filtering-manual)

`TransferFucntions` provides some functions for discrete filtering. Both convolution and correlation have a dedicated
function. This function accepts an image `A`, kernel `K` and an optional [border](@ref "Border Types") and an output
type `T`. There are both mutating and non mutating versions of the functions.

```@docs
conv(::AbstractArray, ::AbstractArray, ::Vararg)
conv!
corr
corr!
```

<!-- TODO: Fill the docs when transfer functions are implemented in full <14-08-25> -->
