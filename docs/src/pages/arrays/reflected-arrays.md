```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

```@setup reflected-arrays
using TransferFunctions, TestImages
using MakieMaestro.Recipes
using OffsetArrays
using OffsetArrays: no_offset_view
A = testimage("mandril_color")
```

# Reflected Arrays

A `ReflectedArray` is a light wrapper on a `parent::AbstractArray` that reflects the array in all dimensions.

```@docs
ReflectedArray
ReflectedMatrix
ReflectedVector
reflect
```

It is useful for using discrete convolution instead of correlation in the general `filter` function.

```@example reflected-arrays
A_reflected = reflect(A)
nothing # hide
```

```@makie reflected-arrays; basename="reflected_arrays", formats=:png
f,_,_ = Recipes.mosaic(Recipes.image!, A', collect(no_offset_view(A_reflected))'; axis=(;title=["Image", "Reflected"], yreversed=true),
linkaxes=false)
f
```
