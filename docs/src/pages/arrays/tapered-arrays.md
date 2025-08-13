```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

```@setup tapered-arrays
using TransferFunctions, TestImages
using MakieMaestro.Recipes
using OffsetArrays: no_offset_view
A = testimage("mandril_color")
```

# Tapered Arrays

A `TaperedArray` is a light wrapper on a `parent::AbstractArray` that applies a windowing functions to the arrays edges
to taper them down progressively.

```@docs
TransferFunctions.TaperedArray
taperedges
```

A `TaperedArray` is useful for reducing the boundary effects in of the Fourier images of an array. It is desirable for
making calculations in the Fourier domain to determine some value or a parameter. 

!!! tip
    When you want to pad the array for the purpose of filtering a [`BorderArray`](@ref) is relevant to you.

```@example tapered-arrays
A_tapered_inner = taperedges(A, 30) # edges of size 30 are tapered
@assert size(A_tapered_inner) == size(A) # edges are not extended
A_tapered_padded = taperedges(A, 30, :replicate)
@assert size(A_tapered_padded) == size(A) .+ (60,60) # edges are extended at both sides
nothing # hide
```

```@makie tapered-arrays; basename="tapered_arrays_mosaic", formats=:png
f,_,_ = Recipes.mosaic(Recipes.image!,collect(no_offset_view(A_tapered_inner))', collect(no_offset_view(A_tapered_padded))'; axis=(;title=["Inner Taper", "Extended Taper"], yreversed=true), linkaxes=false)
f
```

You can also select an [apodization function](@ref Apodization) to use for the tapering.
```@example tapered-arrays
using TransferFunctions: Apodization as Apo
A_tapered_hann = taperedges(Apo.Hann(), A, 30) 
nothing # hide
```

```@makie tapered-arrays; basename="tapered_arrays_hann", formats=:png
Recipes.image(collect(no_offset_view(A_tapered_hann))'; axis=(;yreversed=true))
```
