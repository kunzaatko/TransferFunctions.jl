```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

```@setup sampled-arrays
using TransferFunctions, TestImages
using MakieMaestro.Recipes
A = testimage("mandril_color")
```

# Sampled Arrays

A `SampledArray` is a light wrapper on a `parent::AbstractArray` that holds the sampling rate in every dimension. 

!!! warning "Sampling"
    Only equidistant sampling is supported.

```@docs
SampledArray
SampledMatrix
SampledVector
```

Most relevant for microscopy images is the `SpatialArray` which is a `SampledArray` with [`Length`](@extref Unitful
`Length`) being the dimension of the sampling.

```@docs
SpatialArray
SpatialMatrix
SpatialVector
```

You can construct a `SampledArray` by providing the `parent` and the `sampling`
```@example sampled-arrays
A_sampled = SpatialArray(A, 61u"nm")
A_nonequally_sampled = SpatialArray(A, (61u"nm", 54u"nm"))
nothing # hide
```

```@makie sampled-arrays; basename="sampled_arrays", formats=:png
f,ax,_ = Recipes.image(A_sampled'; axis=(;yreversed=true))
scalebar!(ax, A_sampled; textcolor=:white, linecolor=:yellow)
f
```


## Methods on `SampledArray`s

To get the sampling vertices, you can use [`TransferFunctions.sample_vertices`](@ref)

```@docs
TransferFunctions.sample_vertices(::SampledArray)
```
