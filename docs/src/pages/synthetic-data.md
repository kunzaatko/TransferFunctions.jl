```@meta
CurrentModule = TransferFunctions.SyntheticData
CollapsedDocStrings = true
```

```@setup synthetic-data
using TransferFunctions
using MakieMaestro.Recipes
```

# Synthetic Data

For generating synthetic data there is a module `SyntheticData` that is exported from `TransferFunctions`.

```@docs
SyntheticData
```

Synthetic data models are a subtype of `SyntheticModel`
```@docs
SyntheticModel
```

A ground truth sampled image from the model can then be generated using `groundtruth`
```@docs
groundtruth
```

## Sub-diffraction Beads

```@docs
Beads
bead
beads
```

Generating ``N=100`` beads with a diameter of ``100\;\mathrm{nm}`` in an image of the size ``256 \times 256`` pixels

```@example synthetic-data
beads = SyntheticData.beads(100, 100u"nm", 40u"nm", (256,256)) 
GT = SyntheticData.groundtruth(beads)
nothing # hide
```

```@makie synthetic-data; basename="synthetic_data_beads"
f,ax,_ = Recipes.image(Makie.convert_arguments(Makie.ImageLike(), GT)...)
scalebar!(ax,GT)
f
```
