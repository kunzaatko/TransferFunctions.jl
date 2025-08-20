```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

# Point Spread Functions

To implement your own subtype of a `PointSpreadFunction` you need to decide if the type you want to implement is a model
(`PSFModel`) or something else.
```@docs
TransferFunctions.PSFModel
```

If there is some symmetry in the PSF, you should implement the [`TransferFunctions.symmetry`](@ref) method which may
lead to optimization of the methods that are defined on the abstract types.

```@docs
TransferFunctions.symmetry
```

There are two types of symmetries as of now

```@docs
TransferFunctions.PointSpreadFunctionSymmetry
TransferFunctions.NoSymmetry
TransferFunctions.ZAxisRadialSymmetry
```

Next you have to implement the [`TransferFunctions.intensity`](@ref) method for the correct arguments which you have to
determine from the type of [`TransferFunctions.PointSpreadFunctionSymmetry`](@ref) that is used for the type of PSF that
you want to implement.

```@docs
TransferFunctions.intensity
```

For a subtype that has parameters that may be optimized to fit the PSF to data, you should implement the
[`TransferFunctions.params`](@ref) method 

```@docs
TransferFunctions.params
```
