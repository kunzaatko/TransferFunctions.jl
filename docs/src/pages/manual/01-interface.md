```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

# Transfer Functions

Any _transfer function_ is a subtype of the abstract [`TransferFunction`](@ref)

```@docs
TransferFunction
```

A transfer function can be either linear or non-linear and shift invariant or shift variant.

```@docs
LinearTransferFunction
NonLinearTransferFunction
```

Any fully prescribed transfer function can be used to perform the forward pass (i.e. optical transfer)

```@docs
TransferFunctions.transfer(::TransferFunction, ::SpatialMatrix{<:Real})
```

The reverse is much harder and we need estimation and inverse modeling methods for it

```@docs
restore(::TransferFunction, ::SpatialMatrix{<:Real})
```
