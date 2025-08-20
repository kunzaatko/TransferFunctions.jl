```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

# Filtering Arrays and Circulant Tensors

- [`TransferFunctions.kern_padding`](@ref) is used to determine the padding necessary to perform the full domain
  filtering of an array.

```@docs
TransferFunctions.kern_padding
```

- [`TransferFunctions.ind2sub`](@ref) is used in the [`FilteringMatrix`](@ref) to determine the "`kernel`" and
  "`parent`" indices from the index to the `FilteringMatrix`.

```@docs
TransferFunctions.ind2sub
```

# Border Arrays

- The error [`TransferFunctions.InvalidBorderExtent`](@ref) is used when the border type does not support the extent
  requested in the [`BorderArray`](@ref) constructor. If you implement your border type, you must implement
  [`TransferFunctions.validextension`](@ref)

```@docs
TransferFunctions.validextension
```

- If the border type that you want to implement is determined by an index map on the parent array, you should implement
  a subtype of [`TransferFunctions.IndexMapBorder`](@ref) or [`TransferFunctions.IndependentIndexMapBorder`](@ref)

```@docs
TransferFunctions.IndependentIndexMapBorder
TransferFunctions.IndexMapBorder
```

To implement the `IndexBorderMap`, you need to implement the method

```@docs
Base.getindex(::TransferFunctions.IndexMapBorder{T}, ::AbstractArray{T,N}, ::Vararg{Int64,N}) where {T,N}
```
