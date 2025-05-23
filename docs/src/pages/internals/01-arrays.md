```@meta
CurrentModule = TransferFunctions
```

# Abstract Array Subtypes

## General Storage Arrays

```@docs; canonical=false
Flattened
```

A container array that is used in the implementation of other array types that holds some dimensions in an inner array
and these 

## Filtering Arrays

```@docs; canonical=false
CirculantTensor
```

```@docs; canonical=false
FilteringMatrix
```

## Data Arrays

For representing data that is sampled at some constant rate, possibly different in every dimension, there are
`SampledArray`s that store the data coupled with its sampling.
```@docs; canonical=false
SampledArray
SampledMatrix
SampledVector
```

Since sampling at spatially defined lattice points is so common in the context of microscopy, constants for
`SampledArray`s with spatial sampling are provided as aliases.
```@docs; canonical=false
SpatialArray
SpatialMatrix
SpatialVector
```
