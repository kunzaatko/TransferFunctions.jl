```@meta
CurrentModule = TransferFunctions
```
# Sampled Transfer Functions
```@docs
SampledTransferFunction
```

## Sampled PSFs
```@docs
SampledPSF
apply(::SampledOTF{N, OTF} where OTF<:TransferFunctions.OpticalTransferFunction{N}, ::AbstractArray{T, N}) where {N, T}
```

## Sampled OTFs
```@docs
SampledOTF
apply(::SampledPSF{N, PSF} where PSF<:TransferFunctions.PointSpreadFunction{N}, ::AbstractArray{T, N}) where {T, N}
```
