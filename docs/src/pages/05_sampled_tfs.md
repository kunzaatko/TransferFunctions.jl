```@meta
CurrentModule = TransferFunctions
```

# Sampled PSFs
```@docs
apply(::SampledOTF{N, OTF} where OTF<:TransferFunctions.OpticalTransferFunction{N}, ::AbstractArray{T, N}) where {N, T}
```

# Sampled OTFs
```@docs
apply(::SampledPSF{N, PSF} where PSF<:TransferFunctions.PointSpreadFunction{N}, ::AbstractArray{T, N}) where {T, N}
```
