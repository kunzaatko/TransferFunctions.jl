```@meta
CurrentModule = TransferFunctions
```

# Unit Types
```@docs
Frequency
FrequencyUnits
FrequencyFreeUnits
```

# Transfer Functions
```@docs
TransferFunction
```

# Linear Transfer Functions
```@docs
LinearShiftInvariantTransferFunction
conv(::LinearShiftInvariantTransferFunction, ::SpatialMatrix{<:Real})
deconv(::LinearShiftInvariantTransferFunction, ::SpatialMatrix{<:Real})
```

# Optical Transfer Functions
```@docs
OpticalTransferFunction
OTFArray
RadialOTF
CircularPupilOTF
```

```@docs
otf
insupport
cutoff
attenuation
```

# Point Spread Functions
```@docs
PointSpreadFunction
MeasuredPSF
PSFArray
ModelPSF
RadialPSF
BornWolf
GibsonLanni
```

```@docs
psf
intensity
fit
FWHM
```


# Non-linear Transfer Functions
```@docs
NonLinearTransferFunction
ImpulseResponseMapping
```


# Internals
## Functions
```@docs
TransferFunctions.roundcenter
TransferFunctions.roundupcenter
TransferFunctions.rounddowncenter
TransferFunctions.exactcenter
TransferFunctions.aroundorigin

TransferFunctions.contained
TransferFunctions.interior

TransferFunctions.fftfreqs
TransferFunctions.posgrid

TransferFunctions.fillsize
```

## Type Aliases
```@docs
TransferFunctions.Size
TransferFunctions.Coordinate
TransferFunctions.PixelSize
```

## To Sort
```@docs
TransferFunctions.IndependentIndexMapBorder
TransferFunctions.IndexMapBorder
TransferFunctions.Reflect
TransferFunctions.kern_padding
TransferFunctions.Fill
TransferFunctions.ind2sub
TransferFunctions.Circular
TransferFunctions.InvalidBorderExtent
TransferFunctions.Replicate
TransferFunctions.Symmetric
```
