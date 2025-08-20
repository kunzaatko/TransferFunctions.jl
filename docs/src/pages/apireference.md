```@meta
CurrentModule = TransferFunctions
```

# Unit Types
```@docs
Frequency
FrequencyUnits
FrequencyFreeUnits
```

# Optical Transfer Functions
```@docs
OTFArray
RadialOTF
```

```@docs
otf
insupport
cutoff
attenuation
```

# Point Spread Functions
```@docs
PSFModel
```

```@docs
psf
TransferFunctions.intensity
fit
FWHM
HWHM
```


# Non-linear Transfer Functions
```@docs
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
TransferFunctions.params
TransferFunctions.PointSpreadFunctionSymmetry
TransferFunctions.NoSymmetry
TransferFunctions.ZAxisRadialSymmetry
```
