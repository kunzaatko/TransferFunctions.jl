```@meta
CurrentModule = TransferFunctions
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
psf
fit
FWHM
HWHM
```

# Non-linear Transfer Functions
```@docs
ImpulseResponseMapping
```

# Internal Utility Functions and Types
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
TransferFunctions.posaxes

TransferFunctions.fillsize

TransferFunctions.inner_axes
```

## Type Aliases
```@docs
TransferFunctions.Size
TransferFunctions.Coordinate
TransferFunctions.PixelSize
TransferFunctions.OneEdge
TransferFunctions.Edges
```

# Unit Types
```@docs
Frequency
FrequencyUnits
FrequencyFreeUnits
```
