```@meta
CurrentModule = TransferFunctions
```

# Unit Types
```@docs
Frequency
FrequencyUnits
FrequencyFreeUnits
```

# Array Types
```@docs
Flattened
flatten
CirculantTensor
circulant

FilteringMatrix
SampledArray
SampledMatrix
SampledVector
SpatialArray
SpatialMatrix
SpatialVector
```

# Transfer Functions
```@docs
TransferFunction
```

# Linear Transfer Functions
```@docs
LinearTransferFunction
conv
deconv
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

# Estimation
```@docs
Estimation.bead
Estimation.LeastSquares
```
# Apodization

```@docs
Apodization

Apodization.apodize
Apodization.taperedges

Apodization.ApodizationFunction

Apodization.Blackman
Apodization.ExactBlackman
Apodization.Connes
Apodization.Cosine
Apodization.Gaussian
Apodization.Hamming
Apodization.Welch
Apodization.BlackmanNuttall
Apodization.PowerCosine
Apodization.Triangular
Apodization.Nuttall
Apodization.SineSum
Apodization.BlackmanHarris
Apodization.FlatTop
Apodization.Hann
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
