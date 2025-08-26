```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

```@setup gaussian
using TransferFunctions
```

# Gaussian Point Spread Function Model

```@docs
IsotropicGaussian
```

```@docs
TransferFunctions.C_AiryDisc_lateral
TransferFunctions.C_AiryDisc_axial
```

# Characteristics

```@docs
TransferFunctions.encircled_energy(::IsotropicGaussian{2}, ::Length)
TransferFunctions.energy_radius(::IsotropicGaussian{2}, ::Real)
```
