```@meta
CurrentModule = TransferFunctions
```
# OTF Models

```@docs
ModelOTF
```

## Circular Pupil Optical Function


The OTF is derived from the diffraction caused by the exit pupil of the system and disregards the effect of the entrance
pupil... thus assumes no reshaping of the wavefronts in the optical system. The exit pupil, being located in the optical
system just before the light reaches the image plane, has a greater effect on the optical system OTF.

The OTF can be written in as a function of the cutoff frequency ``ρ₀`` [goodman2005a](@cite)

```math
    ℋ(ρ) = 
    \begin{cases}
    (2/π) \left\{
        \arccos(ρ/2ρ₀) - (ρ/2ρ₀)\sqrt{1 - (ρ/2ρ₀)²}
    \right\} & \text{ for } ρ ≤ 2ρ₀ \\
        0 & \text{ otherwise}.
    \end{cases}

```

The cutoff frequency can be written in terms of the wavelength ``λ``, distance between entrance pupil and the image plane ``f₂`` and the circular pupil radius ``w`` as

```math
    ρ_0 = w/(λ f₂).
```

```@docs
CircularPupilOTF
```

# PSF Models
```@docs
BornWolf
GibsonLanni
```
