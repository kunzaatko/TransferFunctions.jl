```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

```@setup born-wolf
using TransferFunctions
```

# Born & Wolf Point Spread Function Model

The Born & Wolf model comes from the classical diffraction theory described in [Principles of Optics](@cite born2019)
and is thought of as the baseline diffraction-limited model used both in microscopy and astronomy. 

```@docs
BornWolf
```

The Born & Wolf model is derived for a finite aperture ideal lens of a perfect system. It generalizes the mathematical
formulation of the _Airy pattern_ by allowing defocus and a different aperture shape. It assumes that the only
aberration of the system is due to *defocus*. It assumes monochromatic light of a single wavelength, therefore chromatic
effects are ignored. Modern microscope objectives are designed to provide optimal imaging conditions for sources located
directly on the coverslip, in which case the Born & Wolf model is applicable (if the coverslip and immersion is used as
designed). The model disregards spherical and higher order aberrations that are due to the source of illumination being
shifted from the coverslip boundary. It is meant to be used for low to moderate ``NA`` systems (ideally ``NA <: 0.7``).

```@example born-wolf
λ = 488u"nm"    # wavelength
NA = 1.4        # numerical aperture
n = 1.5         # refractive index of medium
nothing # hide
```

```@example born-wolf
bwpsf = BornWolf(λ, NA, n)
```
