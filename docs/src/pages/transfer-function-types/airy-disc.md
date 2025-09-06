```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

```@setup airy-disc
using TransferFunctions
using TransferFunctions: TransferFunctions as TF
```

# Airy Disc Model

The **Airy disc** describes the diffraction-limited point spread function (PSF) of a circular aperture, which is the
fundamental resolution limit of a conventional microscope. It arises from scalar diffraction theory under the assumption
of a perfectly aligned, aberration-free, incoherent imaging system.

It is parametrized by the wavelength (``\lambda``) of the emitted light, refractive index of the medium (``n``), and the
numerical aperture (``\mathrm{NA}``).

```@docs
AiryDisc
AiryDisc(::Length{A}, ::Real, ::Real) where {A<:Real}
```
In the lateral plane it has the mathematical form of

```math
h(r) \;=\; \left( \frac{2 J_1\!\left(\tfrac{\pi \, \mathrm{NA}}{\lambda} \, r\right)}{\tfrac{\pi \, \mathrm{NA}}{\lambda} \, r} \right)^{\!2},
```

where ``J_1`` is the Bessel function of the first kind, ``\lambda`` is the wavelength of the propagated light in the medium and ``\mathrm{NA} = n \sin \theta`` is the numerical aperture.

```@docs
TransferFunctions.intensity(::AiryDisc{<:Any, T}, ::Length) where {T}
```

In the axial direction it is defined as proportionally to the lateral plane with the ratio

```math
h(z) \;\propto\; \left( \frac{\sin\!\left(\tfrac{\pi \, \mathrm{NA}^2}{\lambda n} \, z\right)}{\tfrac{\pi \, \mathrm{NA}^2}{\lambda n} \, z} \right)^{\!2}.
```

```@docs
TransferFunctions.axialintensity(::AiryDisc{3, T}, ::Length) where {T}
TransferFunctions.intensity(::AiryDisc{3}, ::Length, ::Length)
```

With the parameters 

```@example airy-disc
λ = 488u"nm"    # wavelength
NA = 1.4        # numerical aperture
n = 1.5         # refractive index of medium
nothing # hide
```
you can define the `AiryDisc` PSF model for a 3D system as

```@example airy-disc
airypsf_3d = AiryDisc(λ, NA, n)
```
or for a 2D system as

```@example airy-disc
airypsf_2d = AiryDisc{2}(λ, NA)
```

## Characteristics

Lateral FWHM (in the focal plane) is approximately equal to

```math
\mathrm{FWHM}_\perp \;\approx\; 0.51 \,\frac{\lambda}{\mathrm{NA}}
```

which you can check by

```@example airy-disc
FWHM_lateral = TF.FWHM(airypsf_3d)[1]
round(typeof(1.0u"nm"), FWHM_lateral; digits=2) # hide
```

which gives

```@example airy-disc
FWHM_lateral * NA / λ 
round(FWHM_lateral * NA / λ; digits=3) # hide
```

The Axial FWHM (along the optical axis) is approximately given by

```math
\mathrm{FWHM}_z \;\approx\; 0.885 \,\frac{\lambda n}{\mathrm{NA}^2}
```

and can be obtained with

```@example airy-disc
FWHM_axial = TF.FWHM(airypsf_3d)[3]
round(typeof(1.0u"nm"), FWHM_axial; digits = 2) # hide
```

which gives

```@example airy-disc
FWHM_axial * NA^2 / (λ * n)
round(FWHM_axial * NA^2 / (λ * n); digits = 3) # hide
```

The energy for a given radius of the Airy disc has a closed form expression and can be computed using the method
```@docs
TransferFunctions.encircled_energy(::AiryDisc{2}, ::Length)
```

```@example airy-disc
TF.encircled_energy(airypsf_2d, 300u"nm") 
round(TF.encircled_energy(airypsf_2d, 300u"nm"); digits=3) # hide
```

For a desired contained energy the correct radius can be found by [bisection](@extref `Roots.Bisection`) which and can be computed using
the method

```@docs
TransferFunctions.energy_radius(::AiryDisc{2}, ::Real)
```

```@example airy-disc
TF.energy_radius(airypsf_2d, 0.05)
round(typeof(1.0u"nm"), TF.energy_radius(airypsf_2d, 0.05); digits = 2) # hide
```
