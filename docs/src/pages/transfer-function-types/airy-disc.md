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

```@docs
AiryDisc
AiryDisc(::Length{A}, ::Real, ::Real) where {A<:Real}
```
In the lateral plane it has the mathematical form of

```math
h(r) \;=\; \left( \frac{2 J_1\!\left(\tfrac{\pi \, \mathrm{NA}}{\lambda} \, r\right)}{\tfrac{\pi \, \mathrm{NA}}{\lambda} \, r} \right)^{\!2},
```

where ``J_1`` is the Bessel function of the first kind, ``\lambda`` is the wavelength of the propagated light in the medium and ``\mathrm{NA} = n \sin \theta`` is the numerical aperture.

In the axial direction it is defined as proportionally to the lateral plane with the ratio

```math
h(z) \;\propto\; \left( \frac{\sin\!\left(\tfrac{\pi \, \mathrm{NA}^2}{\lambda n} \, z\right)}{\tfrac{\pi \, \mathrm{NA}^2}{\lambda n} \, z} \right)^{\!2}.
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
```

which gives

```@example airy-disc
FWHM_lateral * NA / λ 
```

The Axial FWHM (along the optical axis) is approximately given by

```math
\mathrm{FWHM}_z \;\approx\; 0.885 \,\frac{\lambda n}{\mathrm{NA}^2}
```

and can be obtained with

```@example airy-disc
FWHM_axial = TF.FWHM(airypsf_3d)[3]
```

which gives

```@example airy-disc
FWHM_axial * NA^2 / (λ * n)
```

## Properties

* Radially symmetric intensity profile with a bright central maximum and concentric rings.
* Sets the fundamental resolution limit for widefield microscopy.
* Depends only on wavelength (``\lambda``), refractive index of the medium (``n``), and the numerical aperture (``\mathrm{NA}``).
* Widely used as a reference PSF model and as a basis for Gaussian approximations.
