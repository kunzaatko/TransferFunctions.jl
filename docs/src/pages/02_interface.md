```@meta
CurrentModule = TransferFunctions
```
There are two main categories of transfer function types that are provided by this package
+ Subtypes of [`TransferFunction`](@ref) which define a general transfer function model/measurement of a microscope
    apparatus and 
+ Subtypes of [`SampledTransferFunction`](@ref) that may are created by providing a sampling scheme (most commonly
    a pixel size) to the [`TransferFunction`](@ref). These are used when you are working with specific acquisition with
    known sampling rates/schemes.

# Transfer Functions
```@docs
TransferFunction
```
You can describe a transfer of a microscope with either the _optical transfer function_ (OTF) the _point spread
    function_ (PSF) or the _generalized pupil function_.
These descriptions are mostly interchangeable/convertible between each other.
!!! note
    It may the case that a model can be in a closed form expression for some of these and must be approximated for the
    others, which may be a reason for the choice between them for a specific model.

```@docs
OpticalTransferFunction
PointSpreadFunction
```

You can get a transfer function of your optical setup by supplying parameters of the apparatus to a model transfer
    function that is developed from the underlining physics of a microscope, then you will use a subtype of
```@docs
ModelTransferFunction
```
Otherwise you can estimate the transfer function, most commonly by the means of an acquisition where the imaged sample
is known such as sub-diffraction sized microspheres of known sizes.
This approach will lead to a [`MeasuredTransferFunction`](@ref)
```@docs
MeasuredTransferFunction
```

# Sampled Transfer Functions
<!-- TODO:  <18-11-24> -->
