```@meta
CurrentModule = TransferFunctions
```
There are two main categories of transfer function types that are provided by this package
+ Subtypes of [`TransferFunction`](@ref) and
+ Subtypes of [`SampledTransferFunction`](@ref).
The former facilitate the models or measured transfer functions that describe the transfer of a microscope apparatus.
The latter can be generated from the former by supplying the parameters of the sampling process, most commonly the pixel
    size.
A [`SampledTransferFunction`](@ref) is most commonly what you would want for processing data from your acquisition.

# Transfer Functions
You can describe a transfer of a microscope with either the _optical transfer function_ (OTF) the _point spread
    function_ (PSF) or the _generalized pupil function_.
These descriptions are mostly interchangeable/convertible between each other.
However note, that it may the case that a model can be in a closed form expression for some of these and must be
    approximated for the others, which may be a reason for a choice between them.

You can get a transfer function of your optical setup by supplying parameters to a model transfer function that is
    developed from the underlining physics of a microscope, then you will use a subtype of
```@docs
ModelTransferFunction
```
or you can estimate the transfer function, most commonly by the means of an acquisition where the imaged sample is known
such as sub-diffraction sized microspheres of known 
```@docs
MeasuredTransferFunction
```

# Sampled Transfer Functions
