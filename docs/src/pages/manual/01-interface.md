```@meta
CurrentModule = TransferFunctions
```

# User Interface

Any _transfer function_ is a subtype of the abstract type [`TransferFunction`](@ref)

```@docs; canonical = false
TransferFunction
```

```@docs; canonical = false
OpticalTransferFunction
PointSpreadFunction
```

You can get a transfer function of your optical setup by supplying parameters of the apparatus to a model transfer
    function that is developed from the underlining physics of a microscope, then you will use a subtype of
Otherwise you can estimate the transfer function, most commonly by the means of an acquisition where the imaged sample
is known such as sub-diffraction sized microspheres of known sizes.
