```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

# Linear Transfer Functions

Linear transfer function is an operator ``T`` that adheres to the superposition principle (i.e. linearity)

```math
T\{ax_1+bx_2\}=aT\{x_1\}+bT\{x_2\}
```

For any two object plane distributions ``x_1`` and ``x_2`` the transfer of the sum of these distributions can be written
as the sum of the transfers of the distributions individually.

Linear transfer functions simplify the analysis because we can use the well established mathematical tools from linear
algebra and linear operator theory on the systems. All though this is true in theory, the whole systems are often too
huge to be used directly in the form of matrices.

You can get a transfer function of your optical setup by supplying parameters of the apparatus to a model transfer
function that is developed from the underlining physics of a microscope, then you will use a subtype of
Otherwise you can estimate the transfer function, most commonly by the means of an acquisition where the imaged sample
is known such as sub-diffraction sized microspheres of known sizes.

## Shift Invariant Linear Transfer Functions

```@docs
LinearShiftInvariantTransferFunction
OpticalTransferFunction
```

```@docs
conv(::SpatialMatrix{<:Real}, ::LinearShiftInvariantTransferFunction)
deconv(::SpatialMatrix{<:Real}, ::LinearShiftInvariantTransferFunction)
```

## Point Spread Functions

A point spread functions is a linear shift invariant transfer function that is defines the response of the system to
a single point light source in the focal plane (or object space if the imaging in multiple 3D).

```@docs
TransferFunctions.PointSpreadFunction
```

The `response` method gives the density of the point spread function at a given location relative to its center.

```@docs
response(::PointSpreadFunction, args...) 
```

```@docs
conv(::SpatialMatrix{<:Real}, ::PointSpreadFunction{2}, ::Any)
```

# Estimation

Estimation methods for transfer functions are implemented in the `Estimation` module of the `TransferFunctions` package.

```@docs
TransferFunctions.Estimation.estimate
```

## Non-Blind Methods

Non-blind are the methods for estimation that require a pair or a set of pairs of ground truth images along with their
corresponding acquired raw images (ones that are noisy ground truth images blurred by the transfer function).

### Ground Truth

The ground truth must often also be estimated using a model of the acquired scene. One option for the pair could be an
acquisition of sub-diffraction microspheres or polymer fluorescent beads along with a model consisting of the beads
known dimensions and estimated positions in the scene.

```@docs
Estimation.LeastSquares
```

