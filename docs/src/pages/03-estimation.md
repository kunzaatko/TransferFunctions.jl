```@meta
CurrentModule = TransferFunctions.Estimation
```

# Estimation

Estimation methods for transfer functions are implemented in the `Estimation` module of the `TransferFucntions` package.

## Non-Blind Methods

Non-blind are the methods for estimation that require a pair or a set of pairs of ground truth images along with their
corresponding acquired raw images (ones that are noisy ground truth images blurred by the transfer function).

### Ground Truth

The ground truth must often also be estimated using a model of the acquired scene. One option for the pair could be an
acquisition of sub-diffraction microspheres or polymer fluorescent beads along with a model consisting of the beads
known dimensions and estimated positions in the scene.

```@docs; canonical=false
bead
```

