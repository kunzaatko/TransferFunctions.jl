```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

```@setup circulant-tensors
using TransferFunctions,  TestImages, Images
using MakieMaestro.Recipes
using OffsetArrays
using OffsetArrays: no_offset_view
A = testimage("mandril_color")
```

# Circulant Tensors

A circulant tensor is a construct similar to a [filtering matrix](@ref "Filtering Matrices"). It creates "views" into
the parent array of the indices of the kernel such that if the circulant tensor is
[contracted](https://en.wikipedia.org/wiki/Tensor_contraction) at the starting dimensions with the kernel it produces the
correlation result.

```@docs
CirculantTensor
circulant
contract
```

```@example circulant-tensors
K = OffsetArray(ones(51,51) ./ 51^2, -25:25, -25:25) # origin (0,0) must be contained in the kernel
A_CT = circulant(A, K, :replicate)
A_corr_K = TransferFunctions.contract(A_CT, K)
nothing # hide
```

```@makie circulant-tensors; basename="circulant_tensors", formats=:png
f,_,_ = Recipes.mosaic(Recipes.image!, A', collect(no_offset_view(A_corr_K))'; axis=(;title=["Original", "Filtered"], yreversed=true))
f
```

We can see that the first `ndims(K)` slices of the circulant tensor are simply patches (or views) of the parent array.

```@example circulant-tensors
A_indices = CartesianIndices(A)[begin:(length(A)÷10):end][2:(end-1)]
patches = [A_CT[:,:, I] for I in A_indices]
nothing # hide
```

```@makie circulant-tensors; basename="circulant_tensors_patches", formats=:png, size=(1,1)
f,_,_ = Recipes.mosaic(Recipes.image!, [collect(no_offset_view(p))' for p in patches]...; axis=(;title=map(string∘Tuple,A_indices), yreversed=true), nrows=3)
f
```

At the patch with index `(511, 256)` you can see that the parent array is extended with `:replicate` borders.
