```@meta
CurrentModule = TransferFunctions
CollapsedDocStrings = true
```

```@setup filtering-matrices
using TransferFunctions,  TestImages, Images
using MakieMaestro.Recipes
using OffsetArrays
A = testimage("mandril_color")
```

# Filtering Matrices
Filtering matrices are a construct that is useful for making [LS estimates](@ref "Estimation") of a PSF from ground truth data. It is
a [Toeplitz matrix](https://en.wikipedia.org/wiki/Toeplitz_matrix) constructed from a parent array with the
knowledge of the axes of the filtering kernel. Discrete correlation with the kernel can then be computed by matrix
multiplication of the kernel with the filtering matrix. For vectors ``x`` and ``h``, the correlation ``y = h \ast x`` can be written as
```math
y =
      \begin{bmatrix}
      0         & 0         & \cdots & 0         & x_0       & x_1       & \cdots & x_{k-1} & x_{k} \\
      0         & 0         & \cdots & x_0       & x_1       & x_2       & \cdots & x_{k}   & x_{k+1} \\
      \vdots    & \vdots    &        & \vdots    & \vdots    & \vdots    &        & \vdots  & \vdots \\
      0         & x_0       & \cdots & x_{n-k-1} & x_{n-k}   & x_{n-k+1} & \cdots & x_{n-1} & x_n \\
      x_0       & x_1       & \cdots & x_{n-k}   & x_{n-k+1} & x_{n-k+2} & \cdots & x_n     & 0 \\
      \vdots    & \vdots    &        & \vdots    & \vdots    & \vdots    &        & \vdots  & \vdots \\
      x_{n-k-1} & x_{n-k}   & \cdots & x_{n-2}   & x_{n-1} & x_{n} & \cdots & 0       & 0 \\
      x_{n-k}   & x_{n-k+1} & \cdots & x_{n-1}   & x_n & 0 & \cdots & 0       & 0 \\
    \end{bmatrix}
    \begin{bmatrix}
        h_{-k} \\
        h_{-k+1} \\
        \vdots \\
        h_0 \\
        \vdots \\
        h_{k-1} \\
        h_{+k}
    \end{bmatrix}
```
where the matrix is the filtering matrix of ``h`` ``F_h`` zero extended. The discrete convolution can be written
analogously by reversing/reflecting ``x`` or ``h``. A filtering matrix for more than one dimensions can be used
similarly with a flattened kernel.

```@docs
TransferFunctions.FilteringMatrix
filtering_matrix
```

```@example filtering-matrices
K = OffsetArray(ones(11,11) ./ 121, -5:5, -5:5) # origin (0,0) must be contained in the kernel
F_A_inner = filtering_matrix(A, K)
A_corr_K_inner_flat = F_A_inner' * K[:]
A_corr_K_inner = reshape(A_corr_K_inner_flat, F_A_inner.interior)
@assert size(A_corr_K_inner) == (502, 502) # output of the correlation is reduced in size by the size of the kernel
nothing # hide
```

We can extend the parent array by a border (see [Border Arrays](@ref border-arrays-manual)) to obtain the whole domain of `A` in the output of the
correlation

```@example filtering-matrices
F_A = filtering_matrix(A, K, :replicate)
A_corr_K_flat = F_A' * K[:]
A_corr_K = reshape(A_corr_K_flat, size(A))
@assert size(A_corr_K) == (512, 512) # whole domain of A is retained
nothing # hide
```

```@makie filtering-matrices; basename="filtering_matrix", formats=:png
f,_,_ = Recipes.mosaic(Recipes.image!, A', A_corr_K'; axis=(;title=["Original", "Filtered"], yreversed=true))
f
```

If we resize the image, we can see the structure of the filtering matrix for a matrix.
```@example filtering-matrices
A_small = imresize(A, (20,20))
K_small = filtering_matrix(A_small, (-10:10, -10:10), :fill)
nothing # hide
```

```@makie filtering-matrices; basename="filtering_matrix_small", formats=:png
f,ax,_ = Recipes.mosaic(Recipes.image!, A_small', K_small'; axis=(;title=["Image", "Filtering Matrix"], yreversed=true),
linkaxes=false)
ax[1].xticklabelsvisible = true
ax[1].xticksvisible = true
ax[1].xticks = ([12.5], ["13"])
f
```

As can become apparent from the image, under some conditions on the border and the link between the length of `A` and `K`, the filtering matrix of a 2D matrix is a block [circulant matrix](https://en.wikipedia.org/wiki/Circulant_matrix).

If we take a column slice from the image and create a filtering matrix from it, we can observe the structure of
a filtering matrix of a vector.

```@example filtering-matrices
A_col = A_small[:, 13]
K_A_col = filtering_matrix(A_col, (-10:10,), :fill)
nothing # hide
```

```@makie filtering-matrices; formats=:png, basename="filtering_matrix_column"
f,_,_ = Recipes.mosaic(Recipes.image!, reshape(A_col, (:, 1))', K_A_col'; axis=(;title=["Column", "Filtering Matrix"], yreversed=true), linkaxes=false)
f
```
