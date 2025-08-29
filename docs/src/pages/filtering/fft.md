```@meta
CurrentModule = TransferFunctions.FFT
CollapsedDocStrings = true
```

# Filtering

To determine how to perform element-wise multiplication of the output types of [`fft`](@extref `AbstractFFTs.fft`) and
[`ifft`](@extref `AbstractFFTs.ifft`) we use the `FFTOut` and `RFFTOut` types which are returned based on the element
type of the input to the [`fft`](@extref `AbstractFFTs.fft`) function. If the element type is real, it is possible to
half the number of operations by using [`rfft`](@extref `AbstractFFTs.rfft`) instead of [`fft`](@extref
`AbstractFFTs.fft`) however this returns only one half of the adjoint-symmetric array. This inconsistency is bridged by
defining multiplication on `FFTOut` and `RFFTOut` types.

```@docs
FFT
FFTOut
RFFTOut
fft
ifft
```
