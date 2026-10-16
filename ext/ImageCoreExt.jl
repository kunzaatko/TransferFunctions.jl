module ImageCoreExt
using ImageCore, FFTViews, FFTW
using TransferFunctions
using TransferFunctions.FFT
using ImageCore: Color1

# NOTE: FFT backends can only plan a transform for arrays whose element type is
# `Float32` or `Float64`. Fixed-point pixel types (e.g. `N0f8`) are storage formats
# that only represent values in `[0, 1]` and therefore cannot serve as the arithmetic
# type of an FFT, whose twiddle factors (e.g. `cispi(2/N)`) leave that range.
#
# `FFTW.plan_rfft` only matches `Float32`/`Float64` arrays, so a call with any other
# real eltype falls through AbstractFFTs dispatch to whichever other backend happens
# to be loaded. E.g. FFTA ≤ 0.3.1 (loaded transitively via KernelDensity.jl) defines a
# broad `plan_rfft` method for any `T <: Real` and then builds its plan and workspaces
# as `Complex{T}` for the raw input eltype, crashing with
# `ArgumentError: N0f8 is an 8-bit type ... cannot represent 2`; its 2-D planner also
# `sort!`s the region in place, so `Tuple` dims fail with
# `MethodError: no method matching sort!(::Tuple{Int64, Int64})`.
#
# Promoting with `float` (N0f8 -> Float32, Int -> Float64) before planning avoids both
# problems: the eltype becomes FFT-compatible and FFTW's own planner methods match.
# `Float16` and other non-standard float types are left for the backend to handle.
@inline _fftfloat(A::AbstractArray{<:AbstractFloat}) = A
@inline _fftfloat(A::AbstractArray) = float.(A)

# TODO: Use `channel_view` and call the underlining `corrfft` method instead of creating a new one <03-09-25> 
function TransferFunctions.corrfft(A::AbstractArray{CT}, K) where {CT<:Colorant}
    Av, dims = channelview_dims(A)
    Kc = kreshape(CT, K)
    # NOTE: FFTA ≤ 0.3.1 `sort!`s its 2-D region in place; a `Tuple` is immutable.
    dims = collect(Int, dims)
    B = rfft(_fftfloat(Av), dims)
    B .*= conj!(rfft(_fftfloat(Kc), dims))
    Avf = irfft(B, length(axes(Av, dims[1])), dims)
    colorview(base_colorant_type(CT){eltype(Avf)}, Avf)
end

function FFT.fft(A::AbstractArray{CT}) where {T<:Real, CT<:Color1{T}} 
    rfft_out = FFTW.rfft(_fftfloat(channelview(A)))
    if eltype(rfft_out) <: Complex 
        FFT.RFFTOut(rfft_out, length(axes(A,1)))
    else # FIX: Does this ever occur even? Is is possible that the `rfft` outputs a real array instead of a complex one? <03-09-25> 
        FFT.RFFTOut(colorview(base_colorant_type(CT){eltype(rfft_out)}, rfft_out), length(axes(A,1)))
    end
end

# FIX: This is currently separated only because `FFT` does not support `dims` yet. <03-09-25> 
function TransferFunctions.corrfft(A::AbstractArray{CT}, K) where {CT<:ImageCore.Color1}
    Av = channelview(A)
    Kc = kreshape(CT, K)
    B = FFT.fft(_fftfloat(Av))
    B .*= conj!(FFT.fft(_fftfloat(Kc)))
    Avf = FFT.ifft(B)
    colorview(base_colorant_type(CT){eltype(Avf)}, Avf)
end

channelview_dims(A::AbstractArray{C,N}) where {C<:Colorant,N} = channelview(A), ntuple(d -> d + 1, Val(N))
channelview_dims(A::AbstractArray{C,N}) where {C<:ImageCore.Color1,N} = channelview(A), ntuple(identity, Val(N))

function kreshape(::Type{C}, K::FFTView) where {C<:Colorant}
    Kp = parent(K)
    FFTView(reshape(K, 1, size(Kp)...))
end
kreshape(::Type{C}, K::FFTView) where {C<:ImageCore.Color1} = K
end
