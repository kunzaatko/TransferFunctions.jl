module ImageCoreExt
using ImageCore, FFTViews, FFTW
using TransferFunctions
using TransferFunctions.FFT
using ImageCore: Color1

# TODO: Use `channel_view` and call the underlining `corrfft` method instead of creating a new one <03-09-25> 
function TransferFunctions.corrfft(A::AbstractArray{CT}, K) where {CT<:Colorant}
    Av, dims = channelview_dims(A)
    Kc = kreshape(CT, K)
    B = rfft(Av, dims)
    B .*= conj!(rfft(Kc, dims))
    Avf = irfft(B, length(axes(Av, dims[1])), dims)
    colorview(base_colorant_type(CT){eltype(Avf)}, Avf)
end

function FFT.fft(A::AbstractArray{CT}) where {T<:Real, CT<:Color1{T}} 
    rfft_out = FFTW.rfft(channelview(A))
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
    B = FFT.fft(Av)
    B .*= conj!(FFT.fft(Kc))
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
