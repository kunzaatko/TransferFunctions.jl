module ImageCoreExt
using ImageCore, FFTViews, FFTW
using TransferFunctions

function TransferFunctions.corrfft(A::AbstractArray{CT}, krn) where {CT<:Colorant}
    Av, dims = channelview_dims(A)
    kernrs = kreshape(CT, krn)
    B = rfft(Av, dims)
    B .*= conj!(rfft(kernrs, dims))
    Avf = irfft(B, length(axes(Av, dims[1])), dims)
    colorview(base_colorant_type(CT){eltype(Avf)}, Avf)
end
channelview_dims(A::AbstractArray{C,N}) where {C<:Colorant,N} = channelview(A), ntuple(d -> d + 1, Val(N))
channelview_dims(A::AbstractArray{C,N}) where {C<:ImageCore.Color1,N} = channelview(A), ntuple(identity, Val(N))

function kreshape(::Type{C}, krn::FFTView) where {C<:Colorant}
    kern = parent(krn)
    kernrs = FFTView(reshape(kern, 1, size(kern)...))
end
kreshape(::Type{C}, krn::FFTView) where {C<:ImageCore.Color1} = krn
end
