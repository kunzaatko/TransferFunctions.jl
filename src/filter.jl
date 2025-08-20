using FFTW, FFTViews
using Base: @propagate_inbounds

# TODO: Allow reusable `fft_plan` <11-08-25> 

# NOTE: Factoring may not be worth it since it results in a 0.3% loss in accuracy and slowdown in tests that were made
# with random kernels and arrays <11-08-25> 

"""
    conv([T], A, K, [border])
Discrete convolution of `A` with `K`.

`A` is extended by the [`border`](@ref "Border Types"). If `T` is given, the output will have the eltype `T`. If the
border is not specified, no border is added, which for the `FFT` behaves the same as if `:circular` border was used.

See also [`conv!`](@ref), [`corr`](@ref), [`corr!`](@ref)
"""
conv(img::AbstractArray, kernel::AbstractArray, args...) = corr(img, reflect(kernel), args...)
conv(::Type{T}, img::AbstractArray, kernel::AbstractArray, args...) where {T} = corr(T, img, reflect(kernel), args...)

"""
    conv!(out, A, K, [border])
Mutating version of [`conv`](@ref).

See [`conv`](@ref) for details.
"""
conv!(out::AbstractArray, img::AbstractArray, kernel::AbstractArray, args...) = corr!(out, img, reflect(kernel), args...)

"""
    corr([T], A, K, [border])
Discrete correlation of `A` with `K`.

`A` is extended by the [`border`](@ref "Border Types"). If `T` is given, the output will have the eltype `T`. If the
border is not specified, no border is added, which for the `FFT` behaves the same as if `:circular` border was used.

See also [`corr!`](@ref), [`conv`](@ref), [`conv!`](@ref)
"""
@inline function corr(img::AbstractArray, kernel, args...)
    corr(corr_outputtype(img, kernel), img, kernel, args...)
end
@inline function corr(::Type{T}, img::AbstractArray, kernel::AbstractArray, args...) where {T}
    corr!(similar(img, T), img, kernel, args...)
end

"""
    corr!(out, A, K, [border])
Mutating version of [`corr`](@ref).

See [`corr`](@ref) for details.
"""
@inline function corr!(out::AbstractArray, img::AbstractArray, kernel::AbstractArray, border)
    corr!(out, border_array(img, border, kern_padding(kernel)), kernel)
end
function corr!(out::AbstractArray{S,N}, A::AbstractArray{T,N}, K::AbstractArray) where {S,T,N}
    krn = FFTView(zeros(eltype(K), map(length, axes(A))))
    for I in CartesianIndices(axes(K))
        krn[I] = K[I]
    end
    Af = corrfft(A, krn)
    if map(first, axes(out)) == map(first, axes(Af))
        R = CartesianIndices(axes(out))
        copyto!(out, R, Af, R)
    else
        # Exploit the periodic boundary conditions of FFTView
        dest = FFTView(out)
        src = view(FFTView(Af), axes(dest)...)
        copyto!(dest, src)
    end
    out
end

function corrfft(A::AbstractArray{ST}, K::AbstractArray{KT}) where {ST<:Union{Real,Complex},KT<:Union{Real,Complex}}
    B = _fft(A)
    C = conj!(_fft(K))
    ifft(B .* C)
end

"""
    RFFTOut{T,N,AA<:AbstractArray{T,N}}  <: AbstractArray{T,N}
A trivial wrapper that allows to specify for [`ifft`](@extref `AbstractFFTs.ifft`) and `mul!` while exploiting the
conjugate symmetry of real arrays under the Fourier transform for real arrays.

See also [`FFTOut`](@ref)
"""
struct RFFTOut{T,N,AA<:AbstractArray{T,N}}  <: AbstractArray{T,N}
    parent::AA
    d::Int
end
Base.parent(a::RFFTOut) = (@inline; a.parent)

"""
    FFTOut{T,N,AA<:AbstractArray{T,N}}  <: AbstractArray{T,N}
A trivial wrapper that allows to specify for [`ifft`](@extref `AbstractFFTs.ifft`) and `mul!` and exploit the conjugate
symmetry of the Fourier transform in the other argument to `mul!` when it is real.

See also [`RFFTOut`](@ref)
"""
struct FFTOut{T,N,AA<:AbstractArray{T,N}}  <: AbstractArray{T,N}
    parent::AA
end
Base.parent(a::FFTOut) = (@inline; a.parent)

for type in (RFFTOut, FFTOut)
    for method in (:size, :axes)
        @eval begin
            @inline Base.$method(a::$type, args...) = $method(parent(a), args...)
        end
    end
    for method in (:getindex, :setindex!)
        @eval begin
            @propagate_inbounds Base.$method(a::$type, I...) = $method(parent(a), I...)
        end
    end
end

# NOTE: FFT followed by IFFT can be optimized using conjugate symmetry for real arrays
@inline _fft(A::AbstractArray{T}) where {T<:Real} = RFFTOut(rfft(A), length(axes(A, 1)))
@inline _fft(A::AbstractArray{T}) where {T<:Complex} = FFTOut(fft(A))
@inline AbstractFFTs.ifft(A::RFFTOut) = irfft(parent(A), A.d)
@inline AbstractFFTs.ifft(A::FFTOut) = ifft(parent(A))

# FIX: This does not work with more than 2 dimensions <12-08-25> 

# NOTE: If for one array, the optimization is used and not for the other, the two arrays do not have the same sizes
# which needs to be dealt with in the element-wise multiplication
@inline function Broadcast.broadcasted(::typeof(*), A_fft::FFTOut, B_fft::RFFTOut)
    A_fft[1, :] .*= B_fft[1, :]
    A_fft[2:(B_fft.d÷2+1), 1] .*= B_fft[2:end, 1]
    A_fft[(B_fft.d÷2+2):end, 1] .*= conj(reverse(B_fft[2:(end-iseven(B_fft.d)), 1]))
    A_fft[2:(B_fft.d÷2+1), 2:end] .*= B_fft[2:end, 2:end]
    A_fft[(B_fft.d÷2+2):end, 2:end] .*= conj(reverse(B_fft[2:(end-iseven(B_fft.d)), 2:end]))
    return A_fft
end
@inline Broadcast.broadcasted(::typeof(*), A_fft::RFFTOut, B_fft::FFTOut) = B_fft .* A_fft
@inline function Broadcast.broadcasted(::typeof(*), A_fft::RFFTOut, B_fft::RFFTOut) 
    @assert A_fft.d == B_fft.d
    RFFTOut(parent(B_fft) .* parent(A_fft), A_fft.d)
end

corr_outputtype(A::AbstractArray{S}, K) where {S} = corr_outputtype(S, K)
corr_outputtype(::Type{S}, kernel::AbstractArray{T}) where {S,T} = typeof(zero(S) * zero(T) + zero(S) * zero(T))
