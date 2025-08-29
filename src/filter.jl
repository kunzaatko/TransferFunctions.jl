using FFTViews

include("fft.jl")

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
conv(A::AbstractArray, K::AbstractArray, args...) = corr(A, reflect(K), args...)
conv(::Type{T}, A::AbstractArray, K::AbstractArray, args...) where {T} = corr(T, A, reflect(K), args...)

"""
    conv!(out, A, K, [border])
Mutating version of [`conv`](@ref).

See [`conv`](@ref) for details.
"""
conv!(out::AbstractArray, A::AbstractArray, K::AbstractArray, args...) = corr!(out, A, reflect(K), args...)

"""
    corr([T], A, K, [border])
Discrete correlation of `A` with `K`.

`A` is extended by the [`border`](@ref "Border Types"). If `T` is given, the output will have the eltype `T`. If the
border is not specified, no border is added, which for the `FFT` behaves the same as if `:circular` border was used.

See also [`corr!`](@ref), [`conv`](@ref), [`conv!`](@ref)
"""
@inline function corr(A::AbstractArray, K, args...)
    corr(corr_outputtype(A, K), A, K, args...)
end
@inline function corr(::Type{T}, A::AbstractArray, K::AbstractArray, args...) where {T}
    corr!(similar(A, T), A, K, args...)
end

"""
    corr!(out, A, K, [border])
Mutating version of [`corr`](@ref).

See [`corr`](@ref) for details.
"""
@inline function corr!(out::AbstractArray, A::AbstractArray, K::AbstractArray, border)
    corr!(out, border_array(A, border, kern_padding(K)), K)
end
function corr!(out::AbstractArray{S,N}, A::AbstractArray{T,N}, K::AbstractArray) where {S,T,N}
    Kv = FFTView(zeros(eltype(K), map(length, axes(A))))
    for I in CartesianIndices(axes(K))
        Kv[I] = K[I]
    end
    Af = corrfft(A, Kv)
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
    B = FFT.fft(A)
    C = conj!(FFT.fft(K))
    FFT.ifft(B .* C)
end

corr_outputtype(A::AbstractArray{S}, K) where {S} = corr_outputtype(S, K)
corr_outputtype(::Type{S}, kernel::AbstractArray{T}) where {S,T} = typeof(zero(S) * zero(T) + zero(S) * zero(T))
