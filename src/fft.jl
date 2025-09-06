"""
    module FFT

Optimized FFT operations exploiting conjugate symmetry for real-valued arrays.

This module provides wrapper types and optimized methods for Fast Fourier Transform operations
that take advantage of the conjugate symmetry property of real arrays under the Fourier transform.
For real-valued input arrays, the FFT output exhibits conjugate symmetry, allowing storage and
computation optimizations by only storing roughly half the frequency domain data.

## Key Features

- **Memory Optimization**: Uses `rfft` (real FFT) for real arrays to reduce memory usage by ~50%
- **Conjugate Symmetry Exploitation**: Automatically handles the symmetric properties of real FFTs
- **Seamless Integration**: Wrapper types behave like regular arrays with transparent optimizations
- **Mixed Operations**: Supports element-wise operations between optimized and standard FFT outputs

## Types

- [`RFFTOut`](@ref): Wrapper for real FFT output with conjugate symmetry optimization
- [`FFTOut`](@ref): Wrapper for complex FFT output compatible with `RFFTOut` operations

## Usage

The module automatically selects the appropriate optimization based on input type.

If you use it on a complex array, then the wrapper [`FFTOut`](@ref) is used and no optimizations are applied.

```jldoctest fft_module; output=false
A_complex = rand(ComplexF64, 64, 64)     # Array{ComplexF64} of size 64×64
FA_complex = TF.FFT.fft(A_complex)       # FFTOut of size 64×64
A_complex_ifft = TF.FFT.ifft(FA_complex) # Array{ComplexF64} of size 64×64
@assert A_complex_ifft ≈ A_complex

# output

```

On real arrays, conjugate symmetry is exploited and the wrapper [`RFFTOut`](@ref) of half the size is returned.

```jldoctest fft_module; output=false
A_real = rand(64, 64)                 # Array{Float64} of size 64×64
FA_real = TF.FFT.fft(A_real)          # RFFTOut of size 33×64
A_real_ifft = TF.FFT.ifft(FA_real)    # Array{Float64} of size 64×64
@assert A_real_ifft ≈ A_real          # size of the array is restored

# output

```

So half of the computations are shaved off.

```jldoctest fft_module; filter = r"\\s.*\\d.*"
julia> @btime TF.FFT.fft(A_real);
  13.225 μs (13 allocations: 33.60 KiB)

julia> @btime FFTW.fft(A_real);
  26.380 μs (11 allocations: 128.50 KiB)
```

These two types play together nicely though. Broadcasting works on the full array sizes if necessary (i.e. `FFTOut` and
`RFFTOut` are interacting).

```jldoctest fft_module; output=false
@assert size(FA_real .* FA_complex) == (64, 64) # FFTOut of size 64×64

# output
```
"""
module FFT
using Base: @propagate_inbounds
using Base.Broadcast: Broadcasted, ArrayStyle, AbstractArrayStyle, DefaultArrayStyle
import Base.Broadcast
using FFTW: FFTW

# PERF: parity of the arrays first dimension stored in the type for specialization on indexing
"""
    RFFTOut{T,N,AA<:AbstractArray{T,N}}  <: AbstractArray{T,N}
An abstract array wrapper for the output of [`rfft`](@extref `AbstractFFTs.rfft`) which is created when calling
[`fft`](@ref) on a real array. 

It allows to handle the array as if it were the full size but exploits conjugate symmetry of the Fourier transform by
running [`AbstractFFTs.rfft`](@extref) instead of [`AbstractFFTs.fft`](@extref) on real arrays.

See also [`FFTOut`](@ref)
"""
struct RFFTOut{T,N,AA<:AbstractArray{T,N},d,Odd}  <: AbstractArray{T,N}
    parent::AA
    RFFTOut(parent::AA, d::Int) where {T,N,AA<:AbstractArray{T,N}} = new{T,N,AA,d,isodd(d)}(parent)
end
@inline firstdim(::Type{<:RFFTOut{<:Any, <:Any, <:Any, d}}) where {d} = d
@inline firstdim(a::RFFTOut) = firstdim(typeof(a))
Base.parent(a::RFFTOut) = (@inline; a.parent)
Base.size(a::RFFTOut) = (@inline; (firstdim(a), size(a.parent)[2:end]...))

@inline oddsize(a::RFFTOut{<:Any, <:Any, <:Any, <:Any, true}) = true
@inline oddsize(a::RFFTOut{<:Any, <:Any, <:Any, <:Any, false}) = false

# NOTE: `rfft` outputs dense arrays with 1-based indexing. It does respect OffsetArrays but only in the input. I.e. the
# FFT procedure is index informed but the output is not. <29-08-25> 
@propagate_inbounds function Base.getindex(a::RFFTOut, I::Vararg{Int,N}) where {N}
    @boundscheck checkbounds(a, I...)
    conjugate = I[1] > size(parent(a), 1)
    Ip = if conjugate  
        I1 = (size(parent(a), 1) - (I[1] % size(parent(a), 1) - (oddsize(a) ? 1 : 0)))
        (I1, map(size(parent(a))[2:end], size(a)[2:end], I[2:end]) do ps, as, i
                i == 1 && return i
                ps - ((i - 1) % ps) + 1
            end...)
    else
        I
    end
    return conjugate ? conj(@inbounds(parent(a)[Ip...])) : @inbounds(parent(a)[Ip...])
end

# NOTE: `rfft` outputs dense arrays with 1-based indexing. It does respect OffsetArrays but only in the input. I.e. the
# FFT procedure is index informed but the output is not. <29-08-25> 
@propagate_inbounds function Base.setindex!(a::RFFTOut, v, I::Vararg{Int,N}) where {N}
    @boundscheck checkbounds(parent(a), I...)
    setindex!(parent(a), v, I...)
end

Base.similar(a::RFFTOut) = RFFTOut(similar(parent(a)), firstdim(a))
Base.similar(a::RFFTOut, ::Type{S}) where {S} = RFFTOut(similar(parent(a), S), firstdim(a))
Base.similar(a::Type{T}, ::Type{S}, sz) where {T <: RFFTOut,S} = RFFTOut(similar(Array{S}, sz), firstdim(T))

Broadcast.BroadcastStyle(a::Type{T}) where {T<:RFFTOut} = ArrayStyle{T}()
# TODO: Add some methods to specialize or throw an error when `N`s are different or `d`s are different <02-09-25> 
Broadcast.BroadcastStyle(::ArrayStyle{A}, ::ArrayStyle{B}) where {B<:RFFTOut, A<:RFFTOut} = ArrayStyle{promote_type(A,B)}()
# NOTE: Necessary because we want to dispatch on the other dimensionalities <03-09-25> 
Broadcast.BroadcastStyle(a::ArrayStyle{S}, ::DefaultArrayStyle{0}) where {T,N,S<:RFFTOut{T, N}} = a # NOTE: T,N need to be here for disambiguation <02-09-25> 
# TODO: Instead of constructing the type of the FFTOut, this should be delegated to a subfunction that converts the type
# of RFFTOut to the matching FFTOut <03-09-25> 
Broadcast.BroadcastStyle(::ArrayStyle{S}, b::DefaultArrayStyle{N}) where {T,N,S<:RFFTOut{T, N}} = Broadcast.BroadcastStyle(ArrayStyle{FFTOut{T, N}}(), b)
Broadcast.BroadcastStyle(::ArrayStyle{S}, b::AbstractArrayStyle{N}) where {T,N,S<:RFFTOut{T, N}} = Broadcast.BroadcastStyle(ArrayStyle{FFTOut{T, N}}(), b)

# FIX: Define conjugate stable functions. I.e. functions that where fn(conj(x)) == fn(x). On these functions, it is
# possible to use create a broadcasted with the `RFFTOut` style. For others it should promote to the `FFTOut` style.
# <02-09-25> 

@inline broadcast_args(args::Tuple) = (broadcast_args(args[1]), broadcast_args(Base.tail(args))...)
@inline broadcast_args(args::NTuple{1}) = (broadcast_args(args[1]),)
@inline broadcast_args(a) = a
@inline broadcast_args(a::RFFTOut) = parent(a)

Base.similar(bc::Broadcasted{<:ArrayStyle{T}}, ::Type{S}) where {T<:RFFTOut, S} = similar(T, S, Broadcast.combine_axes(broadcast_args(bc.args)...))

function Base.copyto!(dest::RFFTOut, bc::Broadcasted{<:ArrayStyle{<:RFFTOut}})
    copyto!(parent(dest), Broadcast.Broadcasted(bc.f, broadcast_args(bc.args)))
    return dest
end

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

Broadcast.BroadcastStyle(T::Type{<:FFTOut}) = ArrayStyle{T}()

# TODO: Consider instead using a `promote_rule` definition to promote `RFFTOut` to `FFTOut` <02-09-25> 
Broadcast.BroadcastStyle(a::ArrayStyle{A}, ::ArrayStyle{B}) where {T1,T2,M,N,A<:FFTOut{T1, M}, B<:RFFTOut{T2,N}} = ArrayStyle{FFTOut{promote_type(T1, T2), max(N,M)}}()
Broadcast.BroadcastStyle(a::ArrayStyle{A}, ::ArrayStyle{B}) where {T1,T2,M,N,A<:FFTOut{T1, M}, B<:FFTOut{T2,N}} = ArrayStyle{FFTOut{promote_type(T1, T2), max(N,M)}}()
Broadcast.BroadcastStyle(a::ArrayStyle{<:FFTOut{T, M}}, ::DefaultArrayStyle{N}) where {T, M, N} = ArrayStyle{FFTOut{T, max(N,M)}}()
Broadcast.BroadcastStyle(a::ArrayStyle{<:FFTOut{T, M}}, ::AbstractArrayStyle{N}) where {T, M, N} = ArrayStyle{FFTOut{T, max(N,M)}}()
Broadcast.BroadcastStyle(::ArrayStyle{A}, ::ArrayStyle{B}) where {B<:FFTOut, A<:FFTOut} = ArrayStyle{promote_type(A,B)}()

Base.similar(bc::Broadcasted{<:ArrayStyle{<:FFTOut}}, ::Type{S}) where {S} = FFTOut(similar(Array{S}, axes(bc)))

for method in (:size, :axes)
    @eval @inline Base.$method(a::FFTOut, args...) = $method(parent(a), args...)
end
for method in (:getindex, :setindex!)
    @eval @propagate_inbounds Base.$method(a::FFTOut, I...) = $method(parent(a), I...)
end

# NOTE: FFT followed by IFFT can be optimized using conjugate symmetry for real arrays
    
"""
    fft(A::AbstractArray)
Compute the FFT while exploiting the conjugate symmetry of real arrays and using the usual FFT for complex arrays.

Returns either an [`RFFTOut`](@ref) or a [`FFTOut`](@ref) depending on the `eltype` of the array. These types are smart
when broadcasted, i.e. if possible `RFFTOut` preserves its symmetry and if not possible, it materializes into the full
`FFTOut`.
"""
@inline fft(A::AbstractArray{T}) where {T<:Real} = RFFTOut(FFTW.rfft(A), length(axes(A, 1)))
@inline fft(A::AbstractArray{T}) where {T<:Complex} = FFTOut(FFTW.fft(A))

"""
    ifft(A::RFFTOut)
    ifft(A::FFTOut)
Computes the iFFT while using the conjugate symmetry that was exploited during the [`fft`](@ref) operation.
"""
@inline ifft(A::RFFTOut) = FFTW.irfft(parent(A), firstdim(A))
@inline ifft(A::FFTOut) = FFTW.ifft(parent(A))

# FIX: Add support for `fft` on only selected dims <03-09-25> 

end # module FFT
