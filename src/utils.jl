using ImageFiltering: padarray
# NOTE: For `taperedges` border <02-09-24>  
@reexport using ImageFiltering: Inner, Pad, Fill
using SpecialFunctions

# NOTE: Currently it is not used anywhere else anyway, so it is quite cheap to just throw it away. <27-08-24> 
# TODO: This should be probably left to the `ImageFiltering.jl` package. It already has a function that pads the
# needed amount <26-08-24> 
# TODO: Add examples with transfer functions. For example an model OTF, when `ifft`ed with padding of the same parity
# should be real...
@doc raw"""
    TransferFunctions.padtosize(a::AbstractArray{T,N}, size...; fourier=false, padvalue=0)
    TransferFunctions.padtosize(a::AbstractArray{T,N}, size; fourier=false, padvalue=0)

Pad `a` to a `size`

Padding is done in a way which ensures that if `a` is symmetric and decreases to zero on the edges (!!) such that
`eltype(fft(a)) <: Real`, then, if possible, output is symmetric under the DFT and `elyptype(fft(a_padded)) <: Real`).
If `a` is in the Fourier domain, then `fourier` should de set to `true`.

# Arguments:
- `size`: size of the output array. If `size isa Integer` then size along all the dimensions is `size`
- `fourier`: `a` is in the Fourier domain (default: `false`).

# Examples:
```jldoctest; setup = :(using TransferFunctions: TransferFunctions as TF; using FFTW)
julia> TF.padtosize(reshape(1:4, 2,2), 3, 3)
3×3 Matrix{Int64}:
 0  0  0
 0  1  3
 0  2  4

julia> TF.padtosize(reshape(1:4, 2,2), 4, 4)
4×4 Matrix{Int64}:
 0  0  0  0
 0  1  3  0
 0  2  4  0
 0  0  0  0

julia> TF.padtosize(reshape(1:4, 2,2), 3, 3; fourier=true)
3×3 Matrix{Int64}:
 4  2  0
 3  1  0
 0  0  0

julia> TF.padtosize(reshape(1:4, 2,2), 4, 4; fourier=true)
4×4 Matrix{Int64}:
 1  0  0  3
 0  0  0  0
 0  0  0  0
 2  0  0  4

julia> ones(1,1) |> fft
1×1 Matrix{ComplexF64}:
 1.0 + 0.0im

julia> TF.padtosize(ones(1,1), 2,2) |> fft
2×2 Matrix{ComplexF64}:
  1.0+0.0im  -1.0+0.0im
 -1.0+0.0im   1.0+0.0im

julia> ones(2,2) |> fft
2×2 Matrix{ComplexF64}:
 4.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im

julia> TF.padtosize(ones(2,2), 3, 3) |> fft # can be made symmetric
3×3 Matrix{ComplexF64}:
  4.0+0.0im  -2.0+0.0im  -2.0+0.0im
 -2.0+0.0im   1.0+0.0im   1.0+0.0im
 -2.0+0.0im   1.0+0.0im   1.0+0.0im

julia> TF.padtosize(ones(2,2), 4, 4) |> fft # cannot be symmetric
4×4 Matrix{ComplexF64}:
  4.0+0.0im  -2.0-2.0im  0.0+0.0im  -2.0+2.0im
 -2.0-2.0im   0.0+2.0im  0.0+0.0im   2.0+0.0im
  0.0+0.0im   0.0+0.0im  0.0+0.0im   0.0+0.0im
 -2.0+2.0im   2.0+0.0im  0.0+0.0im   0.0-2.0im

julia> @assert all(-10^-5 .< imag(fft(TF.padtosize(ones(2,2), 5, 5))) .< 10^-5) # .== 0 (numerical errors)

```
"""
function padtosize(
 a::AbstractArray{T,N}, s::Vararg{Int,N}; fourier=false, padvalue=zero(T)
) where {N,T}
 padxyzpost = Tuple(floor.(Int, (s .- size(a)) ./ 2))
 padxyzpre = padxyzpost .+ isodd.(size(a) .- s)

 # FIX: Is `fourier` correct?! <10-12-23> 
 # PERF: Should be using fftshift and ifftshift views  
 a = fourier ? fftshift(a) : a
 padded = padarray(a, Fill(padvalue, padxyzpre, padxyzpost)).parent
 return fourier ? ifftshift(padded) : padded
end
padtosize(a, size::Integer; vargs...) = padtosize(a, fill(size, ndims(a))...; vargs...)

@doc raw"""
    TransferFunctions.roundupcenter(arr::AbstractArray)

Calculate center index of `arr` rounded-up (i.e. `fft` center)

# Examples
```jldoctest; setup = :(using TransferFunctions: TransferFunctions as TF)
julia> TF.roundupcenter(ones(15,15))
(8, 8)

julia> TF.roundupcenter(ones(16,16))
(9, 9)

julia> TF.roundupcenter(ones(15,15,15))
(8, 8, 8)
```
"""
function roundupcenter(arr::AbstractArray{N})::Coordinate where {N}
 axs = axes(arr)
 l = @. maximum(axs) .- minimum(axs)
 return @. minimum(axs) + round(Int, l ./ 2, RoundUp)
end

@doc raw"""
     TransferFunctions.exactcenter(a::AbstractArray)

Calculate the exact center index of `a`

# Examples
```jldoctest; setup = :(using TransferFunctions: TransferFunctions as TF)
julia> TF.exactcenter(ones(15,15))
(8.0, 8.0)

julia> TF.exactcenter(ones(16,16))
(8.5, 8.5)

julia> TF.exactcenter(ones(15,15,15))
(8.0, 8.0, 8.0)
```
"""
function exactcenter(arr::AbstractArray{N})::Coordinate where {N}
 axs = axes(arr)
 return @. minimum(axs) + (maximum(axs) - minimum(axs)) / 2
end

contained(arr::AbstractArray{T,N}, loc::Coordinate{N}) where {T,N} = all(loc .∈ axes(arr))

# TODO: Documentation <02-09-24> 
# TODO: test <28-08-24> 
# TODO: Should be for N dimensions <28-08-24> 
function freqs(inds::Indices{N}, Δxy::PixelSize{N}, center::Coordinate{N}) where {N}
 fxs, fys = ndgrid(fftfreq(length(inds[1]), 1 / Δxy[1]), fftfreq(length(inds[2]), 1 / Δxy[2]))
 fxs = fxs .- (center[1] - 1) / (Δxy[1] * length(inds[1]))
 fys = fys .- (center[2] - 1) / (Δxy[2] * length(inds[2]))
 return fxs, fys
end

# utils
include("apodization.jl")
