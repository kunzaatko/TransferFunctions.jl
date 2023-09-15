@doc """
  Intensity point spread function i.e. the intensity ratio and phase shift of the sample intensity density
 """ psf

@doc """
psf(tf::TransferFunction, x::Length, y::Length)

Sample the PSF of the transfer function model/data at the point `(x,y)`

```jldoctest
julia> tf = BornWolf(488u"nm", 1.4, 1.7)
BornWolf{Float64}(488.0 nm, 1.4, 1.7)
```

```@repl
tf = BornWolf(488u"nm", 1.4, 1.7) # hide
psf(tf, 0u"nm", 5u"nm")
psf.(tf, 0u"nm", -400u"nm":100u"nm":400u"nm")
```
"""
@inline @traitfn function psf(tf::TF, x::Length, y::Length) where {TF <: ClosedFormPSFModel; SymmetricPupilFunction{TF}}
    psf(tf, hypot(x, y))
end

@traitfn function psf(tf::TF, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length}) where {TF <: ClosedFormPSFModel; !SymmetricPupilFunction{TF}}
    xs = isodd(wh[1]) ? (((-wh[1]-1)÷2):((wh[1]-1)÷2)) .* Δxy[1] : ((-wh[1]-2)÷2):(wh[1]÷2).*Δxy[1]
    ys = isodd(wh[2]) ? (((-wh[2]-1)÷2):((wh[2]-1)÷2)) .* Δxy[2] : ((-wh[2]-2)÷2):(wh[2]÷2).*Δxy[2]
    return centered([psf(tf, x, y) for x in xs, y in ys])
end

@doc """
    psf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})::OffsetMatrix
    psf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length)::OffsetMatrix

Generate the psf with size `wh` with `Δxy` being the distance between the samples in the x and y dimensions

```jldoctest
julia> tf = BornWolf(488u"nm", 1.7, 1.7)
BornWolf{Float64}(488.0 nm, 1.7, 1.7)
```

```@repl
tf = BornWolf(488u"nm", 1.7, 1.7) # hide
julia> psf(tf, (11,11), 60u"nm")
julia> psf(tf, (5,5), (40u"nm", 50u"nm")) # different pixelsizes in x and y direction
```
"""
@traitfn function psf(tf::TF, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length}) where {TF <: ClosedFormPSFModel; SymmetricPupilFunction{TF}}
    s = all(isodd.(wh)) ? wh : (wh .+ 1)
    buf = centered(Matrix{output_type(tf)}(undef, s...)) # the first quadrant
    # PERF: This can be made even faster if the we use the same calculation for a same radii in one quadrant, for 
    # example (1,2) and (2,1) <12-07-23> 
    buf[1:end, -1:-1:begin] .=
        buf[-1:-1:begin, 1:end] .=
            buf[-1:-1:begin, -1:-1:begin] .=
                buf[1:end, 1:end] .=
                    [psf(tf, x, y) for x in (1:(s[1]-1)÷2) .* Δxy[1], y in (1:(s[2]-1)÷2) .* Δxy[2]]
    # PERF: Could be made faster by using the same for the minimum of `wh[1]` and `wh[2]` <12-07-23> 
    buf[1:end, 0] .= buf[-1:-1:begin, 0] .= psf.(tf, (1:(s[1]-1)÷2) .* Δxy[1], Fill(zero(Δxy[1]), (s[1] - 1) ÷ 2))
    buf[0, 1:end] .= buf[0, -1:-1:begin] .= psf.(tf, Fill(zero(Δxy[2]), (s[2] - 1) ÷ 2), (1:(s[2]-1)÷2) .* Δxy[2])
    buf[0, 0] = psf(tf, zero(Δxy[1]), zero(Δxy[2]))

    if all(isodd.(wh))
        return buf
    else
        @warn "If any of the dimensions of `wh` are even, then the PSF will not be symmetric."
        return centered(buf[begin+1:end, begin+1:end])
    end
end

function psf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})
    tf_otf = otf(tf, wh, Δxy)
    tf_psf = centered(fftshift(ifft(tf_otf)))
    # TODO: Is this correct? How about defocused and other aberrations, can they make the intensity in the center lower?
    # is this even true for a PSF at the focal plane?
    return tf_psf ./ tf_psf[0, 0]
end

psf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length) = psf(tf, wh, (Δxy, Δxy))
psf(tf::TransferFunction, wh::Integer, args...) = psf(tf, (wh, wh), args...)
# FIX: This doesn't strictly speaking make sense, since the PSF is used for convolution and not for term-wise
# multiplication <15-07-23> 
psf(tf::TransferFunction, img::AbstractArray, args...) = psf(tf, size(img), args...)

@doc """
Amplitude point spread function
""" apsf

@doc """
    apsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})::OffsetMatrix{<:Real}
    apsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length)::OffsetMatrix{<:Real}
"""
apsf(tf::TransferFunction, args...; varargs...) = imag.(psf(tf, args...; varargs...))

@doc """
Intensity point spread function
""" ipsf

@doc """
     ipsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})::OffsetMatrix{<:Real}
     ipsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length)::OffsetMatrix{<:Real}
 """
ipsf(tf::TransferFunction, args...; varargs...) = real.(psf(tf, args...; varargs...))

function resolution_limit(tf::TF) where {TF<:TransferFunction}
    # TODO: Is this correct?
    all(hasfield.(TF, [:NA, :λ, :nᵢ])) ? tf.λ / (2 * tf.NA * tf.nᵢ) : nothing
end
