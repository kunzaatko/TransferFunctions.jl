struct SampledPSF{N,PSF<:PointSpreadFunction{N}}
    "physical illumination pattern"
    transfer::PSF
    "pixel dimensions"
    Δxy::PixelSize{N}
    "center coordinate of the sample"
    center::Coordinate{N} # NOTE: It is (0,0) by default <26-08-24> 
end

# FIX: Check if PSF correct has dimension... needs to be in constructor <26-08-24> 

# NOTE: step 1 - fill in the pixel-size
SampledPSF(tf::PointSpreadFunction{N}, Δxy::Length, args...) where {N} = SampledPSF(tf, fillsize(Δxy, N), args...)
# NOTE: step 2 - default center (zeros)
SampledPSF(tf::PointSpreadFunction{N}, Δxy::PixelSize{N}) where {N} = SampledPSF(tf, Δxy, ntuple(_ -> 0, Val(N)))


# FIX: This is ambiguous and should be defined for general dimensions <10-12-23> 
# TODO: There are different approaches to normalization. There is the L∞ constraint that ‖psf‖∞ = 1 and the L1
# constraint that assures the preservation of photometry i.e. ‖psf‖₁ = 1 (sum) <10-12-23> 
# TODO: Implement normalizing to sum to 1 <02-10-23> 
# TODO: Implement expectation of Real valued PSF <02-10-23> 
@traitfn function psf(
    OT::Type{<:Number},
    tf::SampledPSF{N,PSF},
    wh::Tuple{Integer,Integer},
) where {N,PSF;!RadiallySymmetric{PSF}}
    xs = isodd(wh[1]) ? (((-wh[1]-1)÷2):((wh[1]-1)÷2)) .* tf.Δxy[1] : ((-wh[1]-2)÷2):(wh[1]÷2).*tf.Δxy[1]
    ys = isodd(wh[2]) ? (((-wh[2]-1)÷2):((wh[2]-1)÷2)) .* tf.Δxy[2] : ((-wh[2]-2)÷2):(wh[2]÷2).*tf.Δxy[2]
    return centered([psf(OT, tf, x, y) for x in xs, y in ys])
end

# FIX: docs and tests <26-08-24> 
@doc """
    psf(tf::SampledPSF, wh::Tuple{Integer,Integer})::OffsetMatrix

Generate the psf with size `wh` with `Δxy` being the distance between the samples in the x and y dimensions

```jldoctest
julia> tf = BornWolf(488u"nm", 1.7, 1.7)
BornWolf{Float64}(488.0 nm, 1.7, 1.7)
```

```@repl
tf = BornWolf(488u"nm", 1.7, 1.7) # hide
psf(tf, (11,11), 60u"nm")
psf(tf, (5,5), (40u"nm", 50u"nm")) # different pixelsizes in x and y direction
```
"""
@traitfn function psf(
    OT::Type{<:Number},
    s_tf::SampledPSF{N,PSF},
    wh::Tuple{Integer,Integer}
) where {N,PSF<:ModelPSF;RadiallySymmetric{PSF}}
    # TODO: Refactor. Make symmetric optimized array generation into its self function <02-10-23> 
    s = all(isodd.(wh)) ? wh : (wh .+ 1)
    buf = centered(Matrix{OT}(undef, s...)) # the first quadrant
    # # PERF: This can be made even faster if the we use the same calculation for a same radii in one quadrant, for 
    # # example (1,2) and (2,1) <12-07-23> 
    # buf[1:end, -1:-1:begin] .=
    #     buf[-1:-1:begin, 1:end] .=
    #         buf[-1:-1:begin, -1:-1:begin] .=
    #             buf[1:end, 1:end] .=
    #                 [psf(tf, x, y) for x in (1:(s[1]-1)÷2) .* Δxy[1], y in (1:(s[2]-1)÷2) .* Δxy[2]]
    # # PERF: Could be made faster by using the same for the minimum of `wh[1]` and `wh[2]` <12-07-23> 
    # buf[1:end, 0] .= buf[-1:-1:begin, 0] .= psf.(tf, (1:(s[1]-1)÷2) .* Δxy[1], Fill(zero(Δxy[1]), (s[1] - 1) ÷ 2))
    # buf[0, 1:end] .= buf[0, -1:-1:begin] .= psf.(tf, Fill(zero(Δxy[2]), (s[2] - 1) ÷ 2), (1:(s[2]-1)÷2) .* Δxy[2])
    # buf[0, 0] = psf(tf, zero(Δxy[1]), zero(Δxy[2]))

    # NOTE: This could be made by concatenating 4 symmetric matrices, that represent a quadrant (optimizes copying and
    # generating #ops) <02-10-23> 

    # NOTE: `IndirectArrays` can be used to make as few computations as necessary <02-10-23, kunzaatko> 
    # https://github.com/JuliaArrays/IndirectArrays.jl

    # TODO: There is probably a method to do this even more efficiently. In general, there, will be 8 pixels with the
    # same distance from origin... Compare it to running on the grid and the goal should be 8× faster <02-10-23> 

    # PERF: This method is comparable to the one above. There is a possibility to optimize it using broadcasting, 
    # probably. <02-10-23> 
    cache = Dict{Length,OT}()
    map!(buf, Tuple.(CartesianIndices(buf))) do (x, y)
        r = hypot(x * s_tf.Δxy[1], y * s_tf.Δxy[2])
        if r ∉ keys(cache)
            cache[r] = psf(OT, s_tf.transfer, x * s_tf.Δxy[1], y * s_tf.Δxy[2])
        end
        cache[r]
    end

    # PERF: No symmetry optimization (only for comparison with optimizations) <02-10-23> 
    # map!(buf, Tuple.(CartesianIndices(buf))) do (x, y)
    #     psf(tf, x * Δxy[1], y * Δxy[2])
    # end

    if all(isodd.(wh))
        return buf
    else
        @warn "If any of the dimensions of `wh` are even, then the PSF will not be symmetric."
        return centered(buf[begin+1:end, begin+1:end])
    end
end

function psf(
    tf::SampledPSF{N,PSF},
    wh::Tuple{Integer,Integer}
) where {N,PSF}
    return psf(preferred_type(PSF), tf, wh)
end

# FIX: Incorrect API architecture <26-08-24> 
# NOTE: Has to be defined for N-dims generally and not specific dimensions because otherwise, there could be ambiguity 
# if a concrete type implements generic N-dim `psf` method <10-12-23> 
function psf(
    OT::Type{<:Number},
    tf::TransferFunction{N},
    wh::NTuple{N,Integer},
    Δxy::NTuple{N,Length}
) where {N}
    tf_otf = otf(tf, wh, Δxy)
    tf_psf = centered(fftshift(ifft(tf_otf)))
    # TODO: Is this correct? How about defocused and other aberrations, can they make the intensity in the center lower?
    # is this even true for a PSF at the focal plane?
    # FIX: This is not correct for a Ndim PSF only for 2D <10-12-23> 
    return tf_psf ./ sum(tf_psf)
end

# psf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length; vargs...) = psf(tf, wh, (Δxy, Δxy); vargs...)
# psf(tf::TransferFunction, wh::Integer, args...; vargs...) = psf(tf, (wh, wh), args...; vargs...)
# FIX: This doesn't strictly speaking make sense, since the PSF is used for convolution and not for term-wise
# multiplication <15-07-23> 
psf(tf::SampledPSF{N,PSF}, img::AbstractArray{OT,N}, args...; vargs...) where {OT,PSF,N} = psf(OT, tf, size(img), args...; vargs...)

# FIX: This is not the correct name!! <30-11-23> 
@doc """
Amplitude point spread function
""" apsf

@doc """
    apsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})::OffsetMatrix{<:Real}
    apsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length)::OffsetMatrix{<:Real}
"""
apsf(OT::Type{<:Real}, tf::SampledPSF, args...; vargs...) = imag(psf(Complex{OT}, tf, args...; vargs...))
apsf(tf::SampledPSF{N,PSF}, args...; vargs...) where {N,PSF} = apsf(preferred_type(PSF), tf, args...; vargs...)

@doc """
Intensity point spread function
""" ipsf

@doc """
     ipsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})::OffsetMatrix{<:Real}
     ipsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length)::OffsetMatrix{<:Real}
 """
ipsf(OT::Type{<:Real}, tf::SampledPSF, args...; vargs...) = real(psf(Complex{OT}, tf, args...; vargs...))
ipsf(tf::SampledPSF{N,PSF}, args...; vargs...) where {N,PSF} = ipsf(preferred_type(PSF), tf, args...; vargs...)

# psf.jl
# psf(tfr::SampledPSF, args...; varargs...) = psf(tfr.transfer, args..., tfr.Δxy; varargs...)
# apsf(tfr::SampledPSF, args...; varargs...) = apsf(tfr.transfer, args..., tfr.Δxy; varargs...)
# ipsf(tfr::SampledPSF, args...; varargs...) = ipsf(tfr.transfer, args..., tfr.Δxy; varargs...)

# TODO: Better Base.show <30-11-23> 
