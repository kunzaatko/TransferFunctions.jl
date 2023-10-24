using Interpolations, OffsetArrays
@doc raw"""
`MeasuredPSF` holds an array of measured data with information about the dimensions of the measurement. This allows
you to make conversions to other representations of a transfer function such as an OTF and to use interpolations to
sample the PSF at arbitrary locations not included in the initial measurement (in general to use the measurements 
similarly to how one could use a model). This can be useful for using the same PSF measurement for an acquisition with 
a different pixelsize or in super resolution applications.

A measurement can be done for example using an image of subresolution microspheres.
"""
struct MeasuredPSF{T<:Real,R<:Real} <: MeasuredTransferFunction
    # TODO: This should be an OffsetArray instead and a conversion should be made if it is not in order to be able to
    # sample at arbitrary locations from the center. Maybe convert by using the maximum pixel or store the center in the
    # struct computed as the maximum of some interpolation? <15-09-23> 
    # TODO: Make a trait applies to a circularly symmetric PSF... That is to assume that the data is symmetric and the
    # non-compliance is caused only by noise. That can then be used for averaging the measurement if a `r::Length` is
    # supplied to the sampling procedure <15-09-23> 
    "array of the measured PSF"
    data::AbstractMatrix{T} # TODO: General dimension array <15-09-23> 
    "dimensions of the `data` array"
    Δxy::Tuple{Length,Length}
    "center of the PSF measurement"
    center::Tuple{R,R}
    # FIX: Check center in-bounds <15-09-23> 
end

MeasuredPSF(data, Δxy::Length, args...) = MeasuredPSF(data, (Δxy, Δxy), args...)
# INFO: infer the center to be the center of the array if missing (non-integer if even array)
MeasuredPSF(data, Δxy::Tuple{Length,Length}) = MeasuredPSF(data, Δxy,)

# TODO: Test <21-09-23> 
@doc raw"""
    psf(tf::MeasuredPSF, wh::Tuple{Integer,Integer} [,Δxy::Tuple{Length,Length}]; intp, extp)::OffsetMatrix

Returns a **centered** PSF interpolated from the measurement.

# Parameters
+ `intp = BSpline(Cubic(Flat(OnGrid())))`: interpolation method (default: cubic interpolation with flat slopes on the edge knots)
+ `extp = zero(eltype(tf.data))`: extrapolation method
  - `zero(eltype(tf.data))` constant extrapolation outside of the bounds (ideal for background subtracted PSF)
  - `Flat()` constant slope beyond the bounds (gives the last knot value and is equivalent to `0` for a background subtracted PSF)
+ `δ = wh .÷ 2`: output center of the PSF (default: center of the output size array rounded down)

# Extended Help
## Implementation Details
+ Interpolation is done using [scaled BSplines](https://juliamath.github.io/Interpolations.jl/latest/control/#Scaled-BSplines).
"""
psf(tf::MeasuredPSF, wh::Tuple{Integer,Integer}) = psf(tf, wh, tf.Δxy)
function psf(tf::MeasuredPSF, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length};
    intp=BSpline(Cubic(Flat(OnGrid()))),
    extp=zero(eltype(tf.data)),
    δ=(wh .÷ 2))

    in_x, in_y = (range(tf.Δxy[i], size(tf.data, i) * tf.Δxy[i]; length=size(tf.data, i)) .- (tf.center[i] + 1) * tf.Δxy[i] for i in 1:2)
    out_x, out_y = (range(Δxy[i], wh[i] * Δxy[i]; length=wh[i]) .- (δ[i] + 1) * Δxy[i] for i in 1:2)

    in_intp = interpolate(tf.data, intp)
    in_intp_s = scale(in_intp, in_x, in_y)
    # NOTE: `extp=0` for the extrapolation means that every out of grid call will be constant 0  
    in_extp_s = extrapolate(in_intp_s, extp)

    # TODO: Determine what type should the output have if we have a non-integer center <21-09-23> 
    @assert all(isinteger.(δ)) "Currently only integer centers of the output are accepted"

    # FIX: The units in this computation make a mess. they should be handled outside the extrapolation <22-09-23> 
    return OffsetMatrix([in_extp_s(x, y) for x in out_x, y in out_y], δ .* (-1))
end

function Base.show(tf::MeasuredPSF)
    resolution = tf.Δxy[1] == tf.Δxy[2] ? "$(tf.Δxy[1])" : "$(tf.Δxy[1])×$(tf.Δxy[2])"
    center = tf.center == tf.data
    t = "MeasuredPSF(Δxy=$resolution, center=$center)"
    return t
end

center(A::AbstractMatrix) = map(ax -> (ax[2] - ax[1]) / 2, extrema.(axes(A)))

# TODO: Implement better `show` for `MeasuredPSF`... The data should be shown outside the dimensions parameters <15-09-23> 
