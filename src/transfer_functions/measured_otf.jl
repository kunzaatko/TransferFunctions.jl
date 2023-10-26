@doc raw"""
`MeasuredOTF` holds an array of measured data with information about the dimensions of the measurement. This allows
you to make conversions to other representations of a transfer function such as a PSF and to use interpolations to
sample the OTF at arbitrary locations not included in the initial measurement (in general to use the measurements 
similarly to how one could use a model). This can be useful for using the same OTF measurement for an acquisition with 
a different pixelsize or in super resolution applications.
"""
struct MeasuredOTF{T<:Real,R<:Real} <: MeasuredTransferFunction
    "array of the measured OTF"
    data::AbstractMatrix{T} # TODO: General dimension array <15-09-23> 
    "dimensions of the `data` array"
    Δxy::Tuple{Length,Length}
    "center of the PSF measurement"

    center::Tuple{R,R}
    # FIX: Check center in-bounds <15-09-23> 
end

