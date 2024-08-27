# TODO: Are OTFs really measurable? Shouldn't it be only for the PSF? Maybe a MeasuredOTF is created as a measured PSF? <26-08-24> 
@doc raw"""
`MeasuredOTF` holds an array of measured data with information about the dimensions of the measurement. This allows
you to make conversions to other representations of a transfer function such as a PSF and to use interpolations to
sample the OTF at arbitrary locations not included in the initial measurement (in general to use the measurements 
similarly to how one could use a model). This can be useful for using the same OTF measurement for an acquisition with 
a different pixelsize or in super resolution applications.
"""
struct MeasuredOTF{T<:Number,N}
    "array of the measured OTF"
    data::AbstractArray{T,N}
    "dimensions of the `data` array"
    Δxy::PixelSize{N}
    "center of the OTF measurement"
    center::Coordinate{N} # FIX: Is this necessary? <26-08-24> 
    function MeasuredOTF(data::AbstractArray{T,N}, Δxy::PixelSize{N}, center::Coordinate{N}) where {N,T<:Number}
        contained(data, center) || throw(DomainError(center, "The center is not within the data bounds: $center ∉ $(axes(data))"))
        new{T,N}(data, Δxy, center)
    end
end

MeasuredOTF(data::AbstractArray{T,N}, Δxy::Length, args...) where {T,N} = MeasuredOTF(data, fillsize(Δxy, N), args...)

# TODO: Add note to documentation that the centre inference prefers integer pixel values <28-11-23> 
# INFO: infer the center to be the center of the array if missing (non-integer if even array)
# TODO: Add info about the selected center. It may not be correctly placed. <26-08-24> 
MeasuredOTF(data::AbstractArray{T,N}, Δxy::PixelSize{N}) where {T,N} = MeasuredOTF(data, Δxy, exactcenter(data))

function Base.show(io::IO, ::MIME"text/plain", tf::MeasuredOTF{T,N}) where {T,N}
    showcenter = tf.center == tf.data .÷ 2
    centerstring = showcenter ? ", center = $(tf.center)" : ""
    print(io, "MeasuredOTF{$N}(Δxy = $(allequal(tf.Δxy) ? tf.Δxy[1] : tf.Δxy)$(centerstring)) with eltype $T with $(join(map(string, size(tf.data)), "×")) points:\n")
    Base.print_array(io, tf.data)
end
