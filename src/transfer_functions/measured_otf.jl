@doc raw"""
`MeasuredOTF` holds an array of measured data with information about the dimensions of the measurement. This allows
you to make conversions to other representations of a transfer function such as a PSF and to use interpolations to
sample the OTF at arbitrary locations not included in the initial measurement (in general to use the measurements 
similarly to how one could use a model). This can be useful for using the same OTF measurement for an acquisition with 
a different pixelsize or in super resolution applications.
"""
struct MeasuredOTF{T<:Real,N} <: MeasuredTransferFunction{N}
    "array of the measured OTF"
    data::AbstractArray{T,N}
    "dimensions of the `data` array"
    Δxy::NTuple{N,Length}
    "center of the PSF measurement"
    center::NTuple{N,Real}
    function MeasuredOTF(data, Δxy, center)
        for (dim, (lims, c)) in enumerate(zip(extrema.(axes(data)), center))
            lims[1] <= c <= lims[2] || throw(DomainError(center, "The center is not within the data bounds for dimension $dim (axes(data, $dim) =  $(axes(data, dim)))"))
        end
        new{eltype(data),ndims(data)}(data, Δxy, center)
    end
end

MeasuredOTF(data, Δxy::Length, args...) = MeasuredOTF(data, tuple(fill(Δxy, ndims(data))...), args...)

# TODO: Add note to documentation that the centre inference prefers integer pixel values <28-11-23> 
# INFO: infer the center to be the center of the array if missing (non-integer if even array)
MeasuredOTF(data, Δxy::NTuple) = MeasuredOTF(data, Δxy, size(data) .÷ 2)

function Base.show(io::IO, ::MIME"text/plain", tf::MeasuredOTF{<:Real,N}) where {N}
    showcenter = tf.center == tf.data .÷ 2
    centerstring = showcenter ? ", center = $(tf.center)" : ""
    print(io, "MeasuredOTF{$N}(Δxy = $(allequal(tf.Δxy) ? tf.Δxy[1] : tf.Δxy)$(centerstring)) with eltype $(eltype(tf.data)) with $(join(map(string, size(tf.data)), "×")) points:\n")
    Base.print_array(io, tf.data)
end
