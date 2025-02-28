# TODO: Reformulate the documentation: This allows you to make conversions to other representations of a transfer
# function such as a PSF and to use interpolations to sample the OTF at arbitrary locations not included in the initial
# measurement (in general to use the measurements similarly to how one could use a model). This can be useful for using
# the same OTF measurement for an acquisition with a different pixel-size or in super resolution applications.

# TODO: Should instead hold a sampled array or instead some derived type of a sampled array that has a frequencies axis <29-12-24> 

# TODO: Are OTFs really measurable? Shouldn't it be only for the PSF? Maybe a MeasuredOTF is created as a measured PSF? Perhaps from the deconvolution (inversion) of the ground truth model problem and with the transferred data (H = I_gt/I_raw) <26-08-24> 
@doc raw"""
    MeasuredOTF{T<:Number, N} <: OpticalTransferFunction{N}
A measurement of an OTF from the acquired image(s). 

# Fields
- `data::AbstractArray{T,N}`
- `Δxy::PixelSize{N}`
- `center::Coordinate{N}`
"""
struct MeasuredOTF{T<:Number,N} <: OpticalTransferFunction{N}
    data::AbstractArray{T,N}
    Δxy::PixelSize{N}
    center::Coordinate{N}
    function MeasuredOTF(data::AbstractArray{T,N}, Δxy::PixelSize{N}, center::Coordinate{N}) where {N,T<:Number}
        contained(data, center) || throw(DomainError(center, "The center is not within the data bounds: $center ∉ $(axes(data))"))
        new{T,N}(data, Δxy, center)
    end
end

# TODO: Should infer the center if an OffsetArray that contains the origin is used in the constructor <29-12-24> 

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
