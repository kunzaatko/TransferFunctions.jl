"""
    OTFArray{T<:Number,N} <: OpticalTransferFunction{N}

# Fields
- `data::AbstractMatrix{T}`
- `center::Coordinate{2}`
"""
struct OTFArray{T<:Number,N, AA<:SampledArray{T, <:Frequency, N}} <: OpticalTransferFunction{N}
    data::AA
    origin::Coordinate{N}
    function OTFArray(data::AA, origin::Coordinate{N}) where {N,AA}
        contained(data, origin) || throw(DomainError(origin, "The center is not within the data bounds: $origin ∉ $(axes(data))"))
        new{eltype(data),N, AA}(data, origin)
    end
end
# OTFArray(data::AbstractMatrix{T}, Δxy::Length, args...) where {T} = OTFArray(data, fillsize(Δxy, 2), args...)
# OTFArray(data::AbstractMatrix{T}, Δxy::PixelSize{2}) where {T} = OTFArray(data, Δxy, exactcenter(data))

# function Base.show(io::IO, ::MIME"text/plain", tf::OTFArray{T}) where {T}
#     showcenter = tf.origin == tf.data .÷ 2
#     centerstring = showcenter ? ", center = $(tf.origin)" : ""
#     print(io, "OTFArray(Δxy = $(allequal(tf.Δxy) ? tf.Δxy[1] : tf.Δxy)$(centerstring)) with eltype $T with $(join(map(string, size(tf.data)), "×")) points:\n")
#     Base.print_array(io, tf.data)
# end

export OTFArray
