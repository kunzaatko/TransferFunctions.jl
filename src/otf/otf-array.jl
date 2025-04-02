"""
    OTFArray{T<:Number} <: OpticalTransferFunction

# Fields
- `data::AbstractMatrix{T}`
- `Δxy::PixelSize{2}`
- `center::Coordinate{2}`
- `extend::ExtensionMethod`
"""
struct OTFArray{T<:Number} <: OpticalTransferFunction
    data::AbstractMatrix{T}
    Δxy::PixelSize{2}
    origin::Coordinate{2}
    extend::ExtensionMethod
    function OTFArray(data::AbstractMatrix{T}, Δxy::PixelSize{2}, origin::Coordinate{2}) where {T<:Number}
        contained(data, origin) || throw(DomainError(origin, "The center is not within the data bounds: $origin ∉ $(axes(data))"))
        new{T}(data, Δxy, origin)
    end
end
OTFArray(data::AbstractMatrix{T}, Δxy::Length, args...) where {T} = OTFArray(data, fillsize(Δxy, 2), args...)
OTFArray(data::AbstractMatrix{T}, Δxy::PixelSize{2}) where {T} = OTFArray(data, Δxy, exactcenter(data))

function Base.show(io::IO, ::MIME"text/plain", tf::OTFArray{T}) where {T}
    showcenter = tf.origin == tf.data .÷ 2
    centerstring = showcenter ? ", center = $(tf.origin)" : ""
    print(io, "OTFArray(Δxy = $(allequal(tf.Δxy) ? tf.Δxy[1] : tf.Δxy)$(centerstring)) with eltype $T with $(join(map(string, size(tf.data)), "×")) points:\n")
    Base.print_array(io, tf.data)
enD
