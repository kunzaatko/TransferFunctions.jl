"""
    PSFArray{T<:Real} <: PointSpreadFunction
A measurement of a PSF from the acquired image(s).

# Fields
- `data::AbstractMatrix{T}`
- `Δ::PixelSize{2}`
- `center::Coordinate{2}`
- `extend::ExtensionMethod`
"""
struct PSFArray{T<:Real} <: MeasuredPSF
    data::AbstractMatrix{T}
    Δ::PixelSize{2}
    center::Coordinate{2}
    extend::ExtensionMethod
    function PSFArray(data::AbstractMatrix{T}, Δ::PixelSize{2}, center::Coordinate{2}) where {T<:Real}
        contained(data, center) || throw(DomainError(center, "The center is not within the data bounds: $center ∉ $(axes(data))"))
        new{T}(data, Δ, center)
    end
end

PSFArray(data::AbstractMatrix{T}, Δ::Length, args...) where {T} = PSFArray(data, fillsize(Δ, 2), args...)
PSFArray(data::AbstractMatrix{T}, Δ::PixelSize{2}) where {T} = PSFArray(data, Δ, exactcenter(data))

function Base.show(io::IO, ::MIME"text/plain", tf::PSFArray{T}) where {T}
    showcenter = tf.center == tf.data .÷ 2
    centerstring = showcenter ? ", center = $(tf.center)" : ""
    print(io, "PSFArray(Δxy = $(allequal(tf.Δ) ? tf.Δ[1] : tf.Δ)$(centerstring)) with eltype $T with $(join(map(string, size(tf.data)), "×")) points:\n")
    Base.print_array(io, tf.data)
end
