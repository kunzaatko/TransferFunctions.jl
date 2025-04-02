"""
    SampledImage{T,N} <: AbstractArray{T,N}
An `N`-dimensional sampled image with data of type `T`and the axis dimensions.

!!! note
    The default implementations rely on the data being stored in the `data` field.
"""
abstract type SampledImage{T,N} <: AbstractArray{T,N} end
Base.size(img::SampledImage) = Base.size(img.data)
Base.getindex(img::SampledImage, args...) = Base.getindex(img.data, args...)

"""
    SpatialImage{T,N} <: SampledImage{T,N} 
An `N`-dimensional sampled image with 
"""
struct SpatialImage{T,N} <: SampledImage{T,N}
    data::AbstractArray{T,N}
    Δxy::PixelSize{N}
end
SpatialImage(data::AbstractArray, Δxy::Length) = SpatialImage(data, fillsize(Δxy, ndims(data)))
