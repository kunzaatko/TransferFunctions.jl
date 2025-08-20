"""
    PSFArray{T<:Real, N} <: PointSpreadFunction{N}
A PSF that is defined by an array of values sampled at some rate. 

The data is stored in the `data` field of the type [`SpatialArray`](@ref) which has to have a dimension of `N`.
"""
struct PSFArray{T, N, P<:SpatialArray{T, N}} <: PointSpreadFunction{N}
    data::P
    function PSFArray(data::SpatialArray)
        new{eltype(data),ndims(data),typeof(data)}(data)
    end
end

export PSFArray
