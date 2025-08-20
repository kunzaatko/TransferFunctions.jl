using Rotations, ComponentArrays, StaticArraysCore

"""
    RotatedPSF{N,T,P} <: PointSpreadFunction{N}
A wrapper around the parent `P <: PointSpreadFunction{N}` that rotates the PSF with [rotation](@extref Rotations
Rotation-Types).
"""
struct RotatedPSF{N,R<:Rotation{N},P<:PointSpreadFunction{N}} <: PointSpreadFunction{N}
  parent::P
  rotation::R
  RotatedPSF(parent::PointSpreadFunction{N}, rotation::Rotation{N}) where {N} = new{N,typeof(rotation),typeof(parent)}(parent, rotation)
end
RotatedPSF(P::Type{<:PointSpreadFunction{3}}, args...; α=0.0, β=0.0, γ=0.0, kwargs...)  = RotatedPSF(P(args...; kwargs...), RotXYZ(α, β, γ))
RotatedPSF(P::Type{<:PointSpreadFunction{2}}, args...; θ=0.0, kwargs...)  = RotatedPSF(P(args...; kwargs...), RotMatrix{2}(θ))
RotatedPSF{N}(args...; kwargs...) where {N} = RotatedPSF(args...; kwargs...)

@inline Base.parent(tf::RotatedPSF) = tf.parent

@inline function response(tf::RotatedPSF{3}, x_r::Length, y_r::Length, z_r::Length=0.0u"nm")
  x_r, y_r, z_r = tf.rotation * SVector(x_r, y_r, z_r)
  response(parent(tf), x_r, y_r, z_r)
end

@inline function response(tf::RotatedPSF{2}, x_r::Length, y_r::Length)
  x_r, y_r = tf.rotation * SVector(x_r, y_r)
  response(parent(tf), x_r, y_r)
end

export RotatedPSF
