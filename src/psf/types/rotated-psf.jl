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
  RotatedPSF{N,R,P}(parent::P, rotation::R) where {N,R,P} = new{N,R,P}(parent, rotation)
end
RotatedPSF{N,R,P}(params::ComponentVector) where {N, R, P} = RotatedPSF{N,R,P}(P(params.parent), R(params.rot...))

RotatedPSF(P::Type{<:PointSpreadFunction{3}}, args...; α=0.0, β=0.0, γ=0.0, kwargs...)  = RotatedPSF(P(args...; kwargs...), RotXYZ(α, β, γ)) # final 2
RotatedPSF(P::Type{<:PointSpreadFunction{2}}, args...; θ=0.0, kwargs...)  = RotatedPSF(P(args...; kwargs...), Angle2d(θ)) # final 3
RotatedPSF{N}(args...; kwargs...) where {N} = RotatedPSF(args...; kwargs...)

@inline Base.parent(tf::RotatedPSF) = tf.parent

@inline function response(tf::RotatedPSF{3}, x_r::Length, y_r::Length, z_r::Length)
  x_r, y_r, z_r = tf.rotation * SVector(x_r, y_r, z_r)
  response(parent(tf), x_r, y_r, z_r)
end

@inline function response(tf::RotatedPSF{2}, x_r::Length, y_r::Length)
  x_r, y_r = tf.rotation * SVector(x_r, y_r)
  response(parent(tf), x_r, y_r)
end

params(tf::RotatedPSF) = ComponentVector(rot=Rotations.params(tf.rotation), parent=params(parent(tf)))

export RotatedPSF
