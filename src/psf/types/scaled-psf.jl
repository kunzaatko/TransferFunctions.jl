using StaticArraysCore

const Scale{N,T} = SVector{N, T} where {T<:Real}

# TODO: Finish the docs with examples <20-08-25> 
"""
    ScaledPSF{N,T,P<:PointSpreadFunction{N}} <: PointSpreadFunction{N}
A wrapper around the parent `P <: PointSpreadFunction{N}` that scales the PSF with `scale`.
"""
struct ScaledPSF{N,T,P<:PointSpreadFunction{N}} <: PointSpreadFunction{N}
    parent::P
    scale::Scale{N,T}
    function ScaledPSF{N,T,P}(parent::P, scale::Scale{N,T}) where {N,T,P}
        all(scale .> zero.(scale)) || throw(ArgumentError("Scale must be positive. Got `scale = $(scale)`"))
        new{N,T,P}(parent, scale)
    end
end
ScaledPSF{N,T,P}(params::ComponentVector) where {N, T, P} = ScaledPSF{N,T,P}(P(params.parent), SVector{N}(T.(params.scale)))

ScaledPSF{N, T}(parent::PointSpreadFunction{N}, scale::Scale{N,T}) where {N,T} = ScaledPSF{N,T,typeof(parent)}(parent, scale) # final 2
function ScaledPSF(parent::PointSpreadFunction{N}, scale::Scale{N}) where {N} # partial 1 -> final 2
    scale = SVector{N}(promote(scale...))
    ScaledPSF{N,eltype(scale)}(parent, scale)
end
ScaledPSF(P::Type{<:PointSpreadFunction{3}}, args...; xscale=1.0, yscale=1.0, zscale=1.0, kwargs...)  = ScaledPSF(P(args...; kwargs...), SVector{3}([xscale, yscale, zscale])) # -> partial 1
ScaledPSF(P::Type{<:PointSpreadFunction{2}}, args...; xscale=1.0, yscale=1.0, kwargs...)  = ScaledPSF(P(args...; kwargs...), SVector{2}([xscale, yscale])) # -> partial 1
ScaledPSF{N}(args...; kwargs...) where {N} = ScaledPSF(args...; kwargs...)

@inline Base.parent(tf::ScaledPSF) = tf.parent

@inline function response(tf::ScaledPSF{3}, x_r::Length, y_r::Length, z_r::Length)
    x_r, y_r, z_r = SVector(x_r, y_r, z_r) ./ tf.scale
  response(parent(tf), x_r, y_r, z_r)
end

@inline function response(tf::ScaledPSF{2}, x_r::Length, y_r::Length)
    x_r, y_r = SVector(x_r, y_r) ./ tf.scale
  response(parent(tf), x_r, y_r)
end

params(tf::ScaledPSF) = ComponentVector(scale=tf.scale, parent=params(parent(tf)))

export ScaledPSF
