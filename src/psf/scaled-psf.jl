const Scale{N,T} = NTuple{N, T} where {T<:Real}
struct ScaledPSF{N,T,P<:PointSpreadFunction{N}} <: PointSpreadFunction{N}
    parent::P
    scale::Scale{N,T}
    function ScaledPSF{N, T}(parent::PointSpreadFunction{N}, scale::Scale{N,T}) where {N,T} 
        all(scale .> zero.(scale)) || throw(ArgumentError("Scale must be positive. Got `scale = $(scale)`"))
        return new{N,T,typeof(parent)}(parent, scale)
    end
end
ScaledPSF(P::Type{<:PointSpreadFunction{3}}, args...; xscale=1.0, yscale=1.0, zscale=1.0, kwargs...)  = ScaledPSF(P(args...; kwargs...), (xscale, yscale, zscale))
ScaledPSF(P::Type{<:PointSpreadFunction{2}}, args...; xscale=1.0, yscale=1.0, kwargs...)  = ScaledPSF(P(args...; kwargs...), (xscale, yscale))
ScaledPSF{N}(args...; kwargs...) where {N} = ScaledPSF(args...; kwargs...)
function ScaledPSF(parent::PointSpreadFunction{N}, scale::Scale{N}) where {N}
    scale = promote(scale...)
    ScaledPSF{N,eltype(scale)}(parent, scale)
end

@inline Base.parent(tf::ScaledPSF) = tf.parent

@inline function response(tf::ScaledPSF{3}, x_r::Length, y_r::Length, z_r::Length=0.0u"nm")
    x_r, y_r, z_r = SVector(x_r, y_r, z_r) ./ tf.scale
  response(parent(tf), x_r, y_r, z_r)
end

@inline function response(tf::ScaledPSF{2}, x_r::Length, y_r::Length)
    x_r, y_r = SVector(x_r, y_r) ./ tf.scale
  response(parent(tf), x_r, y_r)
end

export ScaledPSF
