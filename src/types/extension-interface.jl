abstract type ExtensionMethod end
@require_interface extend(m::ExtensionMethod, ::Coordinate{2})

abstract type InterpolateExtrapolate end
@require_interface extend(m::InterpolateExtrapolate, ::Coordinate{2})
@require_interface inconvexhull(m::InterpolateExtrapolate, ::Coordinate{2})
@require_interface interpolate(m::InterpolateExtrapolate, ::Coordinate{2})
@require_interface extrapolate(m::InterpolateExtrapolate, ::Coordinate{2})
