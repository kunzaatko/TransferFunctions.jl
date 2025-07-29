using InterfaceFunctions

abstract type ExtensionMethod end
@interface extend(m::ExtensionMethod, ::Coordinate{2})

abstract type InterpolateExtrapolate end
@interface extend(m::InterpolateExtrapolate, ::Coordinate{2})
@interface inconvexhull(m::InterpolateExtrapolate, ::Coordinate{2})
@interface interpolate(m::InterpolateExtrapolate, ::Coordinate{2})
@interface extrapolate(m::InterpolateExtrapolate, ::Coordinate{2})
