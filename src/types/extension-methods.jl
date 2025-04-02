abstract type ExtensionMethod end
interpolate(m::ExtensionMethod, ::Coordinate{2}) = no_implemementation_error(typeof(m), :interpolate)
extrapolate(m::ExtensionMethod, ::Coordinate{2}) = no_implemementation_error(typeof(m), :interpolate)

abstract type InterpolateExtrapolate end
extend(m::InterpolateExtrapolate, ::Coordinate{2}) = no_implemementation_error(typeof(m), :extend)
interpolate(m::InterpolateExtrapolate, c::Coordinate{2}) = extend(m, c)
extrapolate(m::InterpolateExtrapolate, c::Coordinate{2}) = extend(m, c)
