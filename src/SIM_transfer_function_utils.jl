using Interpolations

# IDEA!: Shift can function differently based on the type it gets to shift. e.g. "new" OTFEvalution (planned) should 
# handle the ifftshift and shift back and could have different defaults for the boundary conditions of the interpolation
# <24-10-23> 
# TODO: Should accept more dimensions <22-09-23> 
function shift(A::AbstractMatrix, δ::Tuple{T,T};
    intp=BSpline(Cubic(Flat(OnGrid()))),
    extp=zero(eltype(A))
    # IDEA: Add scale argument (as in interpolations... That is add argument for axes of A) <22-09-23> 
) where {T<:Real}
    # TODO: Check if bounds make sense for interpolation (if they are in the bounds of the data) <22-09-23> 

    A_intp = interpolate(A, intp)
    A_extp = extrapolate(A_intp, extp)

    out_x, out_y = (axes(A, i) .- δ[i] for i in 1:2)
    return A_extp(out_x, out_y)
end
