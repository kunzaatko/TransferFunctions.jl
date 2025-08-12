using InterfaceFunctions

"""
    TransferFunction
Super type for all transfer functions

A transfer function in microscopy is an object specifying the response in the image plane to the light coming from the
object plane.
There are several characteristics of a transfer functions which can influence the character of the response. Most
commonly we are concerned with its spatial variance/invariance, i.e. whether the response is dependent of the objects
position in the object plane, and linearity (which is almost always satisfied).

See also [`LinearShiftInvariantTransferFunction`](@ref), [`ImpulseResponseMapping`](@ref)
"""
abstract type TransferFunction end
@interface transfer(t::TransferFunction, ::SpatialMatrix{<:Real})
@interface restore(t::TransferFunction, ::SpatialMatrix{<:Real})

# NOTE: Taken from Distributions.jl <kunzaatko> 
for func in (:(==), :isequal, :isapprox)
    @eval function Base.$func(tf1::A, tf2::B; kwargs...) where {A<:TransferFunction,B<:TransferFunction}
        nameof(A) === nameof(B) || return false
        fields = fieldnames(A)
        fields === fieldnames(B) || return false

        for f in fields
            isdefined(tf1, f) && isdefined(tf2, f) || return false
            # perform equivalence check to support types that have no defined equality, such
            # as `missing`
            getfield(tf1, f) === getfield(tf2, f) || $func(getfield(tf1, f), getfield(tf2, f); kwargs...) || return false
        end

        return true
    end
end

# NOTE: Taken from Distributions.jl <kunzaatko> 
function Base.hash(tf::TransferFunction, h::UInt)
    hashed = hash(TransferFunction, h)
    hashed = hash(nameof(typeof(tf)), hashed)

    for f in fieldnames(typeof(tf))
        hashed = hash(getfield(tf, f), hashed)
    end

    return hashed
end

export transfer, restore
