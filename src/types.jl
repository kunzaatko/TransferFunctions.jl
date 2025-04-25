include("types/sampled-arrays.jl")
include("types/extension-methods.jl")
include("types/circulant-tensor.jl")

"""
    TransferFunction
""" # TODO: Docs <24-04-25> 
abstract type TransferFunction end
@require_interface transfer(t::TransferFunction, ::SpatialArray{<:Real,2})
@require_interface restore(t::TransferFunction, ::SpatialArray{<:Real,2})

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
