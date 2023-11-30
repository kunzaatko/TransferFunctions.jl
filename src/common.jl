# TODO: Documentation <28-11-23> 
abstract type TransferFunction{N} end
const TF{N} = TransferFunction{N}

# FIX: Is this useful? Shouldn't it be implemented default dispatch as in IlluminationPatterns <28-11-23> 
output_type(::TransferFunction{N}) where {N} = Complex{Float64}

# TODO: Update docs <24-10-23> 
@doc raw"""
An abstract type for any transfer function that is based on a physical model of an optical system. Contrary to a 
[`MeasuredTransferFunction`](@ref), a `ModelTransferFunction` must be quantifiable at any point (in either spatial or 
frequency domain). The model can be of an OTF (`ClosedFormOTFModel <: ModelTransferFunction`) or a PSF 
(`ClosedFormPSFModel <: ModelTransferFunction`) or a PupilFunction.
"""
abstract type ModelTransferFunction{N} <: TransferFunction{N} end

# NOTE: Allows broadcasting `func.(tf::ModelTransferFunction, a:b)`where `func` can be any of `psf`,`otf`,`mtf`, etc.
Broadcast.broadcastable(tf::ModelTransferFunction) = Ref(tf)

# FIX: This should be instead type branched as abstract OTF and abstract PSF?! <28-11-23> 

@doc raw"""
An abstract type for the measurement of the transfer function of an optical system. It can be either a PSF measurement 
[`MeasuredPSF`](@ref) or a [`MeasuredOTF`](@ref).
"""
abstract type MeasuredTransferFunction{N} <: TransferFunction{N} end

@doc """
If the pupil function of the system is symmetric, the OTF as well as the PSF are radially symmteric which can be used to optimize the calculations
"""
@traitdef RadiallySymmetric{TF<:TransferFunction}

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
function Base.hash(tf::TF, h::UInt) where {TF<:TransferFunction}
    hashed = hash(TransferFunction, h)
    hashed = hash(nameof(TF), hashed)

    for f in fieldnames(TF)
        hashed = hash(getfield(tf, f), hashed)
    end

    return hashed
end

