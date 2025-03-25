# # TODO: Test this... <11-03-25> 
# TODO: Should this be a thing for the `MeasuredTransferFunction` as well? <26-08-24> 
# NOTE: Allows broadcasting `func.(tf::ModelTransferFunction, a:b)`where `func` can be any of `psf`,`otf`,`mtf`, etc.
Broadcast.broadcastable(tf::ModelTransferFunction) = Ref(tf)

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

# TODO: Util function for printing the resolution if the both the dimensional resolutions are the same i.e. 64 nm
# instead (64 nm, 64 nm). This will be used for both the sampled types. <19-12-24> 
# TODO: Better Base.show using `typeof` and `nameof` <30-11-23> 
function Base.show(io::IO, ::MIME"text/plain", tf::S) where {S<:SampledTransferFunction}
    print(io, nameof(typeof(tf)), "(")
    show(io, MIME("text/plain"), tf.transfer)
    print(io, ") with Δxy=", allequal(tf.Δxy) ? tf.Δxy[1] : tf.Δxy)
    if !all(iszero.(tf.center))
        print(io, ", δ = ", tf.center)
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
