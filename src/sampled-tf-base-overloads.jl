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
