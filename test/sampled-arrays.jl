using TransferFunctions
using TransferFunctions: TransferFunctions as TF

@test TF.SampledArray(ones(10, 10), (20u"m^-1", 20u"m^-1")) isa TF.SampledArray
@test let s = SpatialArray(ones(), ())
    s isa SampledArray && ndims(s) == 0
end
@test let sa = TF.SampledArray(ones(10, 10), 20u"m^-1")
    sa isa TF.SampledArray && TF.sampling(sa) == (20u"m^-1", 20u"m^-1")
end
@test SpatialArray(ones(10, 10), (20u"nm", 20u"nm")) isa SpatialArray
@test SpatialArray(ones(10, 10), (20.0u"nm", 20u"nm")) isa SampledArray{<:Any,typeof(20.0u"nm")} # Mixed types in sampling are promoted
@test SpatialArray(ones(10, 10), (20u"nm", 2e-3u"m")) isa SpatialArray # Mixed units are promoted
@test let s = SpatialArray(ones(10, 10), 20u"nm")
    s isa SpatialArray && allequal(TransferFunctions.sampling(s))  # Single sampling is inferred for all dimensions
end
@test TransferFunctions.sampling(similar(SpatialArray(ones(10, 10), (50u"nm", 50u"nm")))) == (50u"nm", 50u"nm") # `similar` preserves sampling
@test_throws MethodError SpatialArray(ones(10, 10), (20u"nm", 20u"nm", 20u"nm"))
@test_throws MethodError SpatialArray(ones(10, 10), (20u"nm",))

@test SpatialMatrix(ones(10, 10), 20u"nm") isa SpatialArray
@test SpatialVector(ones(10), 20u"nm") isa SpatialArray
