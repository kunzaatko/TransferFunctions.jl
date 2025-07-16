using Base: @propagate_inbounds, Indices
using BlockArrays, ImageFiltering
using TensorOperations, LinearAlgebra

include("sampled-array.jl")
include("flattened-array.jl")
include("filtering-arrays.jl")
include("show.jl")
