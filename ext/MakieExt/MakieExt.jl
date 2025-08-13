module MakieExt
using Makie
using TransferFunctions
using TransferFunctions: Apodization as Apo
using TransferFunctions: SampledMatrix, sampling
import TransferFunctions: scalebar, scalebar!, scalebarformat

include("MakieExt_scalebar.jl")

Makie.convert_arguments(::PointBased, fn::Apo.ApodizationFunction)  = (-1..1,x -> Apo.apodization(fn, x))
Makie.convert_arguments(::ImageLike, A::SampledMatrix)  = (parent(A),)
function Makie.convert_arguments(::Type{<:Scalebar}, A::SampledMatrix) 
    scale = sampling(A)
    @assert allequal(scale) "Sampling must be isotropic otherwise a `scalebar` is not well defined."
    return (first(scale),)
end

end # end module
