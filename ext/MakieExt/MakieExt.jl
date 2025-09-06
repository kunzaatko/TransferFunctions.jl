module MakieExt
using Makie
using Makie: Interval
using TransferFunctions
using TransferFunctions: Apodization as Apo
using TransferFunctions: SampledMatrix, sampling
import TransferFunctions: scalebar, scalebar!, scalebarformat
using TransferFunctions: TransferFunctions as TF

include("MakieExt_scalebar.jl")

Makie.convert_arguments(::PointBased, fn::Apo.ApodizationFunction)  = (-1..1,x -> Apo.apodization(fn, x))
function Makie.convert_arguments(::ImageLike, A::SpatialMatrix)
    return ((Interval(lims...) for lims in extrema.(axes(A)))..., parent(A),)
end
Makie.convert_arguments(il::ImageLike, A::PSFArray) = Makie.convert_arguments(il, parent(A))
function Makie.convert_arguments(::Type{<:Scalebar}, A::SampledMatrix) 
    scale = sampling(A)
    @assert allequal(scale) "Sampling must be isotropic otherwise a `scalebar` is not well defined."
    return (first(scale),)
end

end # end module
