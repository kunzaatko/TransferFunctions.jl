module MakieExt
using Makie
using TransferFunctions
using TransferFunctions: Apodization as Apo
using TransferFunctions: SpatialMatrix
import TransferFunctions: scalebar, scalebar!, scalebarformat

Makie.convert_arguments(::PointBased, fn::Apo.ApodizationFunction)  = (-1..1,x -> Apo.apodization(fn, x))
Makie.convert_arguments(::ImageLike, A::SpatialMatrix)  = (parent(A),)

include("MakieExt_scalebar.jl")

end # end module
