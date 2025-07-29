module MakieExt
using TransferFunctions
using TransferFunctions: Apodization as Apo
using Makie

Makie.convert_arguments(::PointBased, fn::Apo.ApodizationFunction)  = (-1..1,x -> Apo.apodization(fn, x))

end # end module
