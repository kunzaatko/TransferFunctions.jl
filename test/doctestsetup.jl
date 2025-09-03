using TransferFunctions
using TransferFunctions: TransferFunctions as TF
using TransferFunctions.Apodization

using OffsetArrays: OffsetArrays as OAs
using BenchmarkTools

using FFTW

ENV["COLUMNS"] = 100
ENV["LINES"] = 80

using TestImages
filenames = ["moonsurface.tiff"]; # NOTE: This is a fix for failing doctests since on download, there is a print-out <19-12-24> 
testimage.(filenames; download_only=false);

setup_makiemaestro!() = @eval begin
    using MakieMaestro
    using MakieMaestro.Recipes
    image = Recipes.image
    image! = Recipes.image!
    mosaic = Recipes.mosaic
    nothing
end

setup_params!() = @eval begin
    λ = 488u"nm"
    NA = 1.4
    n = 1.5
    Δ = 64u"nm"
    Δx = 64u"nm"
    Δy = 32u"nm"
    Δz = 15u"nm"
    f = 1.3
    nothing
end

return nothing
