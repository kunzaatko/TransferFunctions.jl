using TransferFunctions
using TransferFunctions: TransferFunctions as TF
using TransferFunctions.Apodization

using OffsetArrays: OffsetArrays as OAs

using FFTW

using TestImages
filenames = ["moonsurface.tiff"]; # NOTE: This is a fix for failing doctests since on download, there is a print-out <19-12-24> 
testimage.(filenames; download_only=false);

return nothing
