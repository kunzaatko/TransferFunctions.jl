using Pkg
Pkg.activate(@__DIR__)
using TransferFunctions
# TODO: Warn on uncommitted changes... The workflow should be to commit, run, and commit --amend <19-12-24> 
using Documenter

# NOTE: When updating, must update also in `docs/make.jl` & 'test/runtests.jl'<18-12-24> 
DocMeta.setdocmeta!(TransferFunctions, :DocTestSetup, :(
        using TransferFunctions;
        using TestImages;
        filenames = ["moonsurface.tiff"]; # NOTE: This is a fix for failing doctests since on download, there is a print-out <19-12-24> 
        testimage.(filenames; download_only=false);
        using FFTW;
        using Logging; # NOTE: This does not need to be in the `make.jl` of docs. We want `@warn ` to function there <19-12-24> 
        Logging.disable_logging(Logging.Warn)
    ); recursive=true)
doctest(TransferFunctions; fix = true)

