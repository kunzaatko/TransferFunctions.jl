using TransferFunctions
using Documenter

DocMeta.setdocmeta!(TransferFunctions, :DocTestSetup, :(
        include(joinpath(@__DIR__, "doctestsetup.jl"));
        using Logging;
        # NOTE: Not necessary in `docs/make.jl`. `@warn` should work there <19-12-24>
        # FIX: I cannot get the doctest filtering to work  <17-09-25> 
        Logging.disable_logging(Logging.Warn)
    ); recursive=true)
doctest(TransferFunctions; fix=true)
