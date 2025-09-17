using TransferFunctions
using Aqua, Test, Documenter, CompatHelperLocal

const run_all = isempty(ARGS) ? true : false

skip = Dict{String,Bool}(
    "compat" => !(VERSION >= v"1.9"), # NOTE: `CompatHelperLocal` only compatible with later Julia version <28-02-25> 
    "aqua" => !haskey(ENV, "GITHUB_ACTIONS") && !haskey(ENV, "RUNTESTS_FULL"),
    "doctests" => !haskey(ENV, "RUNTESTS_FULL") && !(haskey(ENV, "RUNNER_OS") && ENV["RUNNER_OS"] == "Linux"),
    "ambiguities" => true # FIX: Fix the ambiguities <24-04-25> 
)

function should_test(arg::String)::Bool
    global run_all
    if run_all
        return !get(skip, arg, false)
    elseif arg in ARGS
        return true
    end
    return false
end

macro cond_testset(name, block)
    quote
        if should_test($name)
            @testset $name begin
                esc($block)
            end
        end
    end
end

@testset "TransferFunctions.jl" begin
    @testset "Code quality" begin
        @cond_testset "aqua" begin
            Aqua.test_all(
                TransferFunctions;
                ambiguities=false,
            )
        end

        @cond_testset "ambiguities" begin
            aqua_ambiguities = false
            if aqua_ambiguities
                Agua.test_ambiguities(TransferFunctions)
            else
                @test length(Test.detect_ambiguities(TransferFunctions)) == 0
            end
        end

        @cond_testset "compat" begin
            @test CompatHelperLocal.check(TransferFunctions; checktest=false)
        end
    end

    @cond_testset "doctests" begin
        # FIX: When running locally, do not ask for SSH key password <10-12-23> 
        # NOTE: Show for `Unitful.jl` does nm⁻¹ on macOS and nm^-1 on Linux. This is necessary, since the `jldoctest` is only one
        # NOTE: Better than doc-testing in `make.jl` because, I can track the coverage and it doesn't take time when
        # building documentation
        # NOTE: When updating, must update also in `docs/make.jl` & `test/fix_doctests.jl` <18-12-24> 
        DocMeta.setdocmeta!(TransferFunctions, :DocTestSetup, :(
                include(joinpath(@__DIR__, "doctestsetup.jl"));
                using Logging;
                # NOTE: Not necessary in `docs/make.jl`. `@warn` should work there <19-12-24>
                # FIX: I cannot get the doctest filtering to work  <17-09-25> 
                Logging.disable_logging(Logging.Warn)
            ); recursive=true)
        !haskey(ENV, "FIX_DOCTESTS") && @info "You can fix doctests by setting `ENV[\"FIX_DOCTESTS\"] = true`."
        doctest(TransferFunctions; fix=ifelse(haskey(ENV, "FIX_DOCTESTS"), true, false))
    end

    @cond_testset "utils" begin
        include("utils.jl")
    end

    @cond_testset "types" begin
        include("types.jl")
    end

    @cond_testset "sampled-arrays" begin
        include("sampled-arrays.jl")
    end

    @cond_testset "circulant-tensors" begin
        include("circulant-tensors.jl")
    end

    @cond_testset "filtering-matrices" begin
        include("filtering-matrices.jl")
    end

    @cond_testset "border-arrays" begin
        include("border-arrays.jl")
    end

    @cond_testset "tapered-arrays" begin
        include("tapered-arrays.jl")
    end

    @cond_testset "reflected-arrays" begin
        include("reflected-arrays.jl")
    end

    @cond_testset "fft" begin
        include("fft.jl")
    end

    @cond_testset "filter" begin
        include("filter.jl")
    end

    @cond_testset "interfaces" begin
        include("interfaces.jl")
    end

    @cond_testset "apodization" begin
        include("apodization.jl")
    end

    @cond_testset "optical-transfer-function" begin
        include("optical-transfer-function.jl")
    end

    @cond_testset "point-spread-function" begin
        include("point-spread-function.jl")
    end

    @cond_testset "point-spread-function-models" begin
        include("point-spread-function-models.jl")
    end

    @cond_testset "synthetic-data" begin
        include("synthetic-data.jl")
    end

    @cond_testset "Estimation" begin
        include("estimation.jl")
    end

    @cond_testset "base-overloads" begin
        include("base-overloads.jl")
    end
end

if run_all
    @warn "Skipped: $(keys(skip))"
else
    @info "Ran: $ARGS"
end
