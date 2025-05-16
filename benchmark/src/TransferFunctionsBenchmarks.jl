# NOTE: Adapted from `JuliaCI/BaseBenchmarks.jl`
module TransferFunctionsBenchmarks

using BenchmarkTools, Dates, Distributed

BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1.0
BenchmarkTools.DEFAULT_PARAMETERS.samples = 10000
BenchmarkTools.DEFAULT_PARAMETERS.time_tolerance = 0.15
BenchmarkTools.DEFAULT_PARAMETERS.memory_tolerance = 0.01

const PARAMS_PATH = normpath(joinpath(dirname(@__FILE__), "..", "params.json"))
const SUITE = BenchmarkGroup()
const MODULES = Dict(
    "array" => :ArrayBenchmarks,
)

load!(id::AbstractString; kwargs...) = load!(SUITE, id; kwargs...)
function load!(group::BenchmarkGroup, id::AbstractString; tune::Bool=true)
    modsym = MODULES[id]
    modpath = joinpath(dirname(@__FILE__), id, "$(modsym).jl")
    Core.eval(TransferFunctionsBenchmarks, :(include($modpath)))
    mod = Core.eval(TransferFunctionsBenchmarks, modsym)
    modsuite = @invokelatest getglobal(mod, :SUITE)
    group[id] = modsuite
    if tune
        results = BenchmarkTools.load(PARAMS_PATH)[1]
        haskey(results, id) && loadparams!(modsuite, results[id], :evals)
    end
    return group
end

loadall!(; kwargs...) = loadall!(SUITE; kwargs...)

function loadall!(group::BenchmarkGroup; verbose::Bool=true, tune::Bool=true)
    for id in keys(MODULES)
        if verbose
            print("loading group $(repr(id))... ")
            time = @elapsed load!(group, id, tune=false)
            println("done (took $time seconds)")
        else
            load!(group, id, tune=false)
        end
    end
    if tune
        results = BenchmarkTools.load(PARAMS_PATH)[1]
        for (id, suite) in group
            haskey(results, id) && loadparams!(suite, results[id], :evals)
        end
    end
    return group
end

function tune!(; verbose=true)
    addprocs(1)
    @info "Loading all benchmarks..."
    TransferFunctionsBenchmarks.loadall!(tune=false)
    @info "Warming up..."
    warmup(SUITE)
    @info "Tuning the parameters..."
    BenchmarkTools.tune!(SUITE; verbose)
    @info "Saving benchmarks parameters to $PARAMS_PATH"
    BenchmarkTools.save(PARAMS_PATH, params(SUITE))
end

const BENCHMARK_PATH = normpath(joinpath(dirname(@__FILE__), "..", "benchmarks.json"))
const NEW_BENCHMARKS = Ref{Union{Missing,BenchmarkGroup}}(missing)
const CURRENT_BENCHMARKS = Ref{Union{Missing,BenchmarkGroup}}(isfile(BENCHMARK_PATH) ? BenchmarkTools.load(BENCHMARK_PATH)[1] : missing)

function compare(heuristic=minimum)
    base = CURRENT_BENCHMARKS[]
    if ismissing(base)
        @warn "No benchmarks are saved!"
        return nothing
    end
    new = NEW_BENCHMARKS[]
    if ismissing(new)
        @warn "No benchmarks where run yet! You must load the needed benchmarks and run them first!"
        return nothing
    end
    return BenchmarkTools.judge(heuristic(new), heuristic(base))
end

function runall!(; verbose=true, tune=false)
    loadall!(; verbose)
    if tune
        tune!()
        loadall!()
    end
    run!(; verbose)
end

function run!(group=SUITE; verbose=true)
    @info "Running loaded benchmarks..."
    new = run(group; verbose)
    NEW_BENCHMARKS[] = new
end

function archive!()
    archive_dir = normpath(joinpath(dirname(@__FILE__), "..", "archive", Dates.format(now(), "yyyymmdd_HHMM")))
    if !isdir(archive_dir)
        @info "Creating archive directory at `$archive_dir`"
        mkpath(archive_dir)
    end
    archive_benchmarks_path = joinpath(archive_dir, "benchmarks.json")
    archive_params_path = joinpath(archive_dir, "params.json")
    if isfile(BENCHMARK_PATH) && isfile(PARAMS_PATH)
        @info "Copying current benchmarks from `$BENCHMARK_PATH` to `$archive_benchmarks_path`, and params from `$PARAMS_PATH` to `$archive_params_path`"
        cp(BENCHMARK_PATH, archive_benchmarks_path)
        cp(PARAMS_PATH, archive_params_path)
    end
    return nothing
end

# FIX: Saved benchmarks should always be up to date with the saved parameters!!! This should be checked when saving <16-05-25> 
function save!()
    new = NEW_BENCHMARKS[]
    if ismissing(new)
        @warn "No benchmarks were run! You must `run!` the benchmarks first..."
        return nothing
    end
    @info "Saving benchmarks to `$BENCHMARK_PATH`"
    BenchmarkTools.save(BENCHMARK_PATH, new)
end

function archive_and_save!()
    @info "Archiving the old benchmarks..."
    archive!()
    @info "Running all benchmarks..."
    runall!(; tune=true)
    @info "Saving the new benchmarks..."
    save!()
end

# TODO: Add some instructions to the info and warn messages... <16-05-25> 

# FIX: This is a bad system! Ideally I would like to save it into the archive automatically at save and only add a link
# to the current benchmarks to the base directory... <16-05-25> 

end # module
