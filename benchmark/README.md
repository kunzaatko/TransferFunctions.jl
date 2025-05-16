# TransferFunctions.jl Benchmarks

This directory contains benchmarking code for the TransferFunctions.jl package. It is built on top of the
[BenchmarkTools.jl](https://github.com/JuliaCI/BenchmarkTools.jl) package and is adapted from
[BaseBenchmarks.jl](https://github.com/JuliaCI/BaseBenchmarks.jl).

## Setup

1. Navigate to this directory and instantiate the benchmark environment:

```bash
cd /path/to/TransferFunctions.jl/benchmark
julia --project -e 'using Pkg; Pkg.instantiate()'
```

## Basic usage of the `TransferFunctionsBenchmarks` Module

### Running Benchmarks

Start Julia with the benchmark project:

```bash
julia --project
```

Then run benchmarks using the provided functions:

```julia
using TransferFunctionsBenchmarks

# Tune the benchmark parameters (recommended before first run)
TransferFunctionsBenchmarks.tune!()
```

Then you can run all the benchmarks with
```julia
# Load and run all benchmarks
B = TransferFunctionsBenchmarks.runall!()
```
or run only selected benchmarks
```julia
TransferFunctionsBenchmarks.load!("array")
TransferFunctionsBenchmarks.run!("array")
```
or
```julia
TransferFunctionsBenchmarks.loadall!()
TransferFunctionsBenchmarks.run!(TransferFunctionsBenchmarks.SUITE[@tagged "matmul"])
```

If you already have benchmarks saved from previous runs, you can analyse the results and compare them with the current.

### Archiving and Managing Benchmark Results

When you are satisfied with the improvements, you should archive the current baseline benchmarks for future reference
and save the current ones.

Benchmark results are saved in the following locations:
- Current benchmark results: `benchmark/benchmarks.json`
- Benchmark parameters: `benchmark/params.json`
- Archived benchmarks: `benchmark/archive/YYYYMMDD_HHMM/benchmarks.json`
- Archived parameters: `benchmark/archive/YYYYMMDD_HHMM/params.json`

To archive the old results and create a new baseline, run

```julia
# Save benchmark results
TransferFunctionsBenchmarks.archive_and_save!()
```

This function:
1. Archives the current benchmarks and parameters to a timestamped directory
2. Tunes the benchmark parameters and saves them
3. Loads and runs all benchmarks
4. Saves the new benchmark results

### Loading Specific Benchmark Groups

You can also load and run specific benchmark groups:

```julia
# Load just the array benchmarks
TransferFunctionsBenchmarks.load!("array")

# Run only loaded benchmarks
TransferFunctionsBenchmarks.run!()
```

### Running with Tag Filters

You can use BenchmarkTools tag filtering to run specific benchmarks:

```julia
using BenchmarkTools
using TransferFunctionsBenchmarks

# Load all benchmarks
TransferFunctionsBenchmarks.loadall!()

# Run only Array construction benchmarks
TransferFunctionsBenchmarks.run!(TransferFunctionsBenchmarks.SUITE["array"][@tagged "construction"])

# Run only 2D benchmarks
TransferFunctionsBenchmarks.run!(TransferFunctionsBenchmarks.SUITE[@tagged "2D"])
```

## Comparing Benchmark Results

After running benchmarks, you can compare them to previously saved results:

```julia
# Compare new benchmarks to previously saved ones
judgement = TransferFunctionsBenchmarks.compare()

regs = regressions(judgement)
pairs = leaves(regs) # an array of (ID, `TrialJudgement`) pairs

imps = improvements(judgement)
```

The `compare` function accepts a heuristic function that's applied to benchmark results (default is `minimum`). You can also use `median` or other BenchmarkTools statistics:

```julia
# Compare using median instead of minimum
judgement = TransferFunctionsBenchmarks.compare(median)
```

To see where the bottlenecks are you can profile code in the benchmark
```julia
using ProfileCanvas # or ProfileView, PProf etc.
@profile run(BaseBenchmarks.SUITE[["FilteringMatrix", "operations", "2D", "matmul", "FM'*FM", "K_5×5:A_50×50", "inner"]])
```

For further comparison and analysis options [BenchmarkTools.jl
documentation](https://github.com/JuliaCI/BenchmarkTools.jl).


## Adding New Benchmarks

To add new benchmark categories:

1. Create a new module in `benchmark/src/your_category/YourCategoryBenchmarks.jl`
2. Define a `SUITE` BenchmarkGroup in that module
3. Add your module to the `MODULES` dictionary in `benchmark/src/TransferFunctionsBenchmarks.jl`

Your module should follow this structure:

```julia
module YourCategoryBenchmarks
using BenchmarkTools
using TransferFunctions

const SUITE = BenchmarkGroup()

# Add your benchmarks here
g = addgroup!(SUITE, "your_group")
g["benchmark_name"] = @benchmarkable some_function(x) setup = begin
    x = rand(100)
end

end # module
```

## Contributing Benchmarks

When adding benchmarks, please follow these guidelines:

- Make sure your benchmarks are deterministic and not affected by global state
- Use appropriate setup blocks to isolate benchmark preparation
- Add benchmarks for all performance-critical code paths
- Use meaningful tags and groups to organize your benchmarks
- Test benchmarks with a variety of input sizes representative of real use
