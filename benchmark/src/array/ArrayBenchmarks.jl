module ArrayBenchmarks
using TransferFunctions
using TransferFunctions: FilteringMatrix, CirculantTensor

using ImageFiltering: ImageFiltering as IF

using OffsetArrays
using OffsetArrays: OffsetArrays as OAs
using OffsetArrays: OffsetArray as OA

include(joinpath(dirname(@__FILE__), "..", "utils", "RandUtils.jl"))

using .RandUtils
using BenchmarkTools

const SUITE = BenchmarkGroup()

g = addgroup!(SUITE, "types")

include("create_kernel.jl")

size_str(s) = join(map(string, s), "×")
dim_sizes = Dict( # Realistic application sizes
    1 => Dict(
        :Asize => [(100,), (1000,)],
        :Ksize => [(5,), (15,)]
    ),
    2 => Dict(
        :Asize => [(50, 50), (500, 500)],
        :Ksize => [(5, 5), (15, 15)]
    ),
    3 => Dict(
        :Asize => [(30, 30, 15), (200, 200, 100)],
        :Ksize => [(5, 5, 3), (15, 15, 7)]
    )
)
for d in (1, 2, 3)
    for Asize in dim_sizes[d][:Asize]
        for Ksize in dim_sizes[d][:Ksize]
            for (ptag, p) in [("inner", (x, A, K) -> :($x($A, $K))), ("padded", (x, A, K) -> :($x($A, $K, "replicate")))]
                for (t, tn) in [(circulant, "CirculantTensor"), (FilteringMatrix, "FilteringMatrix")]
                    g[tn, "construction", "$(d)D", "K_$(size_str(Ksize)):A_$(size_str(Asize))", ptag] = @benchmarkable $p($t, A, K) setup = begin
                        A = $samerand($(Asize)...)
                        K = $centered_monotone_kernel($(Ksize)...)
                    end
                end
                # NOTE: Until optimized, I need to reduce the benchmarking only to arrays which are small <13-05-25> 
                if prod(Ksize)^2 * prod(Asize) <= 80e6
                    g["CirculantTensor", "operations", "$(d)D", "contraction", "CT*K", "K_$(size_str(Ksize)):A_$(size_str(Asize))", ptag] = @benchmarkable begin
                        conv(CT, K)
                    end setup = begin
                        A = $samerand($(Asize)...)
                        K = $centered_monotone_kernel($(Ksize)...)
                        CT = eval($p($circulant, A, K))
                    end
                    g["FilteringMatrix", "operations", "$(d)D", "matmul", "FM*K", "K_$(size_str(Ksize)):A_$(size_str(Asize))", ptag] = @benchmarkable begin
                        OAs.no_offset_view(FM) * OAs.no_offset_view(K)[:]
                    end setup = begin
                        A = $samerand($(Asize)...)
                        K = $centered_monotone_kernel($(Ksize)...)
                        FM = eval($p($FilteringMatrix, A, K))
                    end
                end
                if (prod(Ksize) * prod(Asize))^2 <= 1e10
                    g["FilteringMatrix", "operations", "$(d)D", "matmul", "FM'*FM", "K_$(size_str(Ksize)):A_$(size_str(Asize))", ptag] = @benchmarkable begin
                        OAs.no_offset_view(FM)' * OAs.no_offset_view(FM)
                    end setup = begin
                        A = $samerand($(Asize)...)
                        K = $centered_monotone_kernel($(Ksize)...)
                        FM = eval($p($FilteringMatrix, A, K))
                    end
                end
            end
        end
    end
end

end # module
