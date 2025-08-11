using TransferFunctions: TransferFunctions as TF
using Combinatorics
using OffsetArrays: OffsetArray as OA
using JET

@testset "trivia $fn" for fn in (TF.corr, TF.conv)
    @test fn(ones(2, 2), ones(2, 2)) == fill(4, (2, 2))
    @test fn(fill(1 + 1im, (2, 2)), ones(2, 2)) == fill(4 + 4im, (2, 2))
    @test fn(fill(1 + 1im, (2, 2)), ones(1, 1)) == fill(1 + 1im, (2, 2))
    @test fn(fill(1 + 1im, (2, 2)), im * ones(1, 1)) == fill(1 - 1im, (2, 2))
    @test fn(ones(2, 2, 2), ones(1, 1, 1)) == ones(2, 2, 2)
    @test_broken fn(ones(2, 2, 2) + im * ones(2, 2, 2), ones(1, 1, 1)) == fill(1 + 1im, (2, 2, 2))
end

@testset "JET: filter $(typeof(A)), $(typeof(B))" for (A, B) in map(Tuple, combinations((ones(3, 3), fill(1 + 1im, (3, 3)), OA(ones(3, 3), -1, -1), OA(fill(1 + 1im, (3, 3)), -1, -1)), 2))
    @test_opt target_modules = (TransferFunctions,) TF.conv(A, B)
    @test_opt target_modules = (TransferFunctions,) TF.corr(A, B)
end
