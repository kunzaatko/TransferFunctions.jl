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
    @test fn(ones(2, 2, 2) + im * ones(2, 2, 2), ones(1, 1, 1)) == fill(1 + 1im, (2, 2, 2))
end

@testset "corr and filtering" begin
    A = reshape(1:16, (4, 4))
    A_fm = filtering_matrix(A, (-1:1, -1:1), :circular)
    fA_fm = reshape(A_fm' * ones(9), size(A))
    fA_corr = TF.corr(A, OA(ones(3, 3), -1:1, -1:1))
    @test fA_fm == fA_corr
end

if VERSION >= v"1.12-rc"
    @testset "JET: filter $(typeof(A)), $(typeof(B))" for (A, B) in map(Tuple, combinations((ones(3, 3), fill(1 + 1im, (3, 3)), OA(ones(3, 3), -1, -1), OA(fill(1 + 1im, (3, 3)), -1, -1)), 2))
        function fn_filter(@nospecialize f)
            f !== Base.materialize
        end
        @test_opt target_modules = (TransferFunctions,) function_filter = fn_filter TF.conv(A, B)
        @test_opt target_modules = (TransferFunctions,) function_filter = fn_filter TF.corr(A, B)
    end
end

@testset "ImageCore" begin
    using ImageCore, TestImages
    using OffsetArrays: OffsetArrays as OAs

    K_rand = OAs.centered(rand(3, 3))
    K_rand ./= sum(K_rand)

    @testset "Gray image" begin
        img_gray = TestImages.testimage("mandril_gray")

        @test TF.corr(img_gray, K_rand) isa AbstractMatrix{<:Gray}
        @test TF.conv(img_gray, K_rand) isa AbstractMatrix{<:Gray}
    end

    @testset "RGB image" begin
        img_rgb = TestImages.testimage("mandril_color")

        @test TF.corr(img_rgb, K_rand) isa AbstractMatrix{<:RGB}
        @test TF.conv(img_rgb, K_rand) isa AbstractMatrix{<:RGB}
    end
end
