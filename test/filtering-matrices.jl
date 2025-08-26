using TransferFunctions
using TransferFunctions: TransferFunctions as TF
using OffsetArrays: OffsetArray as OA
using JET

A = reshape(1:25, (5, 5))
@testset "Construction" begin
    @test TF.FilteringMatrix(A, (-1:1, -1:1)) isa TF.FilteringMatrix
    @test TF.FilteringMatrix(A, OA(ones(3, 3), -1:1, -1:1)) == TF.FilteringMatrix(A, (-1:1, -1:1))
    @test TF.FilteringMatrix(A, (-1:1, -1:1)) == TF.FilteringMatrix(A, (-1:1, -2:0))
    @test TF.FilteringMatrix(A, (-1:1, -1:1)) == TF.FilteringMatrix(A, (1:3, 1:3))
    @test TF.FilteringMatrix(OA(A, -1, -1), (-1:1, -1:1)) == TF.FilteringMatrix(A, (-1:1, -2:0))
end

@testset "AbstractArray interface" begin
    let F = TF.FilteringMatrix(A, (-1:1, -1:1))
        @test axes(F) == (1:9, 1:9)
        @test size(F) == (9, 9)
        @test F.interior == (3, 3)
        @test F[begin] == 1
        @test F[end] == 25
        @test_throws BoundsError F[begin-1]
        @test_throws BoundsError F[end+1]
        @test F[:, 1] == F.parent[1:3, 1:3][:]
    end
    let F = TF.FilteringMatrix(A, (-1:1, -2:2))
        @test axes(F) == (1:15, 1:3)
        @test size(F) == (15, 3)
        @test F.interior == (3, 1)
        @test F[begin] == 1
        @test F[end] == 25
        @test_throws BoundsError F[begin-1]
        @test_throws BoundsError F[end+1]
        @test F[:, 1] == F.parent[1:3, 1:5][:]
    end
end

@testset "filtering" begin
    let F = TF.FilteringMatrix(ones(9, 9), (-1:1, -1:1))
        @test all(==(1), (F' * ones(9) ./ 9))
    end
    @test size(filtering_matrix(A, (-1:1, -1:1)))[1] == 9
    @test size(filtering_matrix(A, (-1:1, -1:1), :circular))[1] == 9

    K_rand = OAs.centered(rand(3, 3))
    K_rand ./= sum(K_rand)
    @testset "corr" begin
        @test let A = ones(10, 10)
            fm_A = TF.filtering_matrix(A, K_rand, :circular)
            fm_corr = reshape(fm_A' * K_rand[:], axes(A))
            fft_corr = TF.corr(A, K_rand)
            isapprox(fm_corr, fft_corr, rtol=1e-3)
        end
        @test let A = rand(10, 10)
            fm_A = TF.filtering_matrix(A, K_rand, :circular)
            fm_corr = reshape(fm_A' * K_rand[:], axes(A))
            fft_corr = TF.corr(A, K_rand)
            isapprox(fm_corr, fft_corr, rtol=1e-3)
        end
    end
    @testset "conv" begin
        @test let A = rand(10, 10)
            fm_A = TF.filtering_matrix(reflect(A), K_rand, :circular)
            fm_conv = reshape(reflect(fm_A' * K_rand[:]), axes(A))
            fft_conv = TF.conv(A, K_rand)
            isapprox(fm_conv, fft_conv, rtol=1e-3)
        end
    end
end

@testset "ImageCore" begin
    using ImageCore, TestImages
    using OffsetArrays: OffsetArrays as OAs
    using ImageFiltering

    @testset "Array types compatibility" begin
        K_rand = OAs.centered(rand(3, 3))
        K_rand ./= sum(K_rand)

        @testset "Gray image" begin # Gray image
            img_gray = TestImages.testimage("mandril_gray")[1:10, 1:10]

            fm_img = TF.filtering_matrix(img_gray, K_rand, :circular)
            fm_corr = reshape(fm_img' * K_rand[:], axes(img_gray))
            fft_corr = TF.corr(img_gray, K_rand)
            isapprox(fm_corr, fft_corr, rtol=1e-3)

        end

        @test begin # RGB image
            img_rgb = TestImages.testimage("mandril_color")[1:10, 1:10]

            fm_img = TF.filtering_matrix(img_rgb, K_rand, :circular)
            fm_corr = reshape(fm_img' * K_rand[:], axes(img_rgb))
            fft_corr = TF.corr(img_rgb, K_rand)
            isapprox(fm_corr, fft_corr, rtol=1e-3)
        end
    end
end

if VERSION <= v"1.12"
    @testset "JET: `getindex`" begin
        let F = TF.FilteringMatrix(ones(9, 9), (-1:1, -1:1))
            @test_opt target_modules = (TransferFunctions,) getindex(F, 1, 1)
            @test_opt target_modules = (TransferFunctions,) getindex(F, 1, :)
        end

        let F = TF.FilteringMatrix(ones(9, 9, 9), (-1:1, -1:1, -1:1))
            @test_opt target_modules = (TransferFunctions,) getindex(F, 1, 1)
            @test_opt target_modules = (TransferFunctions,) getindex(F, 1, :)
        end
    end
end
