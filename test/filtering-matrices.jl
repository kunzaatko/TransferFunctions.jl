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
end

@testset "JET: `getindex`" begin
    let F = TF.FilteringMatrix(ones(9, 9), (-1:1, -1:1))
        @test_opt getindex(F, 1, 1)
        @test_opt getindex(F, 1, :)
    end

    let F = TF.FilteringMatrix(ones(9, 9, 9), (-1:1, -1:1, -1:1))
        @test_opt getindex(F, 1, 1)
        @test_opt getindex(F, 1, :)
    end
end

# using OffsetArrays: OffsetArray as OA
# using OffsetArrays: OffsetArrays as OAs

# K_small = OA(K[-4:4, -4:4], -4:4, -4:4)
# C_4D = circulant(A, K_small)
#
# ## Constructors ##
# @test TF.FilteringMatrix(A, K_small) isa TF.FilteringMatrix
# @test TF.FilteringMatrix(C_4D) isa TF.FilteringMatrix
# @test TF.FilteringMatrix(A, (-4:4, -4:4)) isa TF.FilteringMatrix
#
# FM_2D = TF.FilteringMatrix(A, K_small)
# FM_2D_small = TF.FilteringMatrix(A[1:10, 1:10], (-1:1, -1:1))
# @test FM_2D.Kaxes == C_4D.kern
# @test FM_2D.Aaxes == C_4D.interior
#
# # matmul sizes
# @test (FM_2D * K_small[:]) isa AbstractVector
# @test size(FM_2D_small' * FM_2D_small) == (9, 9)
# @test size(FM_2D_small * FM_2D_small') == (8 * 8, 8 * 8)
# @test length(FM_2D_small * K[-1:1, -1:1][:]) == 8 * 8
#
# FM_2D_pad = TF.FilteringMatrix(A, K_small, "replicate")
# @test (FM_2D_pad * K_small[:]) isa AbstractVector
# @test length((FM_2D_pad * K_small[:])) == length(A)
#
# FM_1D = TF.FilteringMatrix(1:10, (-1:1,))
# @test axes(FM_1D, 2) == -1:1
# @test_throws ArgumentError (FM_1D * OA(ones(3), -1:1)) # offsets are not supported
# @test_throws ArgumentError (FM_1D * ones(3)) # offsets are not supported
# @test OAs.no_offset_view(FM_1D) * ones(3) isa AbstractVector
# @test length(OAs.no_offset_view(FM_1D) * ones(3)) == 8
#
# function centered_monotone_kernel(s...)
#     @assert all(isodd, s)
#     K = OAs.centered(Array{Float64}(undef, s))
#     max_hypot = hypot(maximum.(map(x -> abs.(x), extrema.(axes(K))))...)
#     K .= [cos(hypot(Tuple(i)...) ./ max_hypot) for i in CartesianIndices(K)]
#     K .+= rand(s)
#     K ./= sum(K)
#     return K
# end
# for (Asize, Ksize) in [((10,), (3,)), ((50, 50), (5, 5)), ((12, 12, 4), (3, 3, 3))]
#     for p in [(x, A, K) -> x(A, K), (x, A, K) -> x(A, K, "replicate")]
#         for (t, c) in [(TF.CirculantTensor, circulant), (TF.FilteringMatrix, TF.FilteringMatrix)]
#             A = rand(Asize...)
#             K = centered_monotone_kernel(Ksize...)
#             @test p(c, A, K) isa t
#         end
#         CT = p(circulant, A, K)
#         CT_conv_K = TF.conv(CT, K)
#         @test ndims(CT_conv_K) == ndims(A)
#         FM = p(TF.FilteringMatrix, A, K)
#         if length(Asize) == 1 # offset of 1D filtering matrix makes it incompatible with matrix multiplication
#             FM = OAs.no_offset_view(FM)
#         end
#         K = OAs.no_offset_view(K)
#         @test FM * K[:] isa AbstractVector
#         @test FM' * FM isa AbstractMatrix
#     end
# end
#
# ## Filtering and Correctness ##
#
# A_1D = Vector(1:4)
# CT_1D_1D = circulant(A_1D, (-1:1,))
# @test OAs.no_offset_view(CT_1D_1D) == [1 2 3; 2 3 4]
# @test axes(CT_1D_1D, 1) == 2:3
#
# A_2D = reshape(1:12, 4, 3)
# CT_2D_2D = circulant(A_2D, (0:1, 0:1))
# @test OAs.no_offset_view(CT_2D_2D[1, 1, :, :]) == [1 5; 2 6]
# @test OAs.no_offset_view(CT_2D_2D[1, 2, :, :]) == [5 9; 6 10]
# @test CT_2D_2D[2, 2, 1, 1] == 11 # NOTE: There is an offsetted kernel with indices 0:1×0:1 <05-05-25> 
# @test axes(CT_2D_2D) == (1:3, 1:2, 0:1, 0:1)
#
# K = OA([1 0; 0 0], 0:1, 0:1)
# FM_2D_2D = TF.FilteringMatrix(A_2D, K)
#
# @test OAs.no_offset_view(reshape(FM_2D_2D * K[:], FM_2D_2D.Aaxes)) == A_2D[FM_2D_2D.Aaxes...]
#
# @testset "ImageCore" begin
#     using ImageCore
#
#     @testset "Array types compatibility" begin
#         K_rand = OAs.centered(rand(3, 3))
#         K_rand ./= sum(K_rand)
#
#         @test begin # Gray image
#             img_gray = TestImages.testimage("mandril_gray")[1:10, 1:10]
#
#             fm_img = TF.FilteringMatrix(img_gray, K_rand, "replicate")
#             fm_filt = reshape(fm_img * K_rand[:], fm_img.Aaxes)
#
#             fft_filt = imfilter(img_gray, K_rand)
#
#             fm_filt == fft_filt
#         end
#
#         @test begin # RGB image
#             img_rgb = TestImages.testimage("mandril_color")[1:10, 1:10]
#
#             fm_img = TF.FilteringMatrix(img_rgb, K_rand, "replicate")
#             fm_filt = reshape(fm_img * K_rand[:], fm_img.Aaxes)
#
#             fft_filt = imfilter(img_rgb, K_rand)
#
#             fm_filt == fft_filt
#         end
#     end
# end
