@test TF.SampledArray(Ones(40, 40), (20u"m^-1", 20u"m^-1")) isa TF.SampledArray
@test SpatialArray(Ones(40, 40), (20u"nm", 20u"nm")) isa SpatialArray
@test SpatialArray(Ones(40, 40), 20u"nm") isa SpatialArray
@test_throws DimensionMismatch SpatialArray(Ones(40, 40), (20u"nm", 20u"nm", 20u"nm"))
@test_throws DimensionMismatch SpatialArray(Ones(40, 40), (20u"nm",))

@testset "circulant arrays" begin
    using Base: OneTo

    ## Constructors ##

    O = ones(100, 100, 100)

    @testset "Flattened" begin
        flatten_parent_4d = [view(O, 30:40, 50:60, :), view(O, 20:30, 40:50, :)]
        ## Throws ##

        # Mismatch in inner `ndims`
        @test_throws MethodError TF.flatten(flatten_parent_4d, inner=(2, 3))

        # inner axes mismatch
        @test_throws DimensionMismatch TF.flatten([view(O, 30:39, 50:60, :), view(O, 20:30, 40:50, :)], inner=(2, 3, 4))

        # Non existent dimension
        @test_throws DimensionMismatch TF.flatten(flatten_parent_4d, inner=(2, 3, 5))


        ## Construction ##
        @test TF.flatten(flatten_parent_4d, inner=(2, 3, 4)) isa AbstractArray{<:Any,4}
        @test TF.flatten(flatten_parent_4d, outer=(3,), inner=(1, 2, 4)) isa AbstractArray{<:Any,4}

        # Other dim specializations
        @test_broken TF.flatten(flatten_parent_4d, outer=3)

        # infer non-default inner, outer from single argument
        @test_broken TF.flatten(flatten_parent_4d, inner=(1, 2, 4)) == TF.flatten(flatten_parent_4d, outer=(3,))

        flatten_parent_3d = [ones(10, 10), zeros(10, 10)]

        ## Methods ##
        F = TF.flatten(flatten_parent_3d)
        @test size(F) == (2, 10, 10)
        @test length(F) == 200
        @test axes(F) == (OneTo(2), OneTo(10), OneTo(10))

        ## Field correctness ##
        @test intersect(Set(F.innermap), Set(F.outermap)) |> isempty
        @test setdiff(Set((F.innermap..., F.outermap...)), Set(1:ndims(F))) |> isempty

        ## Array correctness ##

        @test stack(flatten_parent_3d; dims=1) == F
    end

    K = OA(zeros(21, 21), -10:10, -10:10)
    K[0, 0] = 0.5
    K[-1, -1] = K[-1, 1] = K[1, -1] = K[1, 1] = 0.5 / 4

    A = rand(30, 30)

    @testset "CirculantTensor" begin
        @test circulant(A, (-5:5, -5:5)) isa TF.CirculantTensor{<:Any,2,typeof(A)}
        @test circulant(A, (-5:5, -5:3)) isa TF.CirculantTensor{<:Any,2,typeof(A)} # Non-square kernel indices
        @test circulant(A, K) isa TF.CirculantTensor{<:Any,2,typeof(A)}

        @testset "IF.Padded constructors" begin
            local A = reshape(1:(9*9), 9, 9)
            local Kinds = (-2:2, -2:2)
            borders = [
                ("replicate",
                    [
                        1 1 1 10 19;
                        1 1 1 10 19;
                        1 1 1 10 19;
                        2 2 2 11 20;
                        3 3 3 12 21
                    ], axes(A)
                ),
                ("symmetric",
                    [
                        11 2 2 11 20;
                        10 1 1 10 19;
                        10 1 1 10 19;
                        11 2 2 11 20;
                        12 3 3 12 21
                    ], axes(A)
                ),
                ("circular",
                    [
                        71 80 8 17 26;
                        72 81 9 18 27;
                        64 73 1 10 19;
                        65 74 2 11 20;
                        66 75 3 12 21
                    ], axes(A)
                ),
                ("reflect",
                    [
                        21 12 3 12 21;
                        20 11 2 11 20;
                        19 10 1 10 19;
                        20 11 2 11 20;
                        21 12 3 12 21
                    ], axes(A)
                ),
                (IF.Fill(0.0, (2, 2), (2, 2)),
                    [
                        0 0 0 0 0;
                        0 0 0 0 0;
                        0 0 1 10 19;
                        0 0 2 11 20;
                        0 0 3 12 21
                    ], axes(A)
                ),
                (IF.Fill(0.0),
                    [
                        0 0 0 0 0;
                        0 0 0 0 0;
                        0 0 1 10 19;
                        0 0 2 11 20;
                        0 0 3 12 21
                    ], axes(A)
                )
            ]
            for (bord, out, a) in borders
                ct = circulant(A, Kinds, bord)
                @test ct.interior == a
                @test OAs.no_offset_view(ct[1, 1, :, :]) == out
            end
        end

        # Different eltypes
        @test eltype(circulant(ones(Int, 30, 30), K)) == Int

        C_4D = circulant(A, K)

        # Correct output indices
        @test ndims(circulant(A, K)) == 4
        @test axes(C_4D)[3:4] == axes(K) == C_4D.kern
        @test axes(C_4D)[1:2] == C_4D.interior


        @tensor B[a, b] := OAs.no_offset_view(C_4D)[a, b, c, d] * OAs.no_offset_view(K)[c, d]
        @test B isa AbstractMatrix
        @test size(B) == length.(C_4D.interior)
    end

    @testset "FilteringMatrix" begin
        K_small = OA(K[-4:4, -4:4], -4:4, -4:4)
        C_4D = circulant(A, K_small)

        ## Constructors ##
        @test TF.FilteringMatrix(A, K_small) isa TF.FilteringMatrix
        @test TF.FilteringMatrix(C_4D) isa TF.FilteringMatrix
        @test TF.FilteringMatrix(A, (-4:4, -4:4)) isa TF.FilteringMatrix

        FM_2D = TF.FilteringMatrix(A, K_small)
        FM_2D_small = TF.FilteringMatrix(A[1:10, 1:10], (-1:1, -1:1))
        @test FM_2D.Kaxes == C_4D.kern
        @test FM_2D.Aaxes == C_4D.interior

        # matmul sizes
        @test (FM_2D * K_small[:]) isa AbstractVector
        @test size(FM_2D_small' * FM_2D_small) == (9, 9)
        @test size(FM_2D_small * FM_2D_small') == (8 * 8, 8 * 8)
        @test length(FM_2D_small * K[-1:1, -1:1][:]) == 8 * 8

        FM_2D_pad = TF.FilteringMatrix(A, K_small, "replicate")
        @test (FM_2D_pad * K_small[:]) isa AbstractVector
        @test length((FM_2D_pad * K_small[:])) == length(A)

        FM_1D = TF.FilteringMatrix(1:90, (-1:1,))
        @test axes(FM_1D, 2) == -1:1
        @test_throws ArgumentError (FM_1D * OA(ones(3), -1:1)) # offsets are not supported
        @test_throws ArgumentError (FM_1D * ones(3)) # offsets are not supported
        @test OAs.no_offset_view(FM_1D) * ones(3) isa AbstractVector
        @test length(OAs.no_offset_view(FM_1D) * ones(3)) == 88
    end

    @testset "Dimensions" begin
        function centered_monotone_kernel(s...)
            @assert all(isodd, s)
            K = OAs.centered(Array{Float64}(undef, s))
            max_hypot = hypot(maximum.(map(x -> abs.(x), extrema.(axes(K))))...)
            K .= [cos(hypot(Tuple(i)...) ./ max_hypot) for i in CartesianIndices(K)]
            K .+= rand(s)
            K ./= sum(K)
            return K
        end
        for (Asize, Ksize) in [((10,), (3,)), ((50, 50), (5, 5)), ((12, 12, 4), (3, 3, 3))]
            for p in [(x, A, K) -> x(A, K), (x, A, K) -> x(A, K, "replicate")]
                for (t, c) in [(TF.CirculantTensor, circulant), (TF.FilteringMatrix, TF.FilteringMatrix)]
                    A = rand(Asize...)
                    K = centered_monotone_kernel(Ksize...)
                    @test p(c, A, K) isa t
                end
                CT = p(circulant, A, K)
                CT_conv_K = TF.conv(CT, K)
                @test ndims(CT_conv_K) == ndims(A)
                FM = p(TF.FilteringMatrix, A, K)
                if length(Asize) == 1 # offset of 1D filtering matrix makes it incompatible with matrix multiplication
                    FM = OAs.no_offset_view(FM)
                end
                K = OAs.no_offset_view(K)
                @test FM * K[:] isa AbstractVector
                @test FM' * FM isa AbstractMatrix
            end
        end
    end

    ## Filtering and Correctness ##

    A_1D = Vector(1:4)
    CT_1D_1D = circulant(A_1D, (-1:1,))
    @test OAs.no_offset_view(CT_1D_1D) == [1 2 3; 2 3 4]
    @test axes(CT_1D_1D, 1) == 2:3

    A_2D = reshape(1:12, 4, 3)
    CT_2D_2D = circulant(A_2D, (0:1, 0:1))
    @test OAs.no_offset_view(CT_2D_2D[1, 1, :, :]) == [1 5; 2 6]
    @test OAs.no_offset_view(CT_2D_2D[1, 2, :, :]) == [5 9; 6 10]
    @test CT_2D_2D[2, 2, 1, 1] == 11 # NOTE: There is an offsetted kernel with indices 0:1×0:1 <05-05-25> 
    @test axes(CT_2D_2D) == (1:3, 1:2, 0:1, 0:1)

    K = OA([1 0; 0 0], 0:1, 0:1)
    FM_2D_2D = TF.FilteringMatrix(A_2D, K)

    @test OAs.no_offset_view(reshape(FM_2D_2D * K[:], FM_2D_2D.Aaxes)) == A_2D[FM_2D_2D.Aaxes...]

    K_rand = OAs.centered(rand(3, 3))
    K_rand ./= sum(K_rand)

    img = float.(gray.(TestImages.testimage("mandril_gray")))

    fm_img = TF.FilteringMatrix(img, K_rand, "replicate")
    fm_filt = reshape(fm_img * K_rand[:], fm_img.Aaxes)

    fft_filt = imfilter(img, K_rand)

    @test fm_filt == fft_filt
end
