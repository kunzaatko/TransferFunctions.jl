using TransferFunctions
using TransferFunctions: TransferFunctions as TF
using TransferFunctions: Frequency
using IntervalSets, FourierTools, FFTViews, Distributions, FillArrays, TensorOperations, OffsetArrays, ImageFiltering, ImageCore, TestImages, DataStructures
using OffsetArrays: OffsetArray as OA
using OffsetArrays: OffsetArrays as OAs
using ImageFiltering: ImageFiltering as IF
using Aqua, Test, Documenter, CompatHelperLocal

@testset "TransferFunctions.jl" begin
    @testset "Code quality" begin
        ambiguities = false # FIX: Fix the ambiguities <24-04-25> 
        aqua_ambiguities = false
        @testset "Aqua.jl" begin
            if haskey(ENV, "RUNTESTS_FULL") || haskey(ENV, "GITHUB_ACTIONS")
                Aqua.test_all(
                    TransferFunctions;
                    ambiguities=aqua_ambiguities && ambiguities,
                    # ambiguities=VERSION >= v"1.1" ? (; broken=true) : false
                )
            else
                @warn "Skipping Aqua.jl quality tests. For a full run set `ENV[\"RUNTESTS_FULL\"]=true`."
            end
        end
        @testset "Ambiguities" begin
            if !aqua_ambiguities && ambiguities
                @test length(Test.detect_ambiguities(TransferFunctions)) == 0
            end
        end
        if VERSION >= v"1.9" # NOTE: Only works for later Julia due to new version changes in the package <28-02-25> 
            @testset "Compat" begin
                CompatHelperLocal.@check(checktest = false)
            end
        end
    end

    @testset "DocTests" begin
        # FIX: When running locally, do not ask for SSH key password <10-12-23> 
        # NOTE: Show for `Unitful.jl` does nm⁻¹ on macOS and nm^-1 on Linux. This is necessary, since the `jldoctest` is only one
        if haskey(ENV, "RUNTESTS_FULL") || haskey(ENV, "RUNNER_OS") && ENV["RUNNER_OS"] == "Linux"
            # NOTE: Better than doc-testing in `make.jl` because, I can track the coverage
            # NOTE: When updating, must update also in `docs/make.jl` & `test/fix_doctests.jl` <18-12-24> 
            DocMeta.setdocmeta!(TransferFunctions, :DocTestSetup, :(
                    using TransferFunctions;
                    using TestImages;
                    filenames = ["moonsurface.tiff"]; # NOTE: This is a fix for failing doctests since on download, there is a print-out <19-12-24> 
                    testimage.(filenames; download_only=false);
                    using TransferFunctions: FFTW;
                    using Logging; # NOTE: This does not need to be in the `make.jl` of docs. We want `@warn ` to function there <19-12-24> 
                    Logging.disable_logging(Logging.Warn)
                ); recursive=true)
            !haskey(ENV, "FIX_DOCTESTS") && @info "You can fix doctests by setting `ENV[\"FIX_DOCTESTS\"] = true`."
            doctest(TransferFunctions; fix=ifelse(haskey(ENV, "FIX_DOCTESTS"), true, false))
        else
            @warn "Skipping DocTests. For a full run set `ENV[\"RUNTESTS_FULL\"]=true`."
        end
    end

    @testset "utils.jl" begin
        using TransferFunctions: fillsize, roundupcenter, exactcenter, fftfreqs, posgrid, contained, interior, rounddowncenter, roundcenter, aroundorigin
        using TransferFunctions: PixelSize, Coordinate, Frequency, Length, OriginAt
        using Base: CartesianIndex as CI

        @testset "types" begin
            @testset "units" begin
                @test 1 / 32u"nm" isa Frequency
                Δkx = 1 / 61u"nm"
                @test 1 / Δkx isa Length
            end

            @test (31.5u"nm", 40u"nm", 50u"nm") isa PixelSize{3}
            @test (31.5u"nm", 40u"nm") isa PixelSize{2}

            @test (3, 3, 3) isa Coordinate{3}
            @test (3, 3, 3) isa Coordinate{3,Int}
            @test (3.0, 3.0, 3.0) isa Coordinate{3,<:AbstractFloat}
            @test (2.4, 3, 3) isa Coordinate{3,Real} # NOTE: Must accept diverse types if supplied <26-08-24> 
        end


        @testset "utility functions" begin
            @test fillsize(31u"nm", 2) == (31u"nm", 31u"nm")
            @test_throws MethodError fillsize(31u"nm", Val(2)) == (31u"nm", 31u"nm") # NOTE: method with integer should always be used <05-05-25> 

            @test (posgrid((11, 11), (31u"nm", 31u"nm")) |> first |> first) isa Length
            @test posgrid((11, 11), 31u"nm") == posgrid((11, 11), (31u"nm", 31u"nm"))
            @test_throws MethodError posgrid((11, 11, 11), (31u"nm", 31u"nm"))

            @test (fftfreqs((11, 11), (31u"nm", 31u"nm")) |> first |> first) isa Frequency
            @test fftfreqs((11, 11), 31u"nm") == fftfreqs((11, 11), (31u"nm", 31u"nm"))
            @test_throws MethodError fftfreqs((11, 11, 11), (31u"nm", 31u"nm"))

            @test roundcenter(RoundFromZero, Ones(3, 4, 2)) == CI(2, 3, 2)

            @test roundupcenter(Ones(3, 4, 2)) == CI(2, 3, 2)
            @test roundupcenter(OA(Ones(3, 3, 3), -2, -2, -2)) == CI(0, 0, 0)
            @test roundupcenter(OA(Ones(4, 3, 3), -2, -2, -2)) == CI(1, 0, 0)

            @test rounddowncenter(Ones(3, 4, 2)) == CI(2, 2, 1)
            @test rounddowncenter(OA(Ones(3, 3, 3), -2, -2, -2)) == CI(0, 0, 0)
            @test rounddowncenter(OA(Ones(4, 3, 3), -2, -2, -2)) == CI(0, 0, 0)

            @test exactcenter(Ones(3, 4, 2)) == (2.0, 2.5, 1.5)
            @test exactcenter(OA(Ones(3, 3, 3), -2, -2, -2)) == (0.0, 0.0, 0.0)
            @test exactcenter(OA(Ones(4, 3, 3), -2, -2, -2)) == (0.5, 0.0, 0.0)

            @test contained(Ones(3, 4, 2), (2, 3, 1))
            @test !contained(Ones(3, 4, 2), (8, 3, 1))
            @test contained(OA(Ones(3, 3, 3), -2, -2, -2), (-1, 1, 0))
            @test !contained(OA(Ones(3, 3, 3), -2, -2, -2), (-2, 1, 0))

            # TODO: Test other types of axes that may occur in an array that I use (OffsetAxes) <05-05-25> 
            @test interior(1:9, -1:5) == 2:4
            @test interior(Base.OneTo(9), -1:5) == 2:4
            @test interior((0:3, -1:3, -3:1), (-1:1, -1:1, -1:1)) == (1:2, 0:2, -2:0)
            @test interior((Base.OneTo(3), -1:3, -3:1), (0:1, -1:1, -1:1)) == (1:2, 0:2, -2:0)

            @test aroundorigin((-3:3, -3:5), (2, 1)) == (-1:5, -2:6)
            @test aroundorigin(-3:4, 4) == 1:8
            @test aroundorigin(-3:4) == -3:4
            @test aroundorigin((3, 3, 3)) == (-1:1, -1:1, -1:1)

            @test OriginAt(CI(2, 2, 2))(Ones(3, 3, 3)) == OA(Ones(3, 3, 3), -2, -2, -2)
        end

        @testset "macros" begin
            @test_throws ["atleast one argument"] @macroexpand TF.@require_interface(function some() end) # No arguments
            @test_throws ["atleast one argument"] @macroexpand TF.@require_interface(some())              # No arguments
            @test_throws ["function or a `:call`"] @macroexpand TF.@require_interface(some = 5)           # Not a call
            @test_throws ["known type"] @macroexpand TF.@require_interface(some(a))                       # No type
            @test_throws ["must be abstract"] @macroexpand TF.@require_interface(some(a::Float64))        # Concrete type

            @test (@macroexpand TF.@require_interface(some(a::AbstractFloat))).args[2].head == Symbol("function")
            @test (@macroexpand TF.@require_interface(some(a::AbstractFloat, b::Int))).args[2].head == Symbol("function")
            @test_broken (@macroexpand TF.@require_interface(some(a::AbstractFloat, b::Int)::Float64)).args[2].head == Symbol("function")
            @test_broken (@macroexpand TF.@require_interface(some(::AbstractFloat, b::Int))).args[2].head == Symbol("function")
            @test_broken (@macroexpand TF.@require_interface(some(::AbstractFloat{A}, b::Int) where {A})).args[2].head == Symbol("function")
        end
    end

    @testset "types.jl" begin
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
    end

    @testset "interfaces" begin
        C_4D = TF.SpatialArray(Ones(40, 40), 20u"nm")

        struct TF_1 <: TF.TransferFunction end
        tf_1 = TF_1()
        @test_throws ["does not implement", r"transfer(.*::TransferFunction, .*::SpatialArray.*)"] transfer(tf_1, C_4D)
        @test_throws ["does not implement", r"restore(.*::TransferFunction, .*::SpatialArray.*)"] restore(tf_1, C_4D)

        struct LTF_1 <: TF.LinearTransferFunction end
        ltf_1 = LTF_1()
        @test_throws ["does not implement", r"conv(.*::LinearTransferFunction, .*::SpatialArray.*)"] TF.conv(ltf_1, C_4D)
        @test_throws ["does not implement", r"deconv(.*::LinearTransferFunction, .*::SpatialArray.*)"] TF.deconv(ltf_1, C_4D)

        struct PSF_1 <: TF.PointSpreadFunction end
        psf_1 = PSF_1()
        @test_throws ["does not implement", r"intensity(.*::PointSpreadFunction, .*::Length, .*::Length)"] intensity(psf_1, 10u"nm", 10u"nm")
        @test_throws MethodError intensity(psf_1, 10u"nm")

        struct PSF_2 <: TF.RadialPSF end
        psf_2 = PSF_2()
        @test_throws ["does not implement", r"intensity(.*::RadialPSF, .*::Length)"] intensity(psf_2, 10u"nm")

        struct OTF_1 <: TF.OpticalTransferFunction end
        otf_1 = OTF_1()
        @test_throws ["does not implement", r"attenuation(.*::OpticalTransferFunction, .*::Frequency, .*::Frequency)"] attenuation(otf_1, 10u"nm^-1", 10u"nm^-1")

        struct OTF_2 <: TF.RadialOTF end
        otf_2 = OTF_2()
        @test_throws ["does not implement", r"attenuation(*::RadialOTF, .*::Frequency)"] attenuation(otf_2, 10u"nm^-1")
    end

    @testset "Apodization" begin
        using TransferFunctions.Apodization
        using TransferFunctions: Apodization as Apo
        # using TransferFunctions: Bartlett

        @testset "Types" begin
            @test Apo.Triangular() isa Apo.ApodizationFunction
            @test Apo.Welch() isa Apo.ApodizationFunction
            @test Apo.Connes() isa Apo.ApodizationFunction
            @test_throws ArgumentError Apo.PowerCosine{1 + 1im}()
            @test Apo.Cosine() isa Apo.PowerCosine{1}
            @test Apo.Hann() isa Apo.PowerCosine{2}
            @test Apo.Hamming() isa Apo.SineSum{2}
            @test Apo.Nuttall() isa Apo.SineSum{4}
            @test Apo.BlackmanNuttall() isa Apo.SineSum{4}
            @test Apo.BlackmanHarris() isa Apo.SineSum{4}
            @test Apo.FlatTop() isa Apo.SineSum{5}
            @test_throws ArgumentError Apo.Blackman(0.5)
            @test_throws ArgumentError Apo.Blackman(0.6)
            @test_throws ArgumentError Apo.Blackman(0)
            @test Apo.Blackman() == Apo.Blackman(0.16)
            @test Apo.Blackman() isa Apo.ApodizationFunction
            @test Apo.ExactBlackman() isa Apo.Blackman
            @test Apo.Gaussian(rand(Uniform(0.01, 0.49))) isa Apo.ApodizationFunction
            @test_throws ArgumentError Apo.Gaussian(0.6)
        end


        # FIX: Most of the apodization / window functions that are defined, should be 0 at the edge (for me the edge
        # is in 1 for implementation reasons) and all should be 1 in the center. That is at 0 in my implementation.
        # <02-09-24> 
        @testset "methods $(typeof(apo))" for apo in [
            Apo.Hamming(), Apo.Hann(), Apo.Welch(), Apo.Connes(), Apo.Cosine(),
            Apo.Nuttall(), Apo.BlackmanNuttall(), Apo.BlackmanHarris(), Apo.FlatTop(), Apo.ExactBlackman(),
            Apo.Blackman(0.4), Apo.Gaussian(0.4),
        ]
            @test Apo.apodization(apo, 0) == 1

            r = rand()
            @test Apo.apodization(apo, r) ≈ Apo.apodization(apo, -r)

            if typeof(apo) ∈ (Apo.Hann, Apo.Triangular, Apo.Cosine, Apo.Connes, Apo.Welch)
                @test Apo.apodization(apo, 1) == Apo.apodization(apo, -1) == 0
            end
        end

        @testset "equavalence $(equivs)" for equivs in [
            (Apo.Hann(), Apo.SineSum((0.5, 0.5))),
            (Apo.PowerCosine{0}(), Apo.SineSum((1.0,))),
            (Apo.PowerCosine{2}(), Apo.SineSum((0.5, 0.5))),
            (Apo.PowerCosine{4}(), Apo.SineSum((0.375, 0.5, 0.125))),
            (Apo.PowerCosine{6}(), Apo.SineSum((0.3125, 0.46875, 0.1875, 0.03125))),
        ]
            @test Apo.apodization.(equivs[1], -1:0.01:1) ≈ Apo.apodization.(equivs[2], -1:0.01:1)
        end

        @testset "internals" begin
            O = Ones(30, 30, 30)
            @test Apo.specify_border(O, IF.Fill(5.0), ((5, 5, 5), (2, 2, 2)), (1, 2, 3)) == IF.Fill(5.0, (5, 5, 5), (2, 2, 2))
            @test Apo.specify_border(O, IF.Fill(5.0), ((5, 5), (2, 2)), (1, 3)) == IF.Fill(5.0, (5, 0, 5), (2, 0, 2))
            @test Apo.specify_border(O, IF.Pad(:replicate), ((5,), (2,)), (1,)) == IF.Pad(:replicate, (5, 0, 0), (2, 0, 0))
            @test Apo.specify_border(O, IF.Inner(), ((5,), (2,)), (1,)) == IF.Inner((5, 0, 0), (2, 0, 0))
            @test_throws ArgumentError Apo.specify_border(O, IF.NoPad(), ((5,), (2,)), (1,)) == IF.Pad(:replicate, (5, 0, 0), (2, 0, 0))
        end

        @testset "methods" begin
            apo = Apo.Cosine()
            @test taperedges(
                      apo, ones(100, 100, 9), ((10, 10), (10, 10)); dims=(1, 2) # All the supplied arguments 
                  ) == taperedges(
                      apo, ones(100, 100, 9), ((10, 10), (10, 10)) # Infer dims
                  ) == taperedges(
                      apo, ones(100, 100, 9), 10; dims=(1, 2) # Supply single unified width
                  ) == taperedges(
                      apo, ones(100, 100, 9), (10, 10); dims=(1, 2) # Supply single width for each dimension
                  ) == taperedges(
                      apo, ones(100, 100, 9), (10, 10) # Infer dims
                  ) == taperedges(
                      ones(100, 100, 9), (10, 10) # default `apo`
                  ) == taperedges(
                      ones(100, 100, 9), (10, 10), "replicate"; dims=(1, 2) # default border 
                  ) == taperedges(
                      ones(100, 100, 9), (10, 10), IF.Pad{0}(:replicate, (), ())  # instantiate border
                  ) == taperedges(
                      ones(100, 100, 9), (10, 10), IF.Pad{3}(:replicate, (10, 10, 0), (10, 10, 0))  # Correct border size
                  )

            @test_throws DimensionMismatch taperedges(apo, ones(100, 100), 10; dims=(1, 2, 3))
            @test_throws DimensionMismatch taperedges(apo, ones(100, 100), 10; dims=(1, 3))

            # `dims`
            @test taperedges(ones(30, 30), 10; dims=:) isa AbstractArray
            @test taperedges(ones(30, 30), 10; dims=2) isa AbstractArray
            @test taperedges(ones(30, 30), 10; dims=Dims((1, 2))) isa AbstractArray

            @test size(taperedges(ones(30, 30), (10, 20))) == (50, 70)

            A_tap = taperedges(Apo.Cosine(), ones(30, 30), 10)

            # RESEARCH: Should this in fact be 0 at the other edge as well? This is done to make the signal periodic
            # so if we taper one edge to 0 and the other edge to the 0 length - 1, we will already have a periodic
            # signal, correct? <10-09-24> 
            @test all(all.([
                (A_tap[-9, :] .== 0),
                (A_tap[40, :] .== 0),
                (A_tap[:, -9] .== 0),
                (A_tap[:, 40] .== 0),
                (A_tap[1:30, 1:30] .== 1),
                (A_tap[-9:0, -9:0] .!= 1),
                (A_tap[31:40, -9:0] .!= 1),
                (A_tap[-9:0, 31:40] .!= 1),
                (A_tap[31:40, 31:40] .!= 1)
            ]))
        end
    end

    @testset "Optical transfer functions" begin
        A = Ones(1024, 1024)

        @testset "OTFArray" begin
            # FIX: @test_throws DomainError OTFArray(ones(3, 3, 3), 32u"nm", (4, 1, 1))
            # FIX: @test OTFArray(ones(3, 3), 32u"nm") isa OTFArray{<:Real}
            # FIX: @test OTFArray(ones(3, 3, 3), 32u"nm") isa OTFArray{<:Real,3}

            # FIX: tf = CircularPupilOTF(488u"nm", 1.4, 1.0, 0.3)
            # FIX: otf_array = otf(tf, 64u"nm", (512, 512))
            # FIX: tf_1 = OTFArray(otf_array, 64u"nm", (2, 3))

            # FIX: @test otf_array == otf(tf_1, 64u"nm", (512, 512))
            # FIX: @test argmax(FFTView(otf(tf_1, (512, 512)))) == CartesianIndex(1, 2)
            # FIX: @test argmax(FFTView(otf(tf_1, (511, 511)))) == CartesianIndex(1, 2)
            # FIX: @test otf(s_tf_1, (512, 512)) ≈ FourierTools.shift(otf(s_tf, (512, 512)), (1, 2))
            # FIX: @test otf(s_tf_1, (511, 511)) ≈ FourierTools.shift(otf(s_tf, (511, 511)), (1, 2))

            # TODO: When there is a model that returns Complex, i.e. the model with a phase shift <26-08-24> 
            # FIX: @test_skip otf(s_tf_1, (512, 512)) ≈ FourierTools.shift(otf(ComplexF32, s_tf, (512, 512)), (1, 2))
            # FIX: @test_skip otf(s_tf_1, (511, 511)) ≈ FourierTools.shift(otf(ComplexF32, s_tf, (511, 511)), (1, 2))

            # FIX: s_tf_2 = SampledOTF(tf, 64u"nm", (1.5, -2.5))
            # FIX: @test otf(s_tf_2, (512, 512)) ≈ real(FourierTools.shift(otf(ComplexF32, s_tf, (512, 512)), (0.5, -3.5)))
            # FIX: @test otf(s_tf_2, (511, 511)) ≈ real(FourierTools.shift(otf(ComplexF32, s_tf, (511, 511)), (0.5, -3.5)))

            # FIX: Probably a numerical error of FFT in the FourierTools package, but not sure <26-08-24> 
            # FIX: @test_broken otf(s_tf_2, (512, 512)) ≈ real(FourierTools.shift(otf(s_tf, (512, 512)), (0.5, -3.5)))
            # FIX: @test_broken otf(s_tf_2, (511, 511)) ≈ real(FourierTools.shift(otf(s_tf, (511, 511)), (0.5, -3.5)))
        end

        @testset "OTF models" begin
            for tf in (
                CircularPupilOTF(488u"nm", 1.4, 1.0, 0.3),
            )
                @testset "$(nameof(typeof(tf)))" begin
                    @test attenuation(tf, 1 // 250u"nm", 1 // 200u"nm") isa AbstractFloat
                    @test attenuation(tf, 250.0u"nm^-1", 200.0u"nm^-1") isa AbstractFloat
                    if tf isa TF.RadialOTF
                        @test attenuation(tf, 1 // 250u"nm") isa Number

                        c = 1 / 200u"nm"
                        @test attenuation(tf, c, 0u"nm^-1") == attenuation(tf, c) ≈
                              attenuation(tf, c / sqrt(10), 3c / sqrt(10)) ≈ attenuation(tf, c / sqrt(2), c / sqrt(2))

                        @test cutoff(tf) isa Frequency
                        @test cutoff(tf, 0.15) isa Frequency

                        ρ_max = cutoff(tf)
                        @test TF.insupport(tf, ρ_max / sqrt(2) / 2) == true
                        @test TF.insupport(tf, ρ_max / 2sqrt(2), ρ_max / 2sqrt(2)) == true
                        @test TF.insupport(tf, ρ_max / sqrt(10), 3ρ_max / 2sqrt(10)) == true
                    end
                    # FIX: @test attenuation(Float32, tf, 1 // 250u"nm") isa Float32
                    # FIX: @test attenuation(ComplexF32, tf, 1 // 250u"nm") isa ComplexF32
                end
            end
        end

        @testset "Sampled OTF" begin
            tf = CircularPupilOTF(488u"nm", 1.4, 1.0, 0.3)
            s_img = SpatialArray(A, 32u"nm")

            ## Methods - Array generation
            @test otf(tf, 60u"nm", (512, 512)) isa Matrix
            @test otf(tf, s_img) isa Matrix
            # FIX: @test otf(ComplexF32, s_tf, (512, 512)) isa Matrix{ComplexF32}

            # NOTE: zeroth frequency is the greatest <26-08-24> 
            @test argmax(FFTView(otf(tf, 60u"nm", (512, 512)))) == CartesianIndex(0, 0)
            @test FFTView(otf(tf, 60u"nm", (512, 512)))[0, 0] == 1

            # FIX: using TransferFunctions: support
            # NOTE: The support calculation should match the generation of the array <26-08-24>
            # FIX: @test support(s_tf, img) isa BitArray
            # FIX: @test (otf(s_tf, (512, 512)) .> 0) == support(s_tf, (512, 512))
            # FIX: @test (otf(s_tf, (512, 512)) .>= 0.15) == support(s_tf, (512, 512), a=0.15)
            # FIX: @test (otf(s_tf_1, (512, 512)) .> 0) == support(s_tf_1, (512, 512))
            # FIX: @test (otf(s_tf_1, (512, 512)) .>= 0.15) == support(s_tf_1, (512, 512), a=0.15)
            # FIX: @test (otf(s_tf_2, (512, 512)) .> 0) == support(s_tf_2, (512, 512))
            # FIX: @test (otf(s_tf_2, (512, 512)) .>= 0.15) == support(s_tf_2, (512, 512), a=0.15)
            # FIX: @test count(support(s_tf, (512, 512), a=0.15)) < count(support(s_tf, (512, 512)))
            # FIX: @test all(otf(s_tf, (512, 512))[support(s_tf, (512, 512)).==false] .== 0)

            # FIX: using TransferFunctions: overlap
            # FIX: @test overlap(s_tf, s_tf, img) isa BitArray
            # FIX: @test overlap(s_tf, s_tf, (512, 512)) == support(s_tf, (512, 512))
            # FIX: @test overlap(s_tf, s_tf, (512, 512); a_1=0.15) == support(s_tf, (512, 512), a=0.15)
            # FIX: @test count(overlap(s_tf_1, s_tf, (512, 512))) < count(support(s_tf, (512, 512)))

            ## Non-methods - Array generation
            @test_throws MethodError otf(tf, 64u"nm", (512.1, 512.4))
        end
    end

    @testset "Point Spread Function" begin
        A = Ones(1024, 1024)

        @testset "MeasuredPSF" begin
            # FIX: @test_throws DomainError PSFArray(ones(3, 3, 3), 32u"nm", (4, 1, 1))
            # FIX: @test PSFArray(ones(3, 3), 32u"nm") isa MeasuredPSF{<:Real,2}
            # FIX: @test MeasuredPSF(ones(3, 3, 3), 32u"nm") isa MeasuredPSF{<:Real,3}
        end

        @testset "PSFArray" begin end

        @testset "ModelPSF" begin
            tf = BornWolf(488u"nm", 1.4, 1.7)

            ## Method Availability
            @test intensity(tf, 250u"nm", 200u"nm") isa Number
            @test intensity(tf, 250u"nm") isa Number
        end

        @testset "Sampled PSF" begin
            s_img = SpatialArray(A, 32u"nm")
            tf = BornWolf(488u"nm", 1.4, 1.7)

            ## Method Availability - Construction
            @test psf(tf, 64u"nm", (512, 512)) isa AbstractMatrix
            @test psf(tf, (64u"nm", 32u"nm"), (512, 512)) isa AbstractMatrix

            ## Methods - Array generation
            @test psf(tf, 60u"nm", (512, 512)) isa OAs.OffsetMatrix
            @test_throws MethodError psf(tf, s_img)

            ## Non-methods - Array generation
            @test_throws MethodError psf(tf, 60u"nm", (512.1, 512.4))
        end
    end

    @testset "Estimation" begin
        using TransferFunctions.Estimation
        @test (bead(100u"nm", 30.5u"nm", intensity=0.75) .<= 0.75) |> all
        @test_throws AssertionError bead(100u"nm", 30.5u"nm"; position=(-1.5, 0.7))
        @test_throws AssertionError bead(100u"nm", 30.5u"nm"; position=(1.5, 0.7))

        @test (bead(100u"nm", 30.5u"nm"; position=(0.5, 0.5)) .== reverse(bead(100u"nm", 30.5u"nm"; position=(-0.5, -0.5)))) |> all
        @test (2bead(100u"nm", 30.5u"nm"; intensity=0.5) .== bead(100u"nm", 30.5u"nm")) |> all
    end

    @testset "Base" begin
        psf1 = BornWolf(488u"nm", 1, 1.7)
        psf2 = BornWolf(488u"nm", 1.0, 1.7)

        @test psf1 == psf2
        @test hash(psf1) == hash(psf2)
        @test isequal(psf1, psf2)
    end
end
