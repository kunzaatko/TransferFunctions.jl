using TransferFunctions
using TransferFunctions: TransferFunctions as TF
using TransferFunctions: Frequency
using IntervalSets, FourierTools, FFTViews, Distributions, FillArrays, TensorOperations, OffsetArrays
using OffsetArrays: OffsetArray as OA
using OffsetArrays: OffsetArrays as OAs
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
                @info "Skipping Aqua.jl quality tests. For a full run set `ENV[\"RUNTESTS_FULL\"]=true`."
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
        if !haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "RUNNER_OS") && ENV["RUNNER_OS"] == "Linux"
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
        end
    end

    @testset "utils.jl" begin
        using TransferFunctions: fillsize, roundupcenter, exactcenter, fftfreqs, posgrid, contained
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
            @test (3.0, 3.0, 3.0) isa Coordinate{3}
            @test (2.4, 3, 3) isa Coordinate{3} # NOTE: Must accept diverse types <26-08-24> 
        end


        @testset "utility functions" begin
            @test fillsize(31u"nm", 2) == (31u"nm", 31u"nm")
            @test fillsize(31u"nm", Val(2)) == (31u"nm", 31u"nm")

            @test (posgrid((11, 11), (31u"nm", 31u"nm")) |> first |> first) isa Length
            @test posgrid((11, 11), 31u"nm") == posgrid((11, 11), (31u"nm", 31u"nm"))
            @test_throws MethodError posgrid((11, 11, 11), (31u"nm", 31u"nm"))

            @test (fftfreqs((11, 11), (31u"nm", 31u"nm")) |> first |> first) isa Frequency
            @test fftfreqs((11, 11), 31u"nm") == fftfreqs((11, 11), (31u"nm", 31u"nm"))
            @test_throws MethodError fftfreqs((11, 11, 11), (31u"nm", 31u"nm"))

            @test roundupcenter(Ones(3, 4, 2)) == CI(2, 3, 2)
            @test roundupcenter(OA(Ones(3, 3, 3), -2, -2, -2)) == CI(0, 0, 0)
            @test roundupcenter(OA(Ones(4, 3, 3), -2, -2, -2)) == CI(1, 0, 0)

            @test exactcenter(Ones(3, 4, 2)) == (2.0, 2.5, 1.5)
            @test exactcenter(OA(Ones(3, 3, 3), -2, -2, -2)) == (0.0, 0.0, 0.0)
            @test exactcenter(OA(Ones(4, 3, 3), -2, -2, -2)) == (0.5, 0.0, 0.0)

            @test contained(Ones(3, 4, 2), (2, 3, 1))
            @test !contained(Ones(3, 4, 2), (8, 3, 1))
            @test contained(OA(Ones(3, 3, 3), -2, -2, -2), (-1, 1, 0))
            @test !contained(OA(Ones(3, 3, 3), -2, -2, -2), (-2, 1, 0))

            @test OriginAt(CI(2, 2, 2))(Ones(3, 3, 3)) == OA(Ones(3, 3, 3), -2, -2, -2)
        end

        @testset "macros" begin
            @test_throws LoadError @macroexpand TF.@require_interface(function some() end) # No arguments
            @test_throws LoadError @macroexpand TF.@require_interface(some())              # No arguments
            @test_throws LoadError @macroexpand TF.@require_interface(some = 5)            # Not a call
            @test_throws LoadError @macroexpand TF.@require_interface(some(a))             # No type
            @test_throws LoadError @macroexpand TF.@require_interface(some(a::Float64))    # Concrete type

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

        @testset "circulant" begin
            using Base: OneTo
            O = Ones(100, 100, 100)

            # Mismatch in `innerdims` length and inner array size.
            @test_throws MethodError TF.OuterInnerArray([view(O, 30:40, 50:60, :), view(O, 20:30, 40:50, :)], (2, 3))

            # Total dimensionality does not match partial dimensionalities
            @test_throws DimensionMismatch TF.OuterInnerArray{Any,3}([view(O, 30:40, 50:60, :), view(O, 20:30, 40:50, :)], (2, 3, 4))

            # inner axes mismatch
            @test_throws DimensionMismatch TF.OuterInnerArray([view(O, 30:39, 50:60, :), view(O, 20:30, 40:50, :)], (2, 3, 4))

            # Non existent dimension
            @test_throws DimensionMismatch TF.OuterInnerArray([view(O, 30:40, 50:60, :), view(O, 20:30, 40:50, :)], (2, 3, 5))

            @test TF.OuterInnerArray([view(O, 30:40, 50:60, :), view(O, 20:30, 40:50, :)], (2, 3, 4)) isa AbstractArray{<:Any,4}
            @test TF.OuterInnerArray{Float64,4}([view(O, 30:40, 50:60, :), view(O, 20:30, 40:50, :)], (2, 3, 4)) isa AbstractArray{Float64,4}

            # Explicit recasting
            @test !(TF.OuterInnerArray{Any,4}([view(O, 30:40, 50:60, :), view(O, 20:30, 40:50, :)], (2, 3, 4)) isa TF.OuterInnerArray{Float64})

            oia = TF.OuterInnerArray([view(O, 30:40, 50:60, :), view(O, 20:30, 40:50, :)], (2, 3, 4))
            @test size(oia) == (2, 11, 11, 100)
            @test axes(oia) == (OneTo(2), OneTo(11), OneTo(11), OneTo(100))

            @test TF.innerdims(oia) == (2, 3, 4)
            @test TF.outerdims(oia) == (1,)

            # TODO: Test offset array compatibility <16-04-25> 

            K = OAs.centered(zeros(21, 21))
            K[0, 0] = 0.5
            K[-1, -1] = K[-1, 1] = K[1, -1] = K[1, 1] = 0.5 / 4

            img = Ones(30, 30)

            @test TF.CirculantTensor(img, (10, 10)) isa TF.CirculantTensor{eltype(img),4}
            @test TF.CirculantTensor(img, (-5:5, -5:5)) isa TF.CirculantTensor{eltype(img),4}
            @test TF.CirculantTensor(img, K) isa TF.CirculantTensor{eltype(img),4}

            # Different eltypes
            @test TF.CirculantTensor(Ones{Int}(30, 30), K) isa TF.CirculantTensor{Int,4}

            A = TF.CirculantTensor(img, K)

            # Correct output indices
            @test ndims(TF.CirculantTensor(img, K)) == 4
            @test axes(A)[3:4] == axes(K) == A.kern
            @test axes(A)[1:2] == A.interior

            @tensor B[a, b] := OAs.no_offset_view(A)[a, b, c, d] * OAs.no_offset_view(K)[c, d]
            @test B isa AbstractMatrix
            @test size(B) == length.(A.interior)
        end
    end

    @testset "interfaces" begin
        A = TF.SpatialArray(Ones(40, 40), 20u"nm")

        struct TF_1 <: TF.TransferFunction end
        tf_1 = TF_1()
        @test_throws ["does not implement", r"transfer(.*::TransferFunction, .*::SpatialArray.*)"] transfer(tf_1, A)
        @test_throws ["does not implement", r"restore(.*::TransferFunction, .*::SpatialArray.*)"] restore(tf_1, A)

        struct LTF_1 <: TF.LinearTransferFunction end
        ltf_1 = LTF_1()
        @test_throws ["does not implement", r"conv(.*::LinearTransferFunction, .*::SpatialArray.*)"] TF.conv(ltf_1, A)
        @test_throws ["does not implement", r"deconv(.*::LinearTransferFunction, .*::SpatialArray.*)"] TF.deconv(ltf_1, A)

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
        using TransferFunctions: apodization, instrument, apodize, taperedges, Apodization
        using TransferFunctions: Triangular, Blackman, Connes, Cosine, Gaussian, Hamming, Hann, Welch, PowerCosine, SineSum, Nuttall, BlackmanNuttall, BlackmanHarris, FlatTop, ExactBlackman
        # using TransferFunctions: Bartlett

        @testset "Types" begin
            @test Triangular() isa Apodization
            @test Welch() isa Apodization
            @test Connes() isa Apodization
            @test_throws ArgumentError PowerCosine{1 + 1im}()
            @test Cosine() isa PowerCosine{1}
            @test Hann() isa PowerCosine{2}
            @test Hamming() isa SineSum{2}
            @test Nuttall() isa SineSum{4}
            @test BlackmanNuttall() isa SineSum{4}
            @test BlackmanHarris() isa SineSum{4}
            @test FlatTop() isa SineSum{5}
            @test_throws ArgumentError Blackman{0.5}()
            @test_throws ArgumentError Blackman{0.6}()
            @test_throws ArgumentError Blackman{0}()
            @test Blackman() == Blackman{0.16}()
            @test Blackman() isa Apodization
            @test ExactBlackman() isa Blackman
            @test Gaussian(rand(Uniform(0.01, 0.49))) isa Apodization
            @test_throws ArgumentError Gaussian(0.6)
        end


        # FIX: Most of the apodization / window functions that are defined, should be 0 at the edge (for me the edge
        # is in 1 for implementation reasons) and all should be 1 in the center. That is at 0 in my implementation.
        # <02-09-24> 
        @testset "methods $(typeof(apo))" for apo in [
            Hamming(), Hann(), Welch(), Connes(), Cosine(),
            Nuttall(), BlackmanNuttall(), BlackmanHarris(), FlatTop(), ExactBlackman(),
            Blackman{0.4}(), Gaussian(0.4),
        ]
            @test apodization(apo, 0) == 1

            r = rand()
            @test apodization(apo, r) ≈ apodization(apo, -r)

            if typeof(apo) ∈ (Hann, Triangular, Cosine, Connes, Welch)
                @test apodization(apo, 1) == apodization(apo, -1) == 0
            end
        end

        @testset "equavalence $(equivs)" for equivs in [
            (Hann(), SineSum{2,(0.5, 0.5)}()),
            (PowerCosine{0}(), SineSum{1,(1,)}()),
            (PowerCosine{2}(), SineSum{2,(0.5, 0.5)}()),
            (PowerCosine{4}(), SineSum{3,(0.375, 0.5, 0.125)}()),
            (PowerCosine{6}(), SineSum{4,(0.3125, 0.46875, 0.1875, 0.03125)}()),
        ]
            @test apodization.(equivs[1], -1:0.01:1) ≈ apodization.(equivs[2], -1:0.01:1)
        end

        @testset "methods" begin
            using ImageFiltering: Pad
            apo = Cosine()
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
                      ones(100, 100, 9), (10, 10), Pad{0}(:replicate, (), ())  # instantiate border
                  ) == taperedges(
                      ones(100, 100, 9), (10, 10), Pad{3}(:replicate, (10, 10, 0), (10, 10, 0))  # Correct border size
                  )

            # TODO: Test with various padding edges... Does the arrays size match? <09-09-24> 

            @test_throws ArgumentError taperedges(apo, ones(100, 100), 10; dims=(1, 2, 3))
            @test_throws ArgumentError taperedges(apo, ones(100, 100), 10; dims=(1, 3))

            # `dims`
            @test_broken taperedges(ones(30, 30), 10; dims=:) isa AbstractArray
            @test_broken taperedges(ones(30, 30), 10; dims=2) isa AbstractArray
            @test taperedges(ones(30, 30), 10; dims=Dims((1, 2))) isa AbstractArray

            @test size(taperedges(ones(30, 30), (10, 20))) == (50, 70)

            A_tap = taperedges(Cosine(), ones(30, 30), 10)

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
        img = Ones(1024, 1024)

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
            s_img = SpatialArray(img, 32u"nm")

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
        img = Ones(1024, 1024)

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
            s_img = SpatialArray(img, 32u"nm")
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
