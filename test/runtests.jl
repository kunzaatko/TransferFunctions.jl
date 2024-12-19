using TransferFunctions
using TransferFunctions: Frequency
using FillArrays, FourierTools, IntervalSets, FFTViews, Distributions
using Aqua, Test, Documenter

@testset "TransferFunctions.jl" begin
    @testset "Code quality" begin
        @testset "Aqua.jl" begin
            if haskey(ENV, "RUNTESTS_FULL") || haskey(ENV, "GITHUB_ACTIONS")
                Aqua.test_all(
                    TransferFunctions;
                    ambiguities=false
                    # ambiguities=VERSION >= v"1.1" ? (; broken=true) : false
                )
            else
                @info "Skipping Aqua.jl quality tests. For a full run set `ENV[\"RUNTESTS_FULL\"]=true`."
            end
        end
        @testset "Ambiguities" begin
            @test length(Test.detect_ambiguities(TransferFunctions)) == 0
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
                    using FFTW;
                    using Logging; # NOTE: This does not need to be in the `make.jl` of docs. We want `@warn ` to function there <19-12-24> 
                    Logging.disable_logging(Logging.Warn)
                ); recursive=true)
            doctest(TransferFunctions)
        end
    end
    @testset "utils.jl + types.jl" begin
        using TransferFunctions: PixelSize, Coordinate, Frequency, Length
        using OffsetArrays
        @test 1 / 32u"nm" isa Frequency
        Δkx = 1 / 61u"nm"
        @test 1 / Δkx isa Length
        @test (31.5u"nm", 40u"nm", 50u"nm") isa PixelSize{3}
        @test TransferFunctions.fillsize(31u"nm", 2) == (31u"nm", 31u"nm")
        @test (3, 3, 3) isa Coordinate{3}
        @test (3.0, 3.0, 3.0) isa Coordinate{3}
        @test (2.4, 3, 3) isa Coordinate{3} # NOTE: Must accept diverse types <26-08-24> 
        @test TransferFunctions.roundupcenter(ones(3, 4, 2)) == (2, 3, 2)
        @test TransferFunctions.roundupcenter(OffsetArray(ones(3, 3, 3), -2, -2, -2)) == Tuple(zeros(3))
        @test TransferFunctions.roundupcenter(OffsetArray(ones(4, 3, 3), -2, -2, -2)) == (1, 0, 0)
        @test TransferFunctions.exactcenter(ones(3, 4, 2)) == (2, 2.5, 1.5)
        @test TransferFunctions.exactcenter(OffsetArray(ones(3, 3, 3), -2, -2, -2)) == Tuple(zeros(3))
        @test TransferFunctions.exactcenter(OffsetArray(ones(4, 3, 3), -2, -2, -2)) == (0.5, 0, 0)
        @test TransferFunctions.contained(ones(3, 4, 2), (2, 3, 1))
        @test TransferFunctions.contained(ones(3, 4, 2), (8, 3, 1)) == false
        @test TransferFunctions.contained(OffsetArray(ones(3, 3, 3), -2, -2, -2), (-1, 1, 0))
        @test TransferFunctions.contained(OffsetArray(ones(3, 3, 3), -2, -2, -2), (-2, 1, 0)) == false

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
    end

    @testset "OTF" begin
        img = Ones(1024, 1024)

        @testset "MeasuredOTF" begin
            @test_throws DomainError MeasuredOTF(ones(3, 3, 3), 32u"nm", (4, 1, 1))
            @test MeasuredOTF(ones(3, 3), 32u"nm") isa MeasuredOTF{<:Real,2}
            @test MeasuredOTF(ones(3, 3, 3), 32u"nm") isa MeasuredOTF{<:Real,3}
        end

        @testset "ModelOTF" begin
            tf = IdealOTFwithCurvature(488u"nm", 1.4, 1.0, 0.3)

            ## Method Availability
            @test attenuation(tf, 1 // 250u"nm", 1 // 200u"nm") isa Number
            @test attenuation(tf, 1 // 250u"nm") isa Number # FIX: This should function only for a RadiallySymmetric psf <26-08-24> 
            @test attenuation(Float32, tf, 1 // 250u"nm") isa Float32
            @test attenuation(ComplexF32, tf, 1 // 250u"nm") isa ComplexF32

            @test cutoff(tf) isa Frequency
            @test cutoff(tf, 0.15) isa Frequency
        end

        @testset "SampledOTF" begin
            tf = IdealOTFwithCurvature(488u"nm", 1.4, 1.0, 0.3)
            psf_tf = BornWolf(488u"nm", 1.4, 1.7)

            ## Method Availability
            @test SampledOTF(tf, 64u"nm") isa SampledOTF # NOTE: Fill in the sizes <26-08-24> 
            @test SampledOTF(tf, (64u"nm", 32u"nm")) isa SampledOTF # NOTE: Non-isometric <26-08-24> 
            @test SampledOTF(tf, 64u"nm", (1, 2)) isa SampledOTF # NOTE: Non-centred <26-08-24> 
            @test SampledOTF(tf, (64u"nm", 32u"nm"), (1, 2)) isa SampledOTF
            @test SampledOTF(tf, 64u"nm", (1.5, 2)) isa SampledOTF # NOTE: Non-integer center <26-08-24> 
            @test SampledOTF(tf, (64u"nm", 32u"nm"), (1.5, 2)) isa SampledOTF

            ## Non-methods
            @test_throws MethodError SampledOTF(psf_tf, 64u"nm") # NOTE: PSF model should not work <26-08-24> 

            s_tf = SampledOTF(tf, 64u"nm")

            ## Methods - Array generation
            @test otf(s_tf, (512, 512)) isa Matrix
            @test otf(s_tf, img) isa Matrix
            @test otf(ComplexF32, s_tf, (512, 512)) isa Matrix{ComplexF32}

            ## Methods - Array generation Consistency
            @test otf(s_tf, img) == (mtf(s_tf, img) .+ im .* ptf(s_tf, img))

            # NOTE: zeroth frequency is the greatest <26-08-24> 
            @test argmax(FFTView(otf(s_tf, (512, 512)))) == CartesianIndex(0, 0)
            @test FFTView(otf(s_tf, (512, 512)))[0, 0] == 1

            s_tf_1 = SampledOTF(tf, 64u"nm", (2, 3))
            # NOTE: In terms of frequencies we have zero based indexing since we need to skip the zeroth frequency (0,0)
            # <26-08-24> 
            @test argmax(FFTView(otf(s_tf_1, (512, 512)))) == CartesianIndex(1, 2)
            @test argmax(FFTView(otf(s_tf_1, (511, 511)))) == CartesianIndex(1, 2)
            @test otf(s_tf_1, (512, 512)) ≈ FourierTools.shift(otf(s_tf, (512, 512)), (1, 2))
            @test otf(s_tf_1, (511, 511)) ≈ FourierTools.shift(otf(s_tf, (511, 511)), (1, 2))

            # TODO: When there is a model that returns Complex, i.e. the model with a phase shift <26-08-24> 
            @test_skip otf(s_tf_1, (512, 512)) ≈ FourierTools.shift(otf(ComplexF32, s_tf, (512, 512)), (1, 2))
            @test_skip otf(s_tf_1, (511, 511)) ≈ FourierTools.shift(otf(ComplexF32, s_tf, (511, 511)), (1, 2))

            s_tf_2 = SampledOTF(tf, 64u"nm", (1.5, -2.5))
            @test otf(s_tf_2, (512, 512)) ≈ real(FourierTools.shift(otf(ComplexF32, s_tf, (512, 512)), (0.5, -3.5)))
            @test otf(s_tf_2, (511, 511)) ≈ real(FourierTools.shift(otf(ComplexF32, s_tf, (511, 511)), (0.5, -3.5)))

            # FIX: Probably a numerical error of FFT in the FourierTools package, but not sure <26-08-24> 
            @test_broken otf(s_tf_2, (512, 512)) ≈ real(FourierTools.shift(otf(s_tf, (512, 512)), (0.5, -3.5)))
            @test_broken otf(s_tf_2, (511, 511)) ≈ real(FourierTools.shift(otf(s_tf, (511, 511)), (0.5, -3.5)))

            using TransferFunctions: support
            # NOTE: The support calculation should match the generation of the array <26-08-24>
            @test support(s_tf, img) isa BitArray
            @test (otf(s_tf, (512, 512)) .> 0) == support(s_tf, (512, 512))
            @test (otf(s_tf, (512, 512)) .>= 0.15) == support(s_tf, (512, 512), a=0.15)
            @test (otf(s_tf_1, (512, 512)) .> 0) == support(s_tf_1, (512, 512))
            @test (otf(s_tf_1, (512, 512)) .>= 0.15) == support(s_tf_1, (512, 512), a=0.15)
            @test (otf(s_tf_2, (512, 512)) .> 0) == support(s_tf_2, (512, 512))
            @test (otf(s_tf_2, (512, 512)) .>= 0.15) == support(s_tf_2, (512, 512), a=0.15)
            @test count(support(s_tf, (512, 512), a=0.15)) < count(support(s_tf, (512, 512)))
            @test all(otf(s_tf, (512, 512))[support(s_tf, (512, 512)).==false] .== 0)

            using TransferFunctions: overlap
            @test overlap(s_tf, s_tf, img) isa BitArray
            @test overlap(s_tf, s_tf, (512, 512)) == support(s_tf, (512, 512))
            @test overlap(s_tf, s_tf, (512, 512); a_1=0.15) == support(s_tf, (512, 512), a=0.15)
            @test count(overlap(s_tf_1, s_tf, (512, 512))) < count(support(s_tf, (512, 512)))

            ## Non-methods - Array generation
            @test_throws MethodError otf(s_tf, (512.1, 512.4)) # NOTE: non-integer image size <27-08-24>
            @test_throws MethodError otf(s_tf, 512) # NOTE: Do not infer size without information <26-08-24> 
        end
        # TODO: Add tests for the particular models <24-10-23> 
    end

    @testset "PSF" begin
        img = Ones(1024, 1024)

        @testset "MeasuredPSF" begin
            @test_throws DomainError MeasuredPSF(ones(3, 3, 3), 32u"nm", (4, 1, 1))
            @test MeasuredPSF(ones(3, 3), 32u"nm") isa MeasuredPSF{<:Real,2}
            @test MeasuredPSF(ones(3, 3, 3), 32u"nm") isa MeasuredPSF{<:Real,3}
        end

        @testset "ModelPSF" begin
            using OffsetArrays
            tf = BornWolf(488u"nm", 1.4, 1.7)

            ## Method Availability
            @test psf(tf, 250u"nm", 200u"nm") isa Number
            @test psf(tf, 250u"nm") isa Number # FIX: This should function only for a RadiallySymmetric psf <26-08-24> 
        end

        @testset "SampledPSF" begin
            using OffsetArrays
            tf = BornWolf(488u"nm", 1.4, 1.7)
            otf_tf = IdealOTFwithCurvature(488u"nm", 1.4, 1.0, 0.3)

            ## Method Availability - Construction
            @test SampledPSF(tf, 64u"nm") isa SampledPSF # NOTE: Fill in the sizes <26-08-24> 
            @test SampledPSF(tf, (64u"nm", 32u"nm")) isa SampledPSF # NOTE: Non-isometry in pixel-sizes <26-08-24> 
            @test SampledPSF(tf, 64u"nm", (1, 2)) isa SampledPSF # NOTE: Non-centred <26-08-24> 
            @test SampledPSF(tf, (64u"nm", 32u"nm"), (1, 2)) isa SampledPSF
            @test SampledPSF(tf, 64u"nm", (1.5, 2)) isa SampledPSF # NOTE: Non-integer center <26-08-24> 
            @test SampledPSF(tf, (64u"nm", 32u"nm"), (1.5, 2)) isa SampledPSF

            ## Non-methods - Construction
            @test_throws MethodError SampledPSF(otf_tf, 64u"nm") # NOTE: OTF model should not work <26-08-24> 

            s_tf = SampledPSF(tf, 64u"nm")

            ## Methods - Array generation
            @test psf(s_tf, (512, 512)) isa OffsetArrays.OffsetMatrix
            @test psf(s_tf, img) isa OffsetArrays.OffsetMatrix

            s_tf_1 = SampledPSF(tf, 64u"nm", (1, 2))
            @test argmax(psf(s_tf_1, img)) == CartesianIndex(1, 2)
            s_tf_2 = SampledPSF(tf, 64u"nm", (-1, -2))
            @test argmax(psf(s_tf_2, img)) == CartesianIndex(-1, -2)
            s_tf_3 = SampledPSF(tf, 64u"nm", (1.98, 1.98))
            @test_broken argmax(psf(s_tf_3, img)) == CartesianIndex(2, 2)

            ## Non-methods - Array generation
            # TODO: Add when shift in generating is implemented <24-10-23>  @test psf(tf, 512, 64u"nm"; δ=(2, 1)) isa Matrix
            @test_throws MethodError psf(s_tf, (512.1, 512.4)) # NOTE: non-integer image size <27-08-24>
            @test_throws MethodError psf(s_tf, 512) # NOTE: Do not infer size without information <26-08-24> 
        end
    end

    @testset "interfaces from Base" begin
        psf1 = BornWolf(488u"nm", 1, 1.7)
        psf2 = BornWolf(488u"nm", 1.0, 1.7)

        # TODO: test `isequal` on missing values and equivalence operator on missing values <19-11-23> 
        # hmissing = Harmonic(missing, π / 4, 2 / 61u"nm", 0)

        @test psf1 == psf2
        @test hash(psf1) == hash(psf2)
        @test isequal(psf1, psf2)
    end
end
