using TransferFunctions
using TransferFunctions: Frequency
using FillArrays
using FourierTools
using FFTViews
using Aqua, Test, Documenter

@testset "TransferFunctions.jl" begin
    if haskey(ENV, "RUNTESTS_FULL") || haskey(ENV, "GITHUB_ACTIONS")
        @testset "Code quality (Aqua.jl)" begin
            Aqua.test_all(
                TransferFunctions;
                ambiguities=(; exclude=VERSION >= v"1.11" ? [checkindex, checkbounds] : [])
                # ambiguities=VERSION >= v"1.1" ? (; broken=true) : false
            )
        end
    else
        @info "Skipping Aqua.jl quality tests. For a full run set `ENV[\"RUNTESTS_FULL\"]=true`."
    end

    # FIX: When running locally, do not ask for SSH key password <10-12-23> 
    if haskey(ENV, "RUNTESTS_FULL") && (!haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "RUNNER_OS") && ENV["RUNNER_OS"] == "Linux")
        @testset "DocTests" begin
            # NOTE: Better than doc-testing in `make.jl` because, I can track the coverage
            DocMeta.setdocmeta!(TransferFunctions, :DocTestSetup, :(using TransferFunctions); recursive=true)
            doctest(TransferFunctions)
        end
    else
        @info "Skipping Documenter.jl doctests. For a full run set `ENV[\"RUNTESTS_FULL\"]=true`."
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
