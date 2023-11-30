using TransferFunctions
using TransferFunctions: Frequency
using FillArrays
using FourierTools
using Aqua, Test

@testset "TransferFunctions.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(TransferFunctions;
            # https://github.com/JuliaArrays/FillArrays.jl/issues/105#issuecomment-1582516319
            ambiguities=VERSION >= v"1.1" ? (; broken=true) : false
        )
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
            @test otf(tf, 1 // 250u"nm", 1 // 200u"nm") isa Number
            @test otf(tf, 1 // 250u"nm") isa Number
            @test otf(tf, 512, 64u"nm") isa Matrix
            @test otf(tf, 512, 64u"nm"; δ=(2, 1)) isa Matrix
            @test otf(tf, 512, 64u"nm"; δ=(2, 1)) == (mtf(tf, 512, 64u"nm"; δ=(2, 1)) .+ im .* ptf(tf, 512, 64u"nm"; δ=(2, 1)))
            @test_throws MethodError otf(tf, (512.1, 512.4), 64u"nm") # NOTE: non-integer image size not possible

            ## Method Consistency
            @test otf(tf, 512, (64u"nm", 64u"nm")) == otf(tf, (512, 512), 64u"nm")
            @test otf(tf, img, 54u"nm") == otf(tf, img, (54u"nm", 54u"nm"))

            @testset "SIM utils" begin
                ## Consistent shift (FourierTools)
                # NOTE: `shift` using Fourier shift theorem not interpolation... When interpolation is used this must be done
                # with centred data
                @test otf(tf, 512, 61u"nm"; δ=(-1, 2)) ≈ real.(FourierTools.shift(ComplexF32.(otf(tf, 512, 61u"nm")), (-1, 2)))
                @test otf(tf, 511, 61u"nm"; δ=(-1, 2)) ≈ real.(FourierTools.shift(ComplexF32.(otf(tf, 511, 61u"nm")), (-1, 2)))
                @test otf(tf, 512, 61u"nm"; δ=(-1.5, 2.5)) ≈ real.(FourierTools.shift(ComplexF32.(otf(tf, 512, 61u"nm")), (-1.5, 2.5)))
                @test otf(tf, 511, 61u"nm"; δ=(-1.5, 2.5)) ≈ real.(FourierTools.shift(ComplexF32.(otf(tf, 511, 61u"nm")), (-1.5, 2.5)))

                using TransferFunctions: otf_support

                ## Method Availability
                @test otf_support(tf, 512, 64u"nm") isa BitMatrix
                @test otf_support(tf, 512, 64u"nm"; ρ=0.5) isa BitMatrix
                @test otf_support(tf, 512, 64u"nm"; ρ=-0.5) isa BitMatrix
                @test otf_support(tf, 512, 64u"nm"; ρ=(0.5, 0.7)) isa BitMatrix
                @test_throws MethodError otf_support(tf, (512.1, 512.4), 64u"nm")

                ## Method Consistency
                @test otf_support(tf, 512, (64u"nm", 64u"nm")) == otf_support(tf, (512, 512), 64u"nm")
                @test otf_support(tf, img, 54u"nm") == otf_support(tf, img, (54u"nm", 54u"nm"))

                ## Theory Consistency
                @test all(otf_support(tf, 512, 64u"nm"; ρ=-0.5) + otf_support(tf, 512, 64u"nm"; ρ=0.5) + otf_support(tf, 512, 64u"nm") .!= 1)
                @test otf_support(tf, 512, 64u"nm"; ρ=-0.5) .+ otf_support(tf, 512, 64u"nm"; ρ=0.5) == otf_support(tf, 512, 64u"nm")
                @test dropdims(any(isone,
                        cat([otf_support(tf, 512, 64u"nm"; ρ=ρ_int) for ρ_int in [0.1, (0.1, 0.5), (0.5, 0.7), -0.3]]..., dims=3),
                        dims=3),
                    dims=3) == otf_support(tf, 512, 64u"nm")
                @test all(otf(tf, 512, 64u"nm")[otf_support(tf, 512, 64u"nm").!=1] .== 0)

                # cutoff_frequency
                @test cutoff_frequency(tf) isa Frequency
            end
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
            @test psf(tf, 250u"nm") isa Number
            @test psf(tf, 512, 64u"nm") isa OffsetArrays.OffsetMatrix
            # TODO: Add when shift in generating is implemented <24-10-23>  @test psf(tf, 512, 64u"nm"; δ=(2, 1)) isa Matrix
            @test_throws MethodError psf(tf, (512.1, 512.4), 64u"nm") # NOTE: non-integer image size not possible

            ## Method Consistency
            @test psf(tf, 512, (64u"nm", 64u"nm")) == psf(tf, (512, 512), 64u"nm")
            @test psf(tf, img, 54u"nm") == psf(tf, img, (54u"nm", 54u"nm"))
        end
    end
end
