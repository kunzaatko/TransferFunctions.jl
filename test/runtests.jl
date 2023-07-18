using TransferFunctions
using TransferFunctions: Frequency
using FillArrays
using Test

@testset "TransferFunctions.jl" begin
    @testset "OTF" begin
        tf = IdealOTFwithCurvature(488u"nm", 1.4, 1.0, 0.3)
        img = Ones(1024, 1024)

        # otf
        ## Method Availability
        @test otf(tf, 1 // 250u"nm", 1 // 200u"nm") isa Number
        @test otf(tf, 1 // 250u"nm") isa Number
        @test otf(tf, 512, 64u"nm") isa Matrix
        @test_throws MethodError otf(tf, (512.1, 512.4), 64u"nm")

        ## Method Consistency
        @test otf(tf, 512, (64u"nm", 64u"nm")) == otf(tf, (512, 512), 64u"nm")
        @test otf(tf, img, 54u"nm") == otf(tf, img, (54u"nm", 54u"nm"))

        # otf_support
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

    @testset "PSF" begin
        # TODO
    end
end
