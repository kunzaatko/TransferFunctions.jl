using TransferFunctions: TransferFunctions as TF
using OffsetArrays: OffsetArrays as OAs

A = ones(100, 100)

λ = 488u"nm"
NA = 1.4
n = 1.5
Δxy = 64u"nm"
Δx = 64u"nm"
Δy = 32u"nm"
Δz = 15u"nm"
f = 1.3

PSF_2D = [AiryDisc{2}(λ, NA), IsotropicGaussian{2}(λ, NA), BornWolf{2}(λ, NA, f)]
@testset "Blanket tests for 2D PSF: $(nameof(typeof(tf)))" for tf in PSF_2D
    @testset "`psf` method" begin
        let psf_array = psf(tf, Δxy, (11, 11))
            @test psf_array isa SpatialMatrix{<:Any,<:OAs.OffsetMatrix}
            @test all(==(Δxy), sampling(psf_array))
        end
        let psf_array = psf(tf, (Δx, Δy), (11, 11))
            @test psf_array isa SpatialMatrix{<:Any,<:OAs.OffsetMatrix}
            @test sampling(psf_array) == (Δx, Δy)
        end
    end
end

PSF_3D = [AiryDisc(λ, NA), IsotropicGaussian(λ, NA)] # , BornWolf(λ, NA, n)]
@testset "Blanket tests for 3D PSF: $(nameof(typeof(tf)))" for tf in PSF_3D
    @testset "`psf` method" begin
        let psf_array = psf(tf, Δxy, (11, 11, 11))
            @test psf_array isa SpatialArray{<:Any,3,<:OAs.OffsetArray}
            @test all(==(Δxy), sampling(psf_array))
        end
        let psf_array = psf(tf, (Δx, Δy, Δz), (11, 11, 11))
            @test psf_array isa SpatialArray{<:Any,3,<:OAs.OffsetArray}
            @test sampling(psf_array) == (Δx, Δy, Δz)
        end
    end
end
