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

PSF_types = [AiryDisc(λ, NA), IsotropicGaussian(λ, NA), BornWolf(λ, NA, n)]
PSF_3D = [AiryDisc(λ, NA), IsotropicGaussian(λ, NA)]

s_img = SpatialArray(A, Δxy)
@testset "`psf` method: $(nameof(typeof(tf)))" for tf in PSF_types
    let psf_array = psf(tf, Δxy, (11, 11))
        @test psf_array isa SpatialMatrix{<:Any,<:OAs.OffsetMatrix}
        @test all(==(Δxy), TF.sampling(psf_array))
    end
    let psf_array = psf(tf, (Δx, Δy), (11, 11))
        @test psf_array isa SpatialMatrix{<:Any,<:OAs.OffsetMatrix}
        @test TF.sampling(psf_array) == (Δx, Δy)
    end

    if tf in PSF_3D
        @testset "`psf` 3D array" begin
            let psf_array = psf(tf, Δxy, (11, 11, 11))
                @test psf_array isa SpatialArray{<:Any,3,<:OAs.OffsetArray}
                @test all(==(Δxy), TF.sampling(psf_array))
            end
            let psf_array = psf(tf, (Δx, Δy, Δz), (11, 11, 11))
                @test psf_array isa SpatialArray{<:Any,3,<:OAs.OffsetArray}
                @test TF.sampling(psf_array) == (Δx, Δy, Δz)
            end
        end
    end
end
