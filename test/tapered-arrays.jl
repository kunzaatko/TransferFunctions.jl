using TransferFunctions
using TransferFunctions: TransferFunctions as TF
using TransferFunctions: Apodization as Apo
using TestImages
using JET

A = ones(50, 50)
A3 = ones(50, 50, 50)

@testset "TaperedArray Construction" begin
    @test TF.TaperedArray(A, Apo.Hann(), (2, 2)) isa TF.TaperedArray
end

@testset "TaperedArray functionality" begin
    let tap = TF.TaperedArray(A, Apo.Hann(), (2, 2))
        @test tap[1, :] == tap[end, :] == tap[:, 1] == tap[:, end] == zeros(length(axes(tap, 1)))
    end
end

@testset "taperedges" begin
    @test taperedges(A, 2) isa TF.TaperedArray
    @test taperedges(A, (2, 2)) isa TF.TaperedArray
    @test taperedges(A, 2, :reflect) isa TF.TaperedArray
    @test let tap = taperedges(A, 2, :reflect)
        parent(tap) isa TF.BorderArray
    end
end

@testset "ImageCore" begin
    img_gray = TestImages.testimage("mandril_gray")
    tap_gray = TF.taperedges(img_gray, 50)

    @test tap_gray[1, 1] == zero(eltype(tap_gray))
    @test allequal(typeof, TF.taperedges(img_gray, 50))

    img_rgb = TestImages.testimage("mandril_color")
    tap_rgb = TF.taperedges(img_rgb, 50)

    @test tap_rgb[1, 1] == zero(eltype(tap_rgb))
    @test allequal(typeof, TF.taperedges(img_rgb, 50))
end

if VERSION <= v"1.12"
    @testset "JET: `getindex`" begin
        ba = TF.TaperedArray(A, Apo.Hann(), 2)
        @test_opt ba[1, 1]
        @test_opt ba[1, :]
    end
end
