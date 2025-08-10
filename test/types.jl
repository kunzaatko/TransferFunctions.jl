using TransferFunctions
using TransferFunctions: TransferFunctions as TF

@testset "PixelSize" begin
    @test (61u"nm", 61u"nm") isa TF.PixelSize{2}
    @test (61, 61) isa TF.Coordinate{2}
    @test (61.5, 61.0) isa TF.Coordinate{2}
    @test (61, 61) isa TF.Size{2}
end
