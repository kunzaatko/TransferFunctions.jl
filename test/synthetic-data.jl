using TransferFunctions
using TransferFunctions.SyntheticData
using Test, Distributions

@testset "beads" begin
  @testset "construction" begin
    @test let bs = SyntheticData.beads(100, 100u"nm", (64u"nm", 64u"nm"), (256, 256))
      bs isa Base.get_extension(TransferFunctions, :TestExt).Beads
    end
    @test_broken let bs = SyntheticData.beads(100, 100u"nm", 64u"nm", (256, 256))
      bs isa Base.get_extension(TransferFunctions, :TestExt).Beads
    end
    @test_throws r"Failed to generate.*decrease `N`.* or the `spacing`.*" SyntheticData.beads(200, 100u"nm", (64u"nm", 64u"nm"), (100, 100); spacing =1000u"nm")
  end
  @testset "groundtruth" begin
    bs = SyntheticData.beads(100, 100u"nm", (64u"nm", 64u"nm"), (256, 256))
    @test SyntheticData.groundtruth(bs) isa SpatialArray
  end
end
