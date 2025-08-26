using TransferFunctions

@testset "beads" begin
  pixel_grid = 3 # PERF: Smaller grid for decreasing the test time

  @testset "bead" begin
    @test (SyntheticData.bead(100u"nm", 30.5u"nm"; intensity=0.75, pixel_grid) .<= 0.75) |> all
    @test_throws AssertionError SyntheticData.bead(100u"nm", 30.5u"nm"; position=(-1.5, 0.7), pixel_grid)
    @test_throws AssertionError SyntheticData.bead(100u"nm", 30.5u"nm"; position=(1.5, 0.7), pixel_grid)

    @test (SyntheticData.bead(100u"nm", 30.5u"nm"; position=(0.5, 0.5), pixel_grid) .== reverse(SyntheticData.bead(100u"nm", 30.5u"nm"; position=(-0.5, -0.5), pixel_grid))) |> all
    @test (2SyntheticData.bead(100u"nm", 30.5u"nm"; intensity=0.5, pixel_grid) .== SyntheticData.bead(100u"nm", 30.5u"nm"; pixel_grid)) |> all
  end

  @testset "construction" begin
    @test let bs = SyntheticData.beads(100, 100u"nm", (64u"nm", 64u"nm"), (256, 256))
      bs isa SyntheticData.Beads
    end
    @test let bs = SyntheticData.beads(100, 100u"nm", 64u"nm", (256, 256))
      bs isa SyntheticData.Beads
    end
    @test_throws r"Failed to generate.*decrease `N`.* or the `spacing`.*" SyntheticData.beads(200, 100u"nm", (64u"nm", 64u"nm"), (100, 100); spacing=1000u"nm")
  end
  @testset "groundtruth" begin
    bs = SyntheticData.beads(100, 100u"nm", (64u"nm", 64u"nm"), (256, 256))
    @test SyntheticData.groundtruth(bs) isa SpatialArray
  end
end
