using TransferFunctions.Estimation

pixel_grid = 3 # PERF: Smaller grid for decreasing the test time

@test (bead(100u"nm", 30.5u"nm"; intensity=0.75, pixel_grid) .<= 0.75) |> all
@test_throws AssertionError bead(100u"nm", 30.5u"nm"; position=(-1.5, 0.7), pixel_grid)
@test_throws AssertionError bead(100u"nm", 30.5u"nm"; position=(1.5, 0.7), pixel_grid)

@test (bead(100u"nm", 30.5u"nm"; position=(0.5, 0.5), pixel_grid) .== reverse(bead(100u"nm", 30.5u"nm"; position=(-0.5, -0.5), pixel_grid))) |> all
@test (2bead(100u"nm", 30.5u"nm"; intensity=0.5, pixel_grid) .== bead(100u"nm", 30.5u"nm"; pixel_grid)) |> all
