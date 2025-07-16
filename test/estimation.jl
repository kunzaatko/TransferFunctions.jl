using TransferFunctions.Estimation

@test (bead(100u"nm", 30.5u"nm", intensity=0.75) .<= 0.75) |> all
@test_throws AssertionError bead(100u"nm", 30.5u"nm"; position=(-1.5, 0.7))
@test_throws AssertionError bead(100u"nm", 30.5u"nm"; position=(1.5, 0.7))

@test (bead(100u"nm", 30.5u"nm"; position=(0.5, 0.5)) .== reverse(bead(100u"nm", 30.5u"nm"; position=(-0.5, -0.5)))) |> all
@test (2bead(100u"nm", 30.5u"nm"; intensity=0.5) .== bead(100u"nm", 30.5u"nm")) |> all
