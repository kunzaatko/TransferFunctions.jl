using TransferFunctions

psf1 = BornWolf(488u"nm", 1, 1.7)
psf2 = BornWolf(488u"nm", 1.0, 1.7)

@test psf1 == psf2
@test hash(psf1) == hash(psf2)
@test isequal(psf1, psf2)
