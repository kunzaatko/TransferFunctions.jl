using TransferFunctions: TransferFunctions as TF

C_4D = TF.SpatialArray(ones(40, 40), 20u"nm")

struct TF_1 <: TF.TransferFunction end
tf_1 = TF_1()
@test_throws ["does not implement", r"transfer(.*::TransferFunction, .*::SpatialMatrix.*)"] transfer(tf_1, C_4D)
@test_throws ["does not implement", r"restore(.*::TransferFunction, .*::SpatialMatrix.*)"] restore(tf_1, C_4D)

struct LTF_1 <: TF.LinearTransferFunction end
ltf_1 = LTF_1()
@test_throws ["does not implement", r"conv(.*::LinearTransferFunction, .*::SpatialMatrix.*)"] TF.conv(ltf_1, C_4D)
@test_throws ["does not implement", r"deconv(.*::LinearTransferFunction, .*::SpatialMatrix.*)"] TF.deconv(ltf_1, C_4D)

struct PSF_1 <: TF.PointSpreadFunction end
psf_1 = PSF_1()
@test_throws ["does not implement", r"intensity(.*::PointSpreadFunction, .*::Length, .*::Length)"] intensity(psf_1, 10u"nm", 10u"nm")
@test_throws MethodError intensity(psf_1, 10u"nm")

struct PSF_2 <: TF.RadialPSF end
psf_2 = PSF_2()
@test_throws ["does not implement", r"intensity(.*::RadialPSF, .*::Length)"] intensity(psf_2, 10u"nm")

struct OTF_1 <: TF.OpticalTransferFunction end
otf_1 = OTF_1()
@test_throws ["does not implement", r"attenuation(.*::OpticalTransferFunction, .*::Frequency, .*::Frequency)"] attenuation(otf_1, 10u"nm^-1", 10u"nm^-1")

struct OTF_2 <: TF.RadialOTF end
otf_2 = OTF_2()
@test_throws ["does not implement", r"attenuation(*::RadialOTF, .*::Frequency)"] attenuation(otf_2, 10u"nm^-1")
