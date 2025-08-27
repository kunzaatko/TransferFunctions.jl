using TransferFunctions: TransferFunctions as TF

# FIX: When there is an way to specify the interfacing type in the `InterfaceFunctions` package, I should test for these
# again <27-08-25> 
C_4D = SpatialArray(ones(40, 40), 20u"nm")

struct TF_1 <: TF.TransferFunction{2} end
tf_1 = TF_1()
@test_skip(@test_throws ["does not implement", r"transfer(.*::TransferFunction, .*::SpatialMatrix.*)"] transfer(C_4D, tf_1))
@test_skip(@test_throws ["does not implement", r"restore(.*::TransferFunction, .*::SpatialMatrix.*)"] restore(C_4D, tf_1))

struct LTF_1 <: TF.LinearShiftInvariantTransferFunction{2} end
ltf_1 = LTF_1()
@test_skip(@test_throws ["does not implement", r"conv(.*::LinearShiftInvariantTransferFunction, .*::SpatialMatrix.*)"] conv(C_4D, ltf_1))
@test_skip(@test_throws ["does not implement", r"deconv(.*::LinearShiftInvariantTransferFunction, .*::SpatialMatrix.*)"] deconv(C_4D, ltf_1))

struct PSF_1 <: TF.PointSpreadFunction{2} end
psf_1 = PSF_1()
@test_throws ["UnimplementedInterface{PointSpreadFunction}", "intensity(psf::PointSpreadFunction{2}, x::Length, y::Length)"] response(psf_1, 10u"nm", 10u"nm") # two argument intensity must be implemented

struct PSF_2 <: TF.PSFModel{2} end
TransferFunctions.symmetry(::PSF_2) = TF.ZAxisRadialSymmetry()
psf_2 = PSF_2()
@test_throws ["UnimplementedInterface{PointSpreadFunction}", "intensity(psf::PointSpreadFunction{2}, r::Length)"] response(psf_2, 10u"nm", 10u"nm") # single argument intensity must be implemented

struct OTF_1 <: TF.OpticalTransferFunction{2} end
otf_1 = OTF_1()
@test_throws ["does not implement", r"attenuation(.*::OpticalTransferFunction, .*::Frequency, .*::Frequency)"] attenuation(otf_1, 10u"nm^-1", 10u"nm^-1")

struct OTF_2 <: TF.RadialOTF{2} end
otf_2 = OTF_2()
@test_throws ["does not implement", r"attenuation(*::RadialOTF, .*::Frequency)"] attenuation(otf_2, 10u"nm^-1")
