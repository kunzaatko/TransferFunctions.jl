using TransferFunctions: TransferFunctions as TF

λ = 488u"nm"
NA = 1.4
n = 1.5
f = 1.3

PSF_radially_symmetric = [AiryDisc, IsotropicGaussian] #, BornWolf]

PSF_2D = [AiryDisc{2}(λ, NA), IsotropicGaussian{2}(λ, NA), BornWolf{2}(λ, NA, f)]
@testset "Blanket tests 2D PSF model: $(nameof(typeof(tf)))" for tf in PSF_2D
  @test tf isa TF.PointSpreadFunction{2}
  @test tf isa TF.TransferFunction{2}
  @test response(tf, 300u"nm", 200u"nm") isa Number
  @test isconcretetype(eltype(response.(tf, range(-400u"nm", 400u"nm", 11), range(-400u"nm", 400u"nm", 11))))
  @test_throws MethodError response(tf, 300u"nm")
  @test_throws MethodError response(tf, 300u"nm", 200u"nm", 200u"nm")

  if any(typeof(tf) isa PSF_radial for PSF_radial in PSF_radially_symmetric)
    @testset "Radial PSF" begin
      @test TF.intensity(tf, 250u"nm") isa Number # radius method existence
      @test response(tf, 300u"nm", 200u"nm") == response(tf, 200u"nm", 300u"nm") # isotropic property
      @test allequal(TF.FWHM(tf)) # same FWHM
    end
  end
end

PSF_3D = [AiryDisc(λ, NA), IsotropicGaussian(λ, NA)] # , BornWolf(λ, NA, n)]
@testset "Blanket tests 3D PSF model: $(nameof(typeof(tf)))" for tf in PSF_3D
  @test tf isa TF.PointSpreadFunction{3}
  @test tf isa TF.TransferFunction{3}
  @test response(tf, 300u"nm", 200u"nm", 100u"nm") isa Number
  @test isconcretetype(eltype(response.(tf, range(-400u"nm", 400u"nm", 11), range(-400u"nm", 400u"nm", 11), range(-100u"nm", 100u"nm", 11))))
  @test_throws MethodError response(tf, 300u"nm", 200u"nm", 100u"nm", 100u"nm")
  @test_throws MethodError response(tf, 300u"nm", 200u"nm")

  if any(typeof(tf) isa PSF_radial for PSF_radial in PSF_radially_symmetric)
    @testset "Radial PSF" begin
      @test TF.intensity(tf, 250u"nm") isa Number # radius method existence
      @test response(tf, 300u"nm", 200u"nm", 100u"nm") == response(tf, 200u"nm", 300u"nm", 100u"nm") # isotropic property
      @test allequal(TF.FWHM(tf)[1:2]) # same FWHM
    end
  end
end

@testset "Airy" begin
  @testset "Constructor Invariants" begin
    @test_throws ArgumentError AiryDisc{1}(λ, NA, n)
    @test_throws DomainError AiryDisc(0u"nm", NA)
    @test_throws DomainError AiryDisc(λ, 0.0)
    @test_throws DomainError AiryDisc(λ, NA, 0.0)
  end

  let airy2d = AiryDisc{2}(λ, NA)
    @test_throws MethodError TF.intensity(airy2d, 0u"nm", 0u"nm", 0u"nm")
    @test TF.intensity(airy2d, 0u"nm") isa Number
    @test response(airy2d, 0u"nm", 0u"nm") isa Number
    @test_throws MethodError response(airy2d, 0u"nm", 0u"nm", 0u"nm")
  end

  let tf = AiryDisc(λ, NA)
    @test tf isa AiryDisc{3}
    FWHM_lateral = TF.C_AiryDisc_lateral * λ / NA
    @test all(≈(FWHM_lateral), TF.FWHM(tf)[1:2])
  end
end

@testset "IsotropicGaussian" begin
  @testset "Constructor Invariants" begin
    @test_throws ArgumentError IsotropicGaussian{1}(λ, NA)
    @test_throws DomainError IsotropicGaussian(0u"nm", NA)
    @test_throws DomainError IsotropicGaussian(λ, 0.0)
    @test_throws DomainError IsotropicGaussian(λ, NA; n=0.0)
    @test_throws DomainError IsotropicGaussian(λ, NA; C_lateral=0.0)
    @test_throws DomainError IsotropicGaussian(λ, NA; C_axial=0.0)
  end

  let gauss2d = IsotropicGaussian{2}(λ, NA)
    @test_throws MethodError TF.intensity(gauss2d, 0u"nm", 0u"nm", 0u"nm")
    @test TF.intensity(gauss2d, 0u"nm") isa Number
    @test response(gauss2d, 0u"nm", 0u"nm") isa Number
    @test_throws MethodError response(gauss2d, 0u"nm", 0u"nm", 0u"nm")
  end

  let tf = IsotropicGaussian(λ, NA)
    @test tf isa IsotropicGaussian{3}
    airy = AiryDisc(λ, NA)
    @test_broken all(TF.FWHM(airy) .≈ TF.FWHM(tf))
  end
end

@testset "Augmentations" begin
  @testset "ScaledPSF" begin
    @testset "Construtor Invariants" begin
      @test_throws ["ArgumentError", "scale"] ScaledPSF(AiryDisc{2}, λ, NA; xscale=0.0)
      @test_throws ["MethodError", "keyword arguments"] ScaledPSF(AiryDisc{2}, λ, NA; zscale=1.0)
    end
  end

  @testset "RotatedPSF" begin
    @testset "Constructor" begin
      @test TF.RotatedPSF(AiryDisc{3}, λ, NA; α=0.1, β=0.2, γ=0.3) isa TF.RotatedPSF{3}
      @test TF.RotatedPSF(AiryDisc{2}, λ, NA; θ=0.1) isa TF.RotatedPSF{2}
    end

    using Rotations
    x, y = 300u"nm", 200u"nm"

    rot_scaled_airy3d = TF.RotatedPSF(ScaledPSF{3}, AiryDisc{3}, λ, NA; α=π / 3, xscale=0.8)
    scaled_airy3d = ScaledPSF(AiryDisc{3}, λ, NA; xscale=0.8)
    rot_scaled_airy3d_rotation = RotXYZ(π / 3, 0, 0)
    @test rot_scaled_airy3d isa TF.RotatedPSF{3}
    @test response(rot_scaled_airy3d, x, y) == response(scaled_airy3d, (rot_scaled_airy3d_rotation * [x, y, 0u"nm"])...)

    rot_scaled_airy2d = TF.RotatedPSF(ScaledPSF{2}, AiryDisc{2}, λ, NA; θ=π / 3, yscale=0.3)
    scaled_airy2d = ScaledPSF(AiryDisc{2}, λ, NA; yscale=0.3)
    rot_scaled_airy2d_rotation = RotMatrix{2}(π / 3)
    @test rot_scaled_airy2d isa TF.RotatedPSF{2}
    @test response(rot_scaled_airy2d, x, y) == response(scaled_airy2d, (rot_scaled_airy2d_rotation * [x, y])...)
  end

  @testset "Composability" begin
    @test TF.RotatedPSF(ScaledPSF{2}, AiryDisc{2}, λ, NA; θ=π / 3, xscale=0.8) isa TF.RotatedPSF{2}
    @test TF.RotatedPSF(ScaledPSF{3}, AiryDisc{3}, λ, NA; α=π / 3, xscale=0.8) isa TF.RotatedPSF{3}
  end
end
