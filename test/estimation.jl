using TransferFunctions: reflect
using TestImages
using OffsetArrays: no_offset_view

@testset "LeastSquares" begin
  @testset "theory" begin
    tf = AiryDisc{2}(488u"nm", 1.4)
    gt = SpatialMatrix(Float32.(TestImages.testimage("mandril_gray"))[begin:5:end, begin:5:end], 64u"nm")
    indices = (-5:5, -5:5)
    psf_matrix = psf(tf, 64u"nm", (length.(indices)))
    raw = TF.conv(gt, psf_matrix)

    F_A = filtering_matrix(reflect(gt), indices)
    B = raw[TF.inner_axes(raw, indices)...]
    B_from_A = reshape(reflect(F_A' * psf_matrix[:]), axes(B))
    @test isapprox(B_from_A, B, rtol=1e-3)

    LhS = F_A * no_offset_view(reflect(F_A')) * psf_matrix[:]
    RhS = F_A * B[:]
    @test isapprox(LhS, RhS, rtol=1e-3)

    psf_estim = (F_A * no_offset_view(reflect(F_A'))) \ RhS
    @test isapprox(psf_matrix[:], psf_estim[:], rtol=1e-1)
  end

  @testset "estimation" begin
    using Test, Distributions
    tf = AiryDisc{2}(488u"nm", 1.4)
    bs = SyntheticData.beads(300, 100u"nm", (64u"nm", 64u"nm"), (256, 256))
    gt = SyntheticData.groundtruth(bs)
    ls = TF.Estimation.LeastSquares((-5:5, -5:5))
    raw = conv(gt, psf(tf, 64u"nm", (length.(ls.indices))))
    psf_estim = TF.Estimation.estimate(ls, raw, gt)
    psf_true = psf(tf, 64u"nm", (length.(ls.indices)))
    @test parent(psf_estim) ≈ psf_true

    # raw_noisy_1 = similar(raw)
    # raw_noisy_1 .= raw .+ 0.01*randn(size(raw)) / sum(raw.^2)
    # raw_noisy_2 = similar(raw)
    # raw_noisy_2 .= raw .+ 0.05*randn(size(raw)) / sum(raw.^2)
    # raw_noisy_3 = similar(raw)
    # raw_noisy_3 .= raw .+ 0.1*randn(size(raw)) / sum(raw.^2)
    # psf_estim_noisy_1 = TF.Estimation.estimate(ls, raw_noisy_1, gt)
    # psf_estim_noisy_2 = TF.Estimation.estimate(ls, raw_noisy_2, gt)
    # psf_estim_noisy_3 = TF.Estimation.estimate(ls, raw_noisy_3, gt)
  end
end
