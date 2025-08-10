A = ones(1024, 1024)

@testset "MeasuredPSF" begin
    # FIX: @test_throws DomainError PSFArray(ones(3, 3, 3), 32u"nm", (4, 1, 1))
    # FIX: @test PSFArray(ones(3, 3), 32u"nm") isa MeasuredPSF{<:Real,2}
    # FIX: @test MeasuredPSF(ones(3, 3, 3), 32u"nm") isa MeasuredPSF{<:Real,3}
end

@testset "PSFArray" begin end

@testset "ModelPSF" begin
    tf = BornWolf(488u"nm", 1.4, 1.7)

    ## Method Availability
    @test intensity(tf, 250u"nm", 200u"nm") isa Number
    @test intensity(tf, 250u"nm") isa Number
end

@testset "Sampled PSF" begin
    s_img = SpatialArray(A, 32u"nm")
    tf = BornWolf(488u"nm", 1.4, 1.7)

    ## Method Availability - Construction
    @test psf(tf, 64u"nm", (512, 512)) isa AbstractMatrix
    @test psf(tf, (64u"nm", 32u"nm"), (512, 512)) isa AbstractMatrix

    ## Methods - Array generation
    @test psf(tf, 60u"nm", (512, 512)) isa OAs.OffsetMatrix
    @test_throws MethodError psf(tf, s_img)

    ## Non-methods - Array generation
    @test_throws MethodError psf(tf, 60u"nm", (512.1, 512.4))
end
