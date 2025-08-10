using FFTViews, TransferFunctions

A = ones(1024, 1024)

@testset "OTFArray" begin
    # FIX: @test_throws DomainError OTFArray(ones(3, 3, 3), 32u"nm", (4, 1, 1))
    # FIX: @test OTFArray(ones(3, 3), 32u"nm") isa OTFArray{<:Real}
    # FIX: @test OTFArray(ones(3, 3, 3), 32u"nm") isa OTFArray{<:Real,3}

    # FIX: tf = CircularPupilOTF(488u"nm", 1.4, 1.0, 0.3)
    # FIX: otf_array = otf(tf, 64u"nm", (512, 512))
    # FIX: tf_1 = OTFArray(otf_array, 64u"nm", (2, 3))

    # FIX: @test otf_array == otf(tf_1, 64u"nm", (512, 512))
    # FIX: @test argmax(FFTView(otf(tf_1, (512, 512)))) == CartesianIndex(1, 2)
    # FIX: @test argmax(FFTView(otf(tf_1, (511, 511)))) == CartesianIndex(1, 2)
    # FIX: @test otf(s_tf_1, (512, 512)) ≈ FourierTools.shift(otf(s_tf, (512, 512)), (1, 2))
    # FIX: @test otf(s_tf_1, (511, 511)) ≈ FourierTools.shift(otf(s_tf, (511, 511)), (1, 2))

    # TODO: When there is a model that returns Complex, i.e. the model with a phase shift <26-08-24> 
    # FIX: @test_skip otf(s_tf_1, (512, 512)) ≈ FourierTools.shift(otf(ComplexF32, s_tf, (512, 512)), (1, 2))
    # FIX: @test_skip otf(s_tf_1, (511, 511)) ≈ FourierTools.shift(otf(ComplexF32, s_tf, (511, 511)), (1, 2))

    # FIX: s_tf_2 = SampledOTF(tf, 64u"nm", (1.5, -2.5))
    # FIX: @test otf(s_tf_2, (512, 512)) ≈ real(FourierTools.shift(otf(ComplexF32, s_tf, (512, 512)), (0.5, -3.5)))
    # FIX: @test otf(s_tf_2, (511, 511)) ≈ real(FourierTools.shift(otf(ComplexF32, s_tf, (511, 511)), (0.5, -3.5)))

    # FIX: Probably a numerical error of FFT in the FourierTools package, but not sure <26-08-24> 
    # FIX: @test_broken otf(s_tf_2, (512, 512)) ≈ real(FourierTools.shift(otf(s_tf, (512, 512)), (0.5, -3.5)))
    # FIX: @test_broken otf(s_tf_2, (511, 511)) ≈ real(FourierTools.shift(otf(s_tf, (511, 511)), (0.5, -3.5)))
end

@testset "OTF models" begin
    @testset "$tf" for tf in [CircularPupilOTF(488u"nm", 1.4, 1.0, 0.3)]
        @testset "`OpticalTransferFunction` interface" begin
            @test attenuation(tf, 1 // 250u"nm", 1 // 200u"nm") isa AbstractFloat
            @test attenuation(tf, 250.0u"nm^-1", 200.0u"nm^-1") isa AbstractFloat

            if tf isa TF.RadialOTF
                @testset "`symmetry(otf)=Radial` interface" begin
                    @test attenuation(tf, 1 // 250u"nm") isa Number

                    c = 1 / 200u"nm"
                    @test attenuation(tf, c, 0u"nm^-1") == attenuation(tf, c) ≈
                          attenuation(tf, c / sqrt(10), 3c / sqrt(10)) ≈ attenuation(tf, c / sqrt(2), c / sqrt(2))

                    @test cutoff(tf) isa Frequency
                    @test cutoff(tf, 0.15) isa Frequency

                    ρ_max = cutoff(tf)
                    @test TF.insupport(tf, ρ_max / sqrt(2) / 2) == true
                    @test TF.insupport(tf, ρ_max / 2sqrt(2), ρ_max / 2sqrt(2)) == true
                    @test TF.insupport(tf, ρ_max / sqrt(10), 3ρ_max / 2sqrt(10)) == true
                end
            end
        end
        @testset "`LinearTransferFunction` interface" begin
            A = SpatialArray(rand(Float64, 10, 10), 60u"nm")

            @test_broken TF.conv(tf, A) isa SpatialArray
            @test_broken TF.deconv(tf, TF.conv(tf, A)) isa SpatialArray
            # FIX: @test attenuation(Float32, tf, 1 // 250u"nm") isa Float32
            # FIX: @test attenuation(ComplexF32, tf, 1 // 250u"nm") isa ComplexF32
        end
    end
end

@testset "Sampled OTF" begin
    tf = CircularPupilOTF(488u"nm", 1.4, 1.0, 0.3)
    s_img = SpatialArray(A, 32u"nm")

    ## Methods - Array generation
    @test otf(tf, 60u"nm", (512, 512)) isa Matrix
    @test otf(tf, s_img) isa Matrix
    # FIX: @test otf(ComplexF32, s_tf, (512, 512)) isa Matrix{ComplexF32}

    # NOTE: zeroth frequency is the greatest <26-08-24> 
    @test argmax(FFTView(otf(tf, 60u"nm", (512, 512)))) == CartesianIndex(0, 0)
    @test FFTView(otf(tf, 60u"nm", (512, 512)))[0, 0] == 1

    # FIX: using TransferFunctions: support
    # NOTE: The support calculation should match the generation of the array <26-08-24>
    # FIX: @test support(s_tf, img) isa BitArray
    # FIX: @test (otf(s_tf, (512, 512)) .> 0) == support(s_tf, (512, 512))
    # FIX: @test (otf(s_tf, (512, 512)) .>= 0.15) == support(s_tf, (512, 512), a=0.15)
    # FIX: @test (otf(s_tf_1, (512, 512)) .> 0) == support(s_tf_1, (512, 512))
    # FIX: @test (otf(s_tf_1, (512, 512)) .>= 0.15) == support(s_tf_1, (512, 512), a=0.15)
    # FIX: @test (otf(s_tf_2, (512, 512)) .> 0) == support(s_tf_2, (512, 512))
    # FIX: @test (otf(s_tf_2, (512, 512)) .>= 0.15) == support(s_tf_2, (512, 512), a=0.15)
    # FIX: @test count(support(s_tf, (512, 512), a=0.15)) < count(support(s_tf, (512, 512)))
    # FIX: @test all(otf(s_tf, (512, 512))[support(s_tf, (512, 512)).==false] .== 0)

    # FIX: using TransferFunctions: overlap
    # FIX: @test overlap(s_tf, s_tf, img) isa BitArray
    # FIX: @test overlap(s_tf, s_tf, (512, 512)) == support(s_tf, (512, 512))
    # FIX: @test overlap(s_tf, s_tf, (512, 512); a_1=0.15) == support(s_tf, (512, 512), a=0.15)
    # FIX: @test count(overlap(s_tf_1, s_tf, (512, 512))) < count(support(s_tf, (512, 512)))

    ## Non-methods - Array generation
    @test_throws MethodError otf(tf, 64u"nm", (512.1, 512.4))
end
