using TransferFunctions.Apodization
using TransferFunctions: Apodization as Apo
using Distributions
using JET

@testset "ApodizationFuntction Constructors" begin
    @testset "Constructor $(nameof(apo))" for apo in [
        Apo.Triangular,
        Apo.Welch,
        Apo.Connes,
        Apo.Cosine,
        Apo.Hann,
        Apo.Hamming,
        Apo.Nuttall,
        Apo.BlackmanNuttall,
        Apo.BlackmanHarris,
        Apo.FlatTop,
        Apo.Blackman,
        Apo.Blackman,
        Apo.ExactBlackman
    ]
        @test apo() isa Apo.ApodizationFunction
        @test apo{Float64}() isa Apo.ApodizationFunction
    end

    @test Apo.Gaussian(rand(Uniform(0.01, 0.49))) isa Apo.ApodizationFunction

    @testset "Construction exceptions Blackman" begin
        @test_throws ArgumentError Apo.Blackman(0.5)
        @test_throws ArgumentError Apo.Blackman(0.6)
        @test_throws ArgumentError Apo.Blackman(0)
    end

    @testset "Construction exceptions Gaussian" begin
        @test_throws ArgumentError Apo.Gaussian(0.5)
        @test_throws ArgumentError Apo.Gaussian(0)
    end
end



@testset "methods $(typeof(apo))" for apo in [
    Apo.Hamming(),
    Apo.Hann(),
    Apo.Welch(),
    Apo.Connes(),
    Apo.Cosine(),
    Apo.Nuttall(),
    Apo.BlackmanNuttall(),
    Apo.BlackmanHarris(),
    Apo.FlatTop(),
    Apo.ExactBlackman(),
    Apo.Blackman(0.4),
    Apo.Gaussian(0.4),
]
    @test Apo.apodization(apo, 0) ≈ 1

    r = rand()
    @test Apo.apodization(apo, r) ≈ Apo.apodization(apo, -r)

    if typeof(apo) ∈ (Apo.Hann, Apo.Triangular, Apo.Cosine, Apo.Connes, Apo.Welch)
        @test Apo.apodization(apo, 1) ≈ Apo.apodization(apo, -1) == 0
    end
end

@testset "equavalence $(equivs)" for equivs in [
    (Apo.Hann(), Apo.SineSum(0.5, 0.5)),
    (Apo.PowerCosine(0), Apo.SineSum(1.0)),
    (Apo.PowerCosine(2), Apo.SineSum(0.5, 0.5)),
    (Apo.PowerCosine(4), Apo.SineSum(0.375, 0.5, 0.125)),
    (Apo.PowerCosine(6), Apo.SineSum(0.3125, 0.46875, 0.1875, 0.03125)),
]
    @test Apo.apodization.(equivs[1], -1:0.01:1) ≈ Apo.apodization.(equivs[2], -1:0.01:1)
end

if VERSION <= v"1.12"
    @testset "JET: `apodization` $apo" for apo in [
        Apo.Triangular,
        Apo.Welch,
        Apo.Connes,
        Apo.Cosine,
        Apo.Hann,
        Apo.Hamming,
        Apo.Nuttall,
        Apo.BlackmanNuttall,
        Apo.BlackmanHarris,
        Apo.FlatTop,
        Apo.Blackman,
        Apo.Blackman,
        Apo.ExactBlackman,
    ]
        for T1 in [Float32, Float64, ComplexF32]
            for T2 in [Float32, Float64]
                a = apo{T1}()
                @test_opt Apo.apodization(a, T2(0.5))
            end
        end
    end
end
