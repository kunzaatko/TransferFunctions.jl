using TransferFunctions.Apodization
using TransferFunctions: Apodization as Apo

@testset "Types" begin
    @test Apo.Triangular() isa Apo.ApodizationFunction
    @test Apo.Welch() isa Apo.ApodizationFunction
    @test Apo.Connes() isa Apo.ApodizationFunction
    @test_throws ArgumentError Apo.PowerCosine{1 + 1im}()
    @test Apo.Cosine() isa Apo.PowerCosine{1}
    @test Apo.Hann() isa Apo.PowerCosine{2}
    @test Apo.Hamming() isa Apo.SineSum{2}
    @test Apo.Nuttall() isa Apo.SineSum{4}
    @test Apo.BlackmanNuttall() isa Apo.SineSum{4}
    @test Apo.BlackmanHarris() isa Apo.SineSum{4}
    @test Apo.FlatTop() isa Apo.SineSum{5}
    @test_throws ArgumentError Apo.Blackman(0.5)
    @test_throws ArgumentError Apo.Blackman(0.6)
    @test_throws ArgumentError Apo.Blackman(0)
    @test Apo.Blackman() == Apo.Blackman(0.16)
    @test Apo.Blackman() isa Apo.ApodizationFunction
    @test Apo.ExactBlackman() isa Apo.Blackman
    @test Apo.Gaussian(rand(Uniform(0.01, 0.49))) isa Apo.ApodizationFunction
    @test_throws ArgumentError Apo.Gaussian(0.6)
end


# FIX: Most of the apodization / window functions that are defined, should be 0 at the edge (for me the edge
# is in 1 for implementation reasons) and all should be 1 in the center. That is at 0 in my implementation.
# <02-09-24> 
@testset "methods $(typeof(apo))" for apo in [
    Apo.Hamming(), Apo.Hann(), Apo.Welch(), Apo.Connes(), Apo.Cosine(),
    Apo.Nuttall(), Apo.BlackmanNuttall(), Apo.BlackmanHarris(), Apo.FlatTop(), Apo.ExactBlackman(),
    Apo.Blackman(0.4), Apo.Gaussian(0.4),
]
    @test Apo.apodization(apo, 0) == 1

    r = rand()
    @test Apo.apodization(apo, r) ≈ Apo.apodization(apo, -r)

    if typeof(apo) ∈ (Apo.Hann, Apo.Triangular, Apo.Cosine, Apo.Connes, Apo.Welch)
        @test Apo.apodization(apo, 1) == Apo.apodization(apo, -1) == 0
    end
end

@testset "equavalence $(equivs)" for equivs in [
    (Apo.Hann(), Apo.SineSum((0.5, 0.5))),
    (Apo.PowerCosine{0}(), Apo.SineSum((1.0,))),
    (Apo.PowerCosine{2}(), Apo.SineSum((0.5, 0.5))),
    (Apo.PowerCosine{4}(), Apo.SineSum((0.375, 0.5, 0.125))),
    (Apo.PowerCosine{6}(), Apo.SineSum((0.3125, 0.46875, 0.1875, 0.03125))),
]
    @test Apo.apodization.(equivs[1], -1:0.01:1) ≈ Apo.apodization.(equivs[2], -1:0.01:1)
end

@testset "internals" begin
    O = Ones(30, 30, 30)
    @test Apo.specify_border(O, IF.Fill(5.0), ((5, 5, 5), (2, 2, 2)), (1, 2, 3)) == IF.Fill(5.0, (5, 5, 5), (2, 2, 2))
    @test Apo.specify_border(O, IF.Fill(5.0), ((5, 5), (2, 2)), (1, 3)) == IF.Fill(5.0, (5, 0, 5), (2, 0, 2))
    @test Apo.specify_border(O, IF.Pad(:replicate), ((5,), (2,)), (1,)) == IF.Pad(:replicate, (5, 0, 0), (2, 0, 0))
    @test Apo.specify_border(O, IF.Inner(), ((5,), (2,)), (1,)) == IF.Inner((5, 0, 0), (2, 0, 0))
    @test_throws ArgumentError Apo.specify_border(O, IF.NoPad(), ((5,), (2,)), (1,)) == IF.Pad(:replicate, (5, 0, 0), (2, 0, 0))
end

@testset "methods" begin
    apo = Apo.Cosine()
    @test taperedges(
              apo, ones(100, 100, 9), ((10, 10), (10, 10)); dims=(1, 2) # All the supplied arguments 
          ) == taperedges(
              apo, ones(100, 100, 9), ((10, 10), (10, 10)) # Infer dims
          ) == taperedges(
              apo, ones(100, 100, 9), 10; dims=(1, 2) # Supply single unified width
          ) == taperedges(
              apo, ones(100, 100, 9), (10, 10); dims=(1, 2) # Supply single width for each dimension
          ) == taperedges(
              apo, ones(100, 100, 9), (10, 10) # Infer dims
          ) == taperedges(
              ones(100, 100, 9), (10, 10) # default `apo`
          ) == taperedges(
              ones(100, 100, 9), (10, 10), "replicate"; dims=(1, 2) # default border 
          ) == taperedges(
              ones(100, 100, 9), (10, 10), IF.Pad{0}(:replicate, (), ())  # instantiate border
          ) == taperedges(
              ones(100, 100, 9), (10, 10), IF.Pad{3}(:replicate, (10, 10, 0), (10, 10, 0))  # Correct border size
          )

    @test_throws DimensionMismatch taperedges(apo, ones(100, 100), 10; dims=(1, 2, 3))
    @test_throws DimensionMismatch taperedges(apo, ones(100, 100), 10; dims=(1, 3))

    # `dims`
    @test taperedges(ones(30, 30), 10; dims=:) isa AbstractArray
    @test taperedges(ones(30, 30), 10; dims=2) isa AbstractArray
    @test taperedges(ones(30, 30), 10; dims=Dims((1, 2))) isa AbstractArray

    @test size(taperedges(ones(30, 30), (10, 20))) == (50, 70)

    A_tap = taperedges(Apo.Cosine(), ones(30, 30), 10)

    # RESEARCH: Should this in fact be 0 at the other edge as well? This is done to make the signal periodic
    # so if we taper one edge to 0 and the other edge to the 0 length - 1, we will already have a periodic
    # signal, correct? <10-09-24> 
    @test all(all.([
        (A_tap[-9, :] .== 0),
        (A_tap[40, :] .== 0),
        (A_tap[:, -9] .== 0),
        (A_tap[:, 40] .== 0),
        (A_tap[1:30, 1:30] .== 1),
        (A_tap[-9:0, -9:0] .!= 1),
        (A_tap[31:40, -9:0] .!= 1),
        (A_tap[-9:0, 31:40] .!= 1),
        (A_tap[31:40, 31:40] .!= 1)
    ]))
end
