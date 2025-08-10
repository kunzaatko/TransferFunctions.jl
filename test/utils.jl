using TransferFunctions: fillsize, roundupcenter, exactcenter, fftfreqs, posgrid, contained, interior, rounddowncenter, roundcenter, aroundorigin
using TransferFunctions: PixelSize, Coordinate, Frequency, Length, OriginAt
using OffsetArrays: OffsetArray as OA
using Base: CartesianIndex as CI

@testset "types" begin
    @testset "units" begin
        @test 1 / 32u"nm" isa Frequency
        Δkx = 1 / 61u"nm"
        @test 1 / Δkx isa Length
    end

    @test (31.5u"nm", 40u"nm", 50u"nm") isa PixelSize{}
    @test (31.5u"nm", 40u"nm") isa PixelSize{2}

    @test (3, 3, 3) isa Coordinate{3}
    @test (3, 3, 3) isa Coordinate{3,Int}
    @test (3.0, 3.0, 3.0) isa Coordinate{3,<:AbstractFloat}
    @test (2.4, 3, 3) isa Coordinate{3,Real} # NOTE: Must accept diverse types if supplied <26-08-24> 
end


@testset "utility functions" begin
    @test fillsize(31u"nm", 2) == (31u"nm", 31u"nm")
    @test_throws MethodError fillsize(31u"nm", Val(2)) == (31u"nm", 31u"nm") # NOTE: method with integer should always be used <05-05-25> 

    @test (posgrid((11, 11), (31u"nm", 31u"nm")) |> first |> first) isa Length
    @test posgrid((11, 11), 31u"nm") == posgrid((11, 11), (31u"nm", 31u"nm"))
    @test_throws MethodError posgrid((11, 11, 11), (31u"nm", 31u"nm"))

    @test (fftfreqs((11, 11), (31u"nm", 31u"nm")) |> first |> first) isa Frequency
    @test fftfreqs((11, 11), 31u"nm") == fftfreqs((11, 11), (31u"nm", 31u"nm"))
    @test_throws MethodError fftfreqs((11, 11, 11), (31u"nm", 31u"nm"))

    @test roundcenter(RoundFromZero, ones(3, 4, 2)) == CI(2, 3, 2)

    @test roundupcenter(ones(3, 4, 2)) == CI(2, 3, 2)
    @test roundupcenter(OA(ones(3, 3, 3), -2, -2, -2)) == CI(0, 0, 0)
    @test roundupcenter(OA(ones(4, 3, 3), -2, -2, -2)) == CI(1, 0, 0)

    @test rounddowncenter(ones(3, 4, 2)) == CI(2, 2, 1)
    @test rounddowncenter(OA(ones(3, 3, 3), -2, -2, -2)) == CI(0, 0, 0)
    @test rounddowncenter(OA(ones(4, 3, 3), -2, -2, -2)) == CI(0, 0, 0)

    @test exactcenter(ones(3, 4, 2)) == (2.0, 2.5, 1.5)
    @test exactcenter(OA(ones(3, 3, 3), -2, -2, -2)) == (0.0, 0.0, 0.0)
    @test exactcenter(OA(ones(4, 3, 3), -2, -2, -2)) == (0.5, 0.0, 0.0)

    @test contained(ones(3, 4, 2), (2, 3, 1))
    @test !contained(ones(3, 4, 2), (8, 3, 1))
    @test contained(OA(ones(3, 3, 3), -2, -2, -2), (-1, 1, 0))
    @test !contained(OA(ones(3, 3, 3), -2, -2, -2), (-2, 1, 0))

    # TODO: Test other types of axes that may occur in an array that I use (OffsetAxes) <05-05-25> 
    @test interior(1:9, -1:5) == 2:4
    @test interior(Base.OneTo(9), -1:5) == 2:4
    @test interior((0:3, -1:3, -3:1), (-1:1, -1:1, -1:1)) == (1:2, 0:2, -2:0)
    @test interior((Base.OneTo(3), -1:3, -3:1), (0:1, -1:1, -1:1)) == (1:2, 0:2, -2:0)

    @test aroundorigin((-3:3, -3:5), (2, 1)) == (-1:5, -2:6)
    @test aroundorigin(-3:4, 4) == 1:8
    @test aroundorigin(-3:4) == -3:4
    @test aroundorigin((3, 3, 3)) == (-1:1, -1:1, -1:1)

    @test OriginAt(CI(2, 2, 2))(ones(3, 3, 3)) == OA(ones(3, 3, 3), -2, -2, -2)
end
