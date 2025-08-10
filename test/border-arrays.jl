using TransferFunctions
using TransferFunctions: TransferFunctions as TF
using OffsetArrays: OffsetArrays as OAs
using JET

A = reshape(float.(1:9), (3, 3))
A3 = reshape(float.(1:27), (3, 3, 3))
@testset "BorderArray constructor $border" for border in (
    (TF.Reflect, TF.Reflect{Float64}, TF.Reflect(), :reflect)...,
    (TF.Symmetric, TF.Symmetric{Float64}, TF.Symmetric(), :symmetric)...,
    (TF.Circular, TF.Circular{Float64}, TF.Circular(), :circular)...,
    (TF.Replicate, TF.Replicate{Float64}, TF.Replicate(), :replicate)...,
    (TF.Fill, TF.Fill{Float64}, TF.Fill(1.0), TF.Fill(Float32(1.0)), :fill)...
)
    for padding in (((2, 2), (2, 2)), (2, 2), 2)
        @test TF.BorderArray(A, border, padding) isa TF.BorderArray
        @test size(TF.BorderArray(A, border, padding)) == (7, 7)
        @test axes(TF.BorderArray(A, border, padding)) == (-1:5, -1:5)
    end
    for padding in (((2, 2), (2, 2), (2, 2)), (2, 2, 2), 2)
        @test TF.BorderArray(A3, border, padding) isa TF.BorderArray
        @test size(TF.BorderArray(A3, border, padding)) == (7, 7, 7)
        @test axes(TF.BorderArray(A3, border, padding)) == (-1:5, -1:5, -1:5)
    end
end

@testset "BorderArray equivalence $border" for (border, template) in (
    (
        TF.Reflect,
        [
            5.0 2.0 2.0 5.0 8.0 8.0 5.0;
            4.0 1.0 1.0 4.0 7.0 7.0 4.0;
            4.0 1.0 1.0 4.0 7.0 7.0 4.0;
            5.0 2.0 2.0 5.0 8.0 8.0 5.0;
            6.0 3.0 3.0 6.0 9.0 9.0 6.0;
            6.0 3.0 3.0 6.0 9.0 9.0 6.0;
            5.0 2.0 2.0 5.0 8.0 8.0 5.0
        ]
    ),
    (
        TF.Symmetric,
        [
            9.0 6.0 3.0 6.0 9.0 6.0 3.0;
            8.0 5.0 2.0 5.0 8.0 5.0 2.0;
            7.0 4.0 1.0 4.0 7.0 4.0 1.0;
            8.0 5.0 2.0 5.0 8.0 5.0 2.0;
            9.0 6.0 3.0 6.0 9.0 6.0 3.0;
            8.0 5.0 2.0 5.0 8.0 5.0 2.0;
            7.0 4.0 1.0 4.0 7.0 4.0 1.0
        ]
    ),
    (
        TF.Circular,
        [
            5.0 8.0 2.0 5.0 8.0 2.0 5.0;
            6.0 9.0 3.0 6.0 9.0 3.0 6.0;
            4.0 7.0 1.0 4.0 7.0 1.0 4.0;
            5.0 8.0 2.0 5.0 8.0 2.0 5.0;
            6.0 9.0 3.0 6.0 9.0 3.0 6.0;
            4.0 7.0 1.0 4.0 7.0 1.0 4.0;
            5.0 8.0 2.0 5.0 8.0 2.0 5.0
        ]
    ),
    (
        TF.Replicate,
        [
            1.0 1.0 1.0 4.0 7.0 7.0 7.0;
            1.0 1.0 1.0 4.0 7.0 7.0 7.0;
            1.0 1.0 1.0 4.0 7.0 7.0 7.0;
            2.0 2.0 2.0 5.0 8.0 8.0 8.0;
            3.0 3.0 3.0 6.0 9.0 9.0 9.0;
            3.0 3.0 3.0 6.0 9.0 9.0 9.0;
            3.0 3.0 3.0 6.0 9.0 9.0 9.0
        ]
    ),
    (
        TF.Fill,
        [
            0.0 0.0 0.0 0.0 0.0 0.0 0.0;
            0.0 0.0 0.0 0.0 0.0 0.0 0.0;
            0.0 0.0 1.0 4.0 7.0 0.0 0.0;
            0.0 0.0 2.0 5.0 8.0 0.0 0.0;
            0.0 0.0 3.0 6.0 9.0 0.0 0.0;
            0.0 0.0 0.0 0.0 0.0 0.0 0.0;
            0.0 0.0 0.0 0.0 0.0 0.0 0.0
        ]
    ),
)
    @test OAs.no_offset_view(TF.BorderArray(A, border, 2)) == template
end

@testset "Exceptions" begin
    @testset "BoundsError" begin
        ba = TF.BorderArray(A, TF.Circular, 2)
        I = (-2, -2)
        @test getindex(ba.border, ba.parent, I...) == 1.0 # Mapped index is in the parent array
        @test_throws BoundsError ba[I...] # but not in the `BorderArray`
    end
    @testset "InvalidBorderExtent" begin
        @test_throws TF.InvalidBorderExtent TF.BorderArray(A, TF.Circular, ((3, 4), (1, 1)))
        @test TF.BorderArray(A, TF.Circular, 3) isa TF.BorderArray
        @test_throws TF.InvalidBorderExtent TF.BorderArray(A, TF.Reflect, ((3, 4), (1, 1)))
        @test TF.BorderArray(A, TF.Reflect, 3) isa TF.BorderArray
        @test_throws TF.InvalidBorderExtent TF.BorderArray(A, TF.Symmetric, 3)
        @test TF.BorderArray(A, TF.Symmetric, 2) isa TF.BorderArray
        @test TF.BorderArray(A, TF.Fill, 100) isa TF.BorderArray
        @test TF.BorderArray(A, TF.Replicate, 100) isa TF.BorderArray
    end
end

@testset "JET: `getindex` $border" for border in (
    TF.Reflect, TF.Symmetric, TF.Circular, TF.Replicate, TF.Fill
)
    ba = TF.BorderArray(A, border, 2)

    @test_opt getindex(ba, 0, 0) # border
    @test_opt getindex(ba, 3, 3) # inside
    @test_opt getindex(ba, 0:2, 0:2) # slices
end
