using TransferFunctions
using TransferFunctions: TransferFunctions as TF
using OffsetArrays: OffsetArrays as OAs
using OffsetArrays: OffsetArray as OA
using ImageFiltering: ImageFiltering as IF
using TensorOperations

K = OA(zeros(9, 9), -4:4, -4:4)
K[0, 0] = 0.5
K[[-1, -1, 1, 1], [-1, 1, -1, 1]] .= 0.5 / 4

A = rand(30, 30)

@testset "CirculantTensor constructor" begin

    @test circulant(A, (-5:5, -5:5)) isa TF.CirculantTensor
    @test circulant(A, (-5:5, -5:3)) isa TF.CirculantTensor
    @test circulant(A, K) isa TF.CirculantTensor

    local A_9_9 = reshape(1:(9*9), 9, 9)
    local Kinds = (-2:2, -2:2)
    @testset "Padded constructors - $bord" for (bord, out) in [
        (
            :replicate, [
                1 1 1 10 19;
                1 1 1 10 19;
                1 1 1 10 19;
                2 2 2 11 20;
                3 3 3 12 21
            ]
        ),
        (
            :symmetric, [
                21 12 3 12 21;
                20 11 2 11 20;
                19 10 1 10 19;
                20 11 2 11 20;
                21 12 3 12 21
            ]
        ),
        (
            :circular, [
                71 80 8 17 26;
                72 81 9 18 27;
                64 73 1 10 19;
                65 74 2 11 20;
                66 75 3 12 21
            ]
        ),
        (
            :reflect, [
                11 2 2 11 20;
                10 1 1 10 19;
                10 1 1 10 19;
                11 2 2 11 20;
                12 3 3 12 21
            ]
        ),
        (
            TF.Fill, [
                0 0 0 0 0;
                0 0 0 0 0;
                0 0 1 10 19;
                0 0 2 11 20;
                0 0 3 12 21
            ]
        )
    ]
        ct = circulant(A_9_9, Kinds, bord)
        @test ct.interior == axes(A_9_9)
        @test OAs.no_offset_view(ct[:, :, 1, 1]) == out
    end

    # Different eltypes
    @test eltype(circulant(ones(Int, 30, 30), K)) == Int

    C_4D = circulant(A, K)

    # Correct output indices
    @test ndims(circulant(A, K)) == 4
    @test axes(C_4D)[1:2] == axes(K)
    @test axes(C_4D)[3:4] == C_4D.interior

    contract(F, K) = dropdims(
        mapslices(F, dims=Dims(1:(ndims(F)÷2))) do S
            sum(S .* K)
        end; dims=Dims(1:(ndims(C_4D)÷2)))

    B = contract(C_4D, K)
    @test B isa AbstractMatrix
    @test size(B) == length.(C_4D.interior)
end
