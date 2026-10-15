using TransferFunctions: TransferFunctions as TF
using OffsetArrays: OffsetArray as OA, OffsetArrays as OAs
using Combinatorics, JET

@testset "trivia $fn" for fn in (TF.corr, TF.conv)
    @test fn(ones(2, 2), ones(2, 2)) == fill(4, (2, 2))
    @test fn(fill(1 + 1im, (2, 2)), ones(2, 2)) == fill(4 + 4im, (2, 2))
    @test fn(fill(1 + 1im, (2, 2)), ones(1, 1)) == fill(1 + 1im, (2, 2))
    @test fn(fill(1 + 1im, (2, 2)), im * ones(1, 1)) == fill(1 - 1im, (2, 2))
    @test fn(ones(2, 2, 2), ones(1, 1, 1)) == ones(2, 2, 2)
    @test fn(ones(2, 2, 2) + im * ones(2, 2, 2), ones(1, 1, 1)) == fill(1 + 1im, (2, 2, 2))
    for T in (Float64, ComplexF64)
        let A = rand(T, 10, 10)
            @test fn(A, OA(ones(1, 1), -1, -1)) ≈ A
        end
    end
    let A = rand(ComplexF64, 10, 10)
        K = OA(ones(1, 1) .+ im .* ones(1, 1), -1, -1)
        @test fn(A, K) ≈ A .- im .* A
    end
end

@testset "filtering_matrix equivalence" begin
    A = reshape(1:16, (4, 4))
    A_fm = filtering_matrix(A, (-1:1, -1:1), :periodic)
    fA_fm = reshape(A_fm' * ones(9), size(A))
    fA_corr = TF.corr(A, OA(ones(3, 3), -1:1, -1:1))
    @test fA_fm == fA_corr
end

@testset "ImageFiltering equivalence" begin
    using ImageFiltering
    # NOTE: `imfilter`'s default FIR algorithm multiplies the kernel unconjugated, whereas `TF.corr`
    # is FFT-based and conjugates the kernel spectrum. To compare the two, either pass an explicitly
    # conjugated kernel to the FIR algorithm, or select the FFT algorithm (which conjugates) with the
    # original kernel. Both algorithms are selected explicitly to pin down this convention.
    for T in (Float64, ComplexF64)
        img = rand(T, 21, 21)
        kernel = OAs.centered(rand(T, 3, 3))
        @test isapprox(imfilter(img, conj.(kernel), "circular", Algorithm.FIR()), TF.corr(img, kernel), rtol=1e-3)
        @test isapprox(imfilter(img, kernel, "circular", Algorithm.FFT()), TF.corr(img, kernel), rtol=1e-3)
    end
end

@testset "corr periodic oracle" begin
    # Direct periodic-sum oracle for correlation:
    #     corr(A, K)[I] = Σ_J A[I + J] * conj(K[J])
    # with periodic (wrapping) indexing over A's axes. The kernel K is centered (offset axes) and
    # the output shares the (possibly offset) axes of A. `TF.corr` is FFT-based, so a small relative
    # tolerance absorbs its floating-point round-off.
    corr_oracle(A, K) = begin
        T = promote_type(eltype(A), eltype(K))
        out = similar(A, T)
        lo = map(a -> first(a), axes(A))
        n = size(A)
        for I in CartesianIndices(A)
            s = zero(T)
            for J in CartesianIndices(K)
                idx = map((l, i, j, nn) -> l + mod(i + j - l, nn), lo, Tuple(I), Tuple(J), n)
                s += A[idx...] * conj(K[J])
            end
            out[I] = s
        end
        out
    end

    A_c = reshape(ComplexF64[1 + 2im, 3 - 1im, 2 + 4im, -1 + 2im, 5, 2 + 1im, -3 + 2im, 1 - 4im, 4 + 2im], (3, 3))
    A_r = reshape(Float64[1, 3, 2, -1, 5, 2, -3, 1, 4], (3, 3))
    K_c = OAs.centered(reshape(ComplexF64[1, 2 + 1im, 1 - 2im, 1im, 1 + 1im, 2, 1 + 3im, -1im, 1], (3, 3)))
    K_r = OAs.centered(reshape(Float64[1, 2, 1, 0, 1, 2, 1, 0, 1], (3, 3)))

    # complex-real and real-complex, dense 1-based input
    let A = A_c, K = K_r
        A_snap, K_snap = copy(A), copy(K)
        @test TF.corr(A, K) ≈ corr_oracle(A, K)
        @test A == A_snap && K == K_snap          # inputs unmutated
    end
    let A = A_r, K = K_c
        A_snap, K_snap = copy(A), copy(K)
        @test TF.corr(A, K) ≈ corr_oracle(A, K)
        @test A == A_snap && K == K_snap
    end

    # offset input → offset output axes
    let A = OAs.OffsetArray(copy(A_c), -1:1, 0:2), K = K_c
        A_snap, K_snap = copy(A), copy(K)
        @test axes(TF.corr(A, K)) == axes(A)
        @test TF.corr(A, K) ≈ corr_oracle(A, K)
        @test A == A_snap && K == K_snap
    end
end

if VERSION >= v"1.12-rc"
    @testset "JET: filter $(typeof(A)), $(typeof(B))" for (A, B) in map(Tuple, combinations((ones(3, 3), fill(1 + 1im, (3, 3)), OA(ones(3, 3), -1, -1), OA(fill(1 + 1im, (3, 3)), -1, -1)), 2))
        function fn_filter(@nospecialize f)
            f !== Base.materialize && f !== conj! && f !== Base.broadcasted
        end
        opt_broken = B isa OAs.OffsetMatrix{<:Real} && A isa Matrix{<:Real}
        @test_opt broken = opt_broken target_modules = (TransferFunctions,) function_filter = fn_filter TF.conv(A, B)
        @test_opt broken = opt_broken target_modules = (TransferFunctions,) function_filter = fn_filter TF.corr(A, B)
    end
end

@testset "ImageCore" begin
    using ImageCore, TestImages

    K_rand = OAs.centered(rand(3, 3))
    K_rand ./= sum(K_rand)

    @testset "Gray image" begin
        img_gray = TestImages.testimage("mandril_gray")

        @test TF.corr(img_gray, K_rand) isa AbstractMatrix{<:Gray}
        @test TF.conv(img_gray, K_rand) isa AbstractMatrix{<:Gray}
    end

    @testset "RGB image" begin
        img_rgb = TestImages.testimage("mandril_color")

        @test TF.corr(img_rgb, K_rand) isa AbstractMatrix{<:RGB}
        @test TF.conv(img_rgb, K_rand) isa AbstractMatrix{<:RGB}
    end
end
