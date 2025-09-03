using TransferFunctions.FFT
using FFTW

const arrays = [rand(10), rand(11), rand(10, 10), rand(11, 11), rand(10, 10, 10), rand(11, 11, 11)]

const real_arrays = arrays
const cmplx_arrays = map(x -> ComplexF64.(x), arrays)

@testset "RFFT for size $(size(A))" for (A, cA) in zip(real_arrays, cmplx_arrays)
    @test FFT.fft(A) isa FFT.RFFTOut
    @test size(FFT.fft(A)) == size(A)
    @test FFT.fft(A) ≈ FFT.fft(cA)
    @test fft(A) ≈ FFT.fft(A)
    @test collect(Base.broadcastable(FFT.fft(A))) == collect(FFT.fft(A))
end

@testset "similar" begin
    fA = FFT.fft(rand(10))
    @test typeof(similar(fA)) == typeof(fA)
end

@testset "Broadcasting" begin
    @testset "Broadcasting RFFTOut" for A in real_arrays
        rout = FFT.fft(A)
        @test rout .* 1 isa FFT.RFFTOut
        @test rout .* rand(size(A)...) isa FFT.FFTOut
        for T in (Float64, Float32, Int32, Int64, ComplexF64, ComplexF32)
            @test eltype(rout .* rand(T, size(A)...)) == promote_type(eltype(rout), T)
        end
    end
    @testset "Broadcasting FFTOut" for A in cmplx_arrays
        cout = FFT.fft(A)
        @test cout .* 1 isa FFT.FFTOut
        @test cout .* rand(size(A)...) isa FFT.FFTOut
        for T in (Float64, Float32, Int32, Int64, ComplexF64, ComplexF32)
            @test eltype(cout .* rand(T, size(A)...)) == promote_type(eltype(cout), T)
        end
    end
    @testset "Broadcasting RFFTOut with FFTOut" for A in real_arrays
        rout = FFT.fft(A)
        CA = rand(ComplexF64, size(A)...)
        cout = FFT.fft(CA)
        @test rout .* cout isa FFT.FFTOut
        @test FFT.ifft(rout .* cout) ≈ ifft(fft(A) .* fft(CA))
    end
    @testset "Broadcasting Operators: $(nameof(typeof(A)))" for (A, cA) in zip(real_arrays, cmplx_arrays)
        fA = FFT.fft(A)
        fcA = FFT.fft(cA)
        @test Broadcast.combine_styles(FFT.fft(A), FFT.fft(cA)) isa Broadcast.ArrayStyle{<:FFT.FFTOut}
        @test conj!(fA) ≈ conj!(fcA)
        @test abs2.(fA) .+ 1 == begin
            A_fA = abs2.(fA)
            A_fA .+ 1
        end
    end
end

@testset "RFFTOut getindex/setindex!" begin
    @testset "setindex! for size $(size(A))" for A in real_arrays
        rout = FFT.fft(A)
        original_parent = deepcopy(parent(rout))
        # Attempt to set an index within the stored part
        rout[ones(Int, ndims(A))...] = 100.0 + 100.0im
        @test parent(rout)[ones(Int, ndims(A))...] ≈ 100.0 + 100.0im
        # Ensure other elements are unchanged
        @test parent(rout)[2, ones(Int, ndims(A) - 1)...] ≈ original_parent[2, ones(Int, ndims(A) - 1)...]

        # Attempt to set an index outside the stored part (should not modify anything and error)
        # The current implementation of `setindex!` only checks bounds of parent(a), so this will error if I is out of bounds of parent(a)
        # This is the correct behaviour, as we only want to modify the stored data.
        if ndims(A) == 2
            first_dim_orig = FFT.firstdim(rout)
            if first_dim_orig > size(parent(rout), 1)
                # throws a bounds error, as `setindex!` is allowed to operate only within the bounds of the parent array
                @test_throws BoundsError rout[size(parent(rout), 1)+1, 1] = 1.0 + 1.0im
            end
        end
    end
end

@testset "ImageCore ext" begin
    using ImageCore, TestImages
    img = testimage("mandril_gray")
    @test FFT.fft(img) isa FFT.RFFTOut
    @test FFT.fft(img) == FFT.fft(float.(gray.(img)))
end
