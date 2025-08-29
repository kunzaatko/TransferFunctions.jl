using TransferFunctions.FFT
using FFTW

const real_arrays = [rand(10), rand(11), rand(10, 10), rand(11, 11), rand(10, 10, 10), rand(11, 11, 11)]

@testset "RFFT for size $(size(A))" for A in real_arrays
    @test FFT.fft(A) isa FFT.RFFTOut
    @test size(FFT.fft(A)) == size(A)
    @test fft(A) ≈ FFT.fft(A)
    @test collect(Base.broadcastable(FFT.fft(A))) == collect(FFT.fft(A))
end

@testset "Broadcasting" begin
    @testset "Broadcasting RFFT" for A in real_arrays
        rout = FFT.fft(A)
        @test rout .* 1 isa FFT.RFFTOut
        @test rout .* rand(size(A)...) isa FFT.FFTOut
    end
    @testset "Broadcasting Combination" for A in real_arrays
        rout = FFT.fft(A)
        CA = rand(ComplexF64, size(A)...)
        cout = FFT.fft(CA)
        @test rout .* cout isa FFT.FFTOut
        @test FFT.ifft(rout .* cout) ≈ ifft(fft(A) .* fft(CA))
    end
end
