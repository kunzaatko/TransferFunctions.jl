using TransferFunctions: TransferFunctions as TF
using OffsetArrays: OffsetArrays as OAs
using JET

V = collect(1:9)
A = reshape(float.(1:9), (3, 3))
A3 = reshape(float.(1:27), (3, 3, 3))
O_V = OAs.OffsetArray(V, -3:5)
O_A = OAs.OffsetArray(A, 0:2, -1:1)
O_A3 = OAs.OffsetArray(A3, 0:2, -1:1, -2:0)

@testset "ReflectedArray constructor" begin
    @test TF.ReflectedArray(V) isa TF.ReflectedVector
    @test TF.ReflectedArray(A) isa TF.ReflectedMatrix
    @test TF.ReflectedArray(A3) isa TF.ReflectedArray{<:Any,3}
    @test TF.ReflectedArray(O_V) isa TF.ReflectedVector
    @test TF.ReflectedArray(O_A) isa TF.ReflectedMatrix
    @test TF.ReflectedArray(O_A3) isa TF.ReflectedArray{<:Any,3}
    @testset "reflect $(typeof(a))" for a in (V, A, A3, O_V, O_A, O_A3)
        @test TF.reflect(a) isa TF.ReflectedArray
    end
end

@testset "AbstractArray interface" begin
    let ra = TF.ReflectedArray(A)
        @test axes(ra) == (-3:-1, -3:-1)
        @test size(ra) == (3, 3)
        @test ra[begin] == 9
        @test ra[end] == 1
        @test_throws BoundsError ra[begin-1]
        @test_throws BoundsError ra[end+1]
    end
    let ra = TF.ReflectedArray(A3)
        @test axes(ra) == (-3:-1, -3:-1, -3:-1)
        @test size(ra) == (3, 3, 3)
        @test ra[begin] == 27
        @test ra[end] == 1
        @test_throws BoundsError ra[begin-1]
        @test_throws BoundsError ra[end+1]
    end
    let ra = TF.ReflectedArray(O_A)
        @test axes(ra) == (-2:0, -1:1)
        @test size(ra) == (3, 3)
        @test ra[begin] == 9
        @test ra[end] == 1
        @test_throws BoundsError ra[begin-1]
        @test_throws BoundsError ra[end+1]
    end
end

if VERSION <= v"1.12"
    @testset "JET: `getindex` $a" for a in (V, A, A3, O_V, O_A, O_A3)
        ra = TF.ReflectedArray(a)
        I = rand.(axes(ra))
        @test_opt getindex(ra, I...)
    end
end
