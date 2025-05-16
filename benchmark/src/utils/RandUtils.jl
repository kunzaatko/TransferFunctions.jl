module RandUtils

using Random
using StableRNGs

samerand(args...) = rand(StableRNG(1), args...)

samerandstring(n) = randstring(StableRNG(1), n)

randvec(T, n) = samerand(T, n)
randvec(n) = samerand(n)

randmat(T, n) = samerand(T, n, n)
randmat(n) = samerand(n, n)

export samerand, samerandstring, randvec, randmat, StableRNGs

end
