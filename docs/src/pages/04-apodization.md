```@meta
CurrentModule = TransferFunctions.Apodization
```

# Apodization

```@example apodization
using TransferFunctions: Apodization as Apo
```

## Apodization functions

```@docs; canonical=false
Blackman
ExactBlackman
BlackmanHarris
BlackmanNuttall
```
```@example apodization
blackman = Apo.Blackman()
exact_blackman = Apo.ExactBlackman()
blackman_harris = Apo.BlackmanHarris()
blackman_nuttall = Apo.BlackmanNuttall()
```
```@makie apodization
f,a,_ = lines(blackman;  label="Blackman")
lines!(a, exact_blackman; linestyle=:dash, label="ExactBlackman")
lines!(a, blackman_harris; label="BlackmanHarris")
lines!(a, blackman_nuttall;linestyle=:dash, label="BlackmanNuttall")
axislegend(a; position=:rt)
f
```
```@docs; canonical=false
Connes
```
```@example apodization
connes = Apo.Connes()
```
```@makie apodization
f,a,_ = lines(connes; label="Connes")
axislegend(a; position=:rt)
f
```
```@docs; canonical=false
Cosine
```
```@example apodization
cosine = Apo.Cosine()
```
```@makie apodization
f,a,_ = lines(cosine; label="Cosine")
axislegend(a; position=:rt)
f
```
```@docs; canonical=false
Gaussian
```
```@example apodization
gaussian = Apo.Gaussian(0.2)
```
```@makie apodization
f,a,_ = lines(gaussian; label="Gaussian")
axislegend(a; position=:rt)
f
```
```@docs; canonical=false
Hamming
```
```@example apodization
hamming = Apo.Hamming()
```
```@makie apodization
f,a,_ = lines(hamming; label="Hamming")
axislegend(a; position=:rt)
f
```
```@docs; canonical=false
Welch
```
```@example apodization
welch = Apo.Welch()
```
```@makie apodization
f,a,_ = lines(welch; label="Welch")
axislegend(a; position=:rt)
f
```
```@docs; canonical=false
PowerCosine
```
```@docs; canonical=false
Triangular
```
```@example apodization
triangular = Apo.Triangular()
```
```@makie apodization
f,a,_ = lines(triangular; label="Triangular")
axislegend(a; position=:rt)
f
```
```@docs; canonical=false
Nuttall
```
```@example apodization
nuttall = Apo.Nuttall()
```
```@makie apodization
f,a,_ = lines(nuttall; label="Nuttall")
axislegend(a; position=:rt)
f
```
```@docs; canonical=false
SineSum
```
```@docs; canonical=false
FlatTop
```
```@example apodization 
flat_top = Apo.FlatTop()
```
```@makie apodization
f,a,_ = lines(flat_top; label="FlatTop")
axislegend(a; position=:rt)
f
```
```@docs; canonical=false
Hann
```
```@example apodization
hann = Apo.Hann()
```
```@makie apodization 
f,a,_ = lines(hann; label="Hann")
axislegend(a; position=:rt)
f
```
