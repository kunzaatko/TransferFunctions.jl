```@meta
CurrentModule = TransferFunctions.Apodization
CollapsedDocStrings = true
```

# Apodization

```@docs
Apodization
ApodizationFunction
```

```@example apodization
using TransferFunctions: Apodization as Apo
```

## Apodization functions

```@docs
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
nothing # hide
```

```@makie apodization; formats=:svg, basename="apodization_blackman"
f,a,_ = lines(blackman;  label="Blackman")
lines!(a, exact_blackman; linestyle=:dash, label="ExactBlackman")
lines!(a, blackman_harris; label="BlackmanHarris")
lines!(a, blackman_nuttall;linestyle=:dash, label="BlackmanNuttall")
axislegend(a; position=:rt)
f
```

```@docs
Connes
```

```@example apodization
connes = Apo.Connes{Float64}()
nothing # hide
```

```@makie apodization; formats=:svg, basename="apodization_connes"
f,a,_ = lines(connes; label="Connes")
axislegend(a; position=:rt)
f
```

```@docs
Cosine
```

```@example apodization
cosine = Apo.Cosine{Float64}()
nothing # hide
```

```@makie apodization; formats=:svg, basename="apodization_cosine"
f,a,_ = lines(cosine; label="Cosine")
axislegend(a; position=:rt)
f
```

```@docs
Gaussian
```

```@example apodization
gaussian = Apo.Gaussian(0.2)
nothing # hide
```

```@makie apodization; formats=:svg, basename="apodization_gaussian"
f,a,_ = lines(gaussian; label="Gaussian")
axislegend(a; position=:rt)
f
```

```@docs
Hamming
```

```@example apodization
hamming = Apo.Hamming{Float64}()
nothing # hide
```

```@makie apodization; formats=:svg, basename="apodization_hamming"
f,a,_ = lines(hamming; label="Hamming")
axislegend(a; position=:rt)
f
```

```@docs
Welch
```

```@example apodization
welch = Apo.Welch{Float64}()
nothing # hide
```

```@makie apodization; formats=:svg, basename="apodization_welch"
f,a,_ = lines(welch; label="Welch")
axislegend(a; position=:rt)
f
```

```@docs
PowerCosine
```

```@docs
Triangular
```

```@example apodization
triangular = Apo.Triangular{Float64}()
nothing # hide
```

```@makie apodization; formats=:svg, basename="apodization_triangular"
f,a,_ = lines(triangular; label="Triangular")
axislegend(a; position=:rt)
f
```

```@docs
Nuttall
```

```@example apodization
nuttall = Apo.Nuttall{Float64}()
nothing # hide
```

```@makie apodization; formats=:svg, basename="apodization_nuttall"
f,a,_ = lines(nuttall; label="Nuttall")
axislegend(a; position=:rt)
f
```

```@docs
SineSum
```

```@docs
FlatTop
```

```@example apodization 
flat_top = Apo.FlatTop{Float64}()
nothing # hide
```

```@makie apodization; formats=:svg, basename="apodization_flat_top"
f,a,_ = lines(flat_top; label="FlatTop")
axislegend(a; position=:rt)
f
```

```@docs
Hann
```

```@example apodization
hann = Apo.Hann{Float64}()
nothing # hide
```

```@makie apodization; formats=:svg, basename="apodization_hann"
f,a,_ = lines(hann; label="Hann")
axislegend(a; position=:rt)
f
```
