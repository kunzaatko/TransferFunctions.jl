include("models/gibson_lanni.jl")
include("models/born_wolf.jl")

# FIX: Add support for N-dims <30-11-23> 
@doc """
  Intensity point spread function i.e. the intensity ratio and phase shift of the sample intensity density
 """ psf

# FIX: @repl blocks not working <15-09-23> 
@doc """
    psf(tf::TransferFunction, x::Length, y::Length)

Sample the PSF of the transfer function model/data at the point `(x,y)`

```jldoctest
julia> tf = BornWolf(488u"nm", 1.4, 1.7)
BornWolf{Float64}(488.0 nm, 1.4, 1.7)
```

```@repl
tf = BornWolf(488u"nm", 1.4, 1.7) # hide
psf(tf, 0u"nm", 5u"nm")
psf.(tf, 0u"nm", -400u"nm":100u"nm":400u"nm")
```
"""
@inline @traitfn function psf(OT::Type{<:Number}, tf::TF, x::Length, y::Length) where {TF <: ModelPSF; RadiallySymmetric{TF}}
    # FIX: The type should be passed to the method of on the model which for it to be able to optimize its run based on
    # it, and instead there should be a catch-all method that implements it if its not implemented in the implementation
    # <26-08-24> 
    OT(psf(tf, hypot(x, y)))
end

@inline @traitfn function psf(tf::TF, x::Length, y::Length) where {TF <: ModelPSF; RadiallySymmetric{TF}}
    psf(preferred_type(TF), tf, x, y)
end

function resolution_limit(tf::TF) where {TF<:TransferFunction}
    # TODO: Is this correct?
    all(hasfield.(TF, [:NA, :λ, :nᵢ])) ? tf.λ / (2 * tf.NA * tf.nᵢ) : nothing
end
