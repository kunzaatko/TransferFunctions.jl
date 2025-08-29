"""
    scalebar(scale; kwargs...)

Add a scale-bar to the plot.

Supports all `lines()` and `text()` attributes, forwarding them to the respective plot calls.
The `position` attribute defines the position of the scale-bar in relative `Axis` coordinates.
The `axfrac` attribute defines the fraction of the axis width the scale-bar should span.
The multiple of `scale` will be chosen automatically (from `multiples`) so that the scale-bar length is closest to
`axfrac`.

Typically, `scale` is a `Unitful` quantity that defines the size of one plot unit.
For example, `scalebar!(1u"mm")` means that the plot units are millimetres.

```jldoctest setup=:(using Makie, GLMakie)
julia> scalebar(84u"nm");

julia> scalebar(84u"nm"; format=q -> q >= 1000u"nm" ? TF.scalebarformat(uconvert(u"μm",q)) : TF.scalebarformat(q))
```
"""
function scalebar end
function scalebar! end

"""
    scalebarformat(qty; sigdigits=3)
Function used to format the quantity length of the scale-bar to a string.

```jldoctest
julia> TransferFunctions.scalebarformat(700.4u"nm")
"700 nm"

julia> TransferFunctions.scalebarformat(70.4u"nm")
"70.4 nm"
```
"""
function scalebarformat end

export scalebar, scalebar!
