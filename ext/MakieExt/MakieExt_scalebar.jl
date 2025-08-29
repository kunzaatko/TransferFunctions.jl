using Unitful, ComputePipeline

function scalebarformat(qty::Quantity; sigdigits=3)
    mul = round(ustrip(qty); sigdigits)
    mul = isinteger(mul) ? Int(mul) : mul
    return "$(mul) $(unit(qty))"
end

# NOTE: Gracefully stolen and adapted from `MakieExtra.jl` all credit goes to @aplavin. Any issues are on me (@kunzaatko) <31-07-25>
@recipe Scalebar (scale,) begin
    begin 
        lines_attributes = Makie.documented_attributes(Lines)
        filter!(a -> a.first ∉ [:color, :linewidth], lines_attributes.d)
        lines_attributes
    end...
    begin
        text_attributes = Makie.documented_attributes(Makie.Text)
        filter!(a -> a.first ∉ ([:position] ∪ keys(Makie.documented_attributes(Lines).d)), text_attributes.d)
        text_attributes
    end...
    
    "Colour of the scale text under the `Scalebar`"
    textcolor =  :white
    "Colour of the `Scalebar` line"
    linecolor =  :yellow
    "Width of the line in during the maximum length of the `Scalebar`"
    linewidth = @inherit linewidth
    "Whether to update the width of the bar of the `Scalebar` with zooming"
    updatewidth = true
    "Fraction of the axis that is filled during the maximum length of the `Scalebar`"
    axfrac = 0.2
    "Position of the bar in the axis"
    position = Point2(0.85, 0.08)
    "Multiples of the quantity that are used for the length of the bar of `Scalebar`"
    multiples = [x*p for p in Real[[10.0^p for p in -50:-1]; [1, 10, 100, 1000, 10000]; [10.0^p for p in 5:50]] for x in [1, 2, 5]]
    "Function to format the quantity representing the length of the scale-bar into a string `(qty::Quantity) -> String`. (By default uses [`scalebarformat`](@ref))"
    format = scalebarformat
    cycle = nothing
end


Makie.data_limits(::Scalebar) = Rect3f(Point3f(NaN), Vec3f(NaN))
Makie.boundingbox(::Scalebar, space::Symbol=:data) = Rect3f(Point3f(NaN), Vec3f(NaN))

# TODO: Add conversion logic for the position attribute to change to a `Vec2` type <22-07-25>

function Makie.plot!(p::Scalebar)
    scene = Makie.parent_scene(p)
    @assert Makie.transform_func(scene)[1] == identity

    add_input!(p.attributes, :projview, Makie.projview_to_2d_limits(p))

    register_computation!(p.attributes, [:projview], [:hlims]) do inputs, changed, cached
        (first.(extrema(inputs.projview)),)
    end

    map!(p.attributes, [:scale, :position, :multiples, :hlims, :axfrac, :linewidth, :updatewidth, :format], [:mul, :barlinepoints, :barlinewidth]) do scale, position, multiples, hlims, axfrac, linewidth, updatewidth, format
        data_units = ustrip(scale) 
        mul  = argmin(multiples) do m
            abs(1 / data_units * m - axfrac * (hlims[2] - hlims[1]))
        end
        length_data = 1 / data_units * mul
        length_ax = length_data / (hlims[2] - hlims[1])
        barlinepoints = [position - Vec2(length_ax/2, 0), position + Vec2(length_ax/2, 0)]

        if updatewidth
            barlinewidth = linewidth * (length_ax/axfrac)
        else
            barlinewidth = linewidth
        end
        (Float64(mul), barlinepoints, barlinewidth)
    end

    map!(p.attributes, [:scale, :mul, :format], :scaletext) do scale, mul, format
        format(mul * unit(scale))
    end

    l = lines!(p, p.attributes, p.barlinepoints, xautolimits=false, yautolimits=false, space=:relative, color=p.linecolor, linewidth=p.barlinewidth)
    t = text!(p, p.attributes, p.scaletext, xautolimits=false, yautolimits=false, align = (:center, :top), space=:relative, color=p.textcolor)
    
    return p
end
