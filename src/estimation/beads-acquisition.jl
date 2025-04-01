using Unitful
using Unitful: Length, Quantity, 𝐋
using OffsetArrays: OffsetMatrix, OffsetArray
using Statistics
using TransferFunctions: PixelSize
const PerLength = Quantity{<:Any,inv(𝐋)}

"""
    bead([T=Float64], d, Δxy, [α=0u"nm^-1"],)
    bead(d, α, Δxy)
Generate a model of a bead with diameter `d` and pixel-size `Δxy`.

# Arguments
- `α = 1u"nm^-1`: evanescent wave attenuation constant in the [TIR-FM microscopy modulation](@cite 2025). For a normal acquisition without total internal reflection 0u"nm^-1" is setting.
- `pixel_grid_length::Int = 10`: length of each pixel in the grid
- `peak_intensity = 1.0`: peak intensity value
- `subpixel_shift = (0.0, 0.0)`

# Examples
```jldoctest; setup = :(using TransferFunctions: Estimation)
julia> Estimation.bead(100u"nm", 30.5u"nm");

julia> Estimation.bead(95u"nm", (30.5u"nm", 25u"nm"));

julia> Estimation.bead(Float32, 0.1u"μm", 30.5u"nm");

julia> Estimation.bead(0.1u"μm", 30.5u"nm"; α=0.01u"nm^-1")
5×5 OffsetArray(::Matrix{Float64}, -2:2, -2:2) with eltype Float64 with indices -2:2×-2:2:
 0.0        0.0       0.0705075  0.0       0.0
 0.0        0.606779  0.892914   0.606779  0.0
 0.0705075  0.892914  1.0        0.892914  0.0705075
 0.0        0.606779  0.892914   0.606779  0.0
 0.0        0.0       0.0705075  0.0       0.0

julia> Estimation.bead(0.1u"μm", 50.5u"nm"; subpixel_shift=(-0.3,-0.2))
3×3 OffsetArray(::Matrix{Float64}, -1:1, -1:1) with eltype Float64 with indices -1:1×-1:1:
 0.333333  0.737374  0.0808081
 0.59596   1.0       0.20202
 0.020202  0.141414  0.0

julia> Estimation.bead(0.1u"μm", 50.5u"nm"; peak_intensity=0.5)
3×3 OffsetArray(::Matrix{Float64}, -1:1, -1:1) with eltype Float64 with indices -1:1×-1:1:
 0.03  0.23  0.03
 0.23  0.5   0.23
 0.03  0.23  0.03
```
"""
function bead(
    T::Type{<:Real},
    d::Length,
    Δxy::PixelSize{2};
    α::PerLength=0u"nm^-1",
    pixel_grid_length=10, peak_intensity=one(T), subpixel_shift::Tuple{Real,Real}=(0.0, 0.0)
)::OffsetMatrix{T}
    @assert all(abs.(subpixel_shift) .< one(eltype(subpixel_shift))) "|subpixel_shift| must be < 1"
    subpx_pad = map(subpixel_shift) do x
        (x < zero(x) ? x : 0, x > zero(x) ? x : 0)  # padding to apply to the buffer to fit the sub-pixel shift
    end
    buf_axes_px = map(Tuple((-x, x) for x in (d ./ (2 .* Δxy)) .- 1 / 2), subpx_pad) do ax_range, pad # center pixel size => -½
        round.(Int, (ax_range .+ pad), RoundFromZero) # range of the axes in px including the sub-pixel padding
    end

    grid = Matrix{T}(undef, (((b - a + 1) * pixel_grid_length) for (a, b) in buf_axes_px)...)
    grid_subpixel_shift = subpixel_shift .* pixel_grid_length
    grid_center = ((1 / 2 .- first.(buf_axes_px)) .* pixel_grid_length) .+ grid_subpixel_shift

    r_grid = map(CartesianIndices(grid)) do xy
        (Tuple(xy) .- grid_center .- 1 / 2) ./ pixel_grid_length .* Δxy |> splat(hypot)
    end # radii from the center
    supp_grid = r_grid .< (d / 2) # support of the bead
    z_grid = similar(r_grid, Union{Missing,Length}) # z-axes offset from the focal plane
    z_grid[supp_grid.==0] .= missing
    z_grid[supp_grid] .= (d / 2) .- sqrt.((d / 2)^2 .- r_grid[supp_grid] .^ 2)
    intensity_grid = map(z_grid) do z
        ismissing(z) ? zero(T) : exp(-α * z)
    end

    buf = Matrix{T}(undef, ((b - a + 1) for (a, b) in buf_axes_px)...)
    map!(buf, CartesianIndices(buf)) do ind
        indsx, indsy = map(Tuple(ind)) do i
            ((i-1)*pixel_grid_length+1):(i*pixel_grid_length)
        end
        mean(intensity_grid[indsx, indsy])
    end
    buf = OffsetArray(buf, (first.(buf_axes_px) .- 1)...)
    buf .*= peak_intensity / buf[0, 0]
    return buf
end
# NOTE: step 1 add default type
bead(d, Δxy; vargs...) = bead(Float64, d, Δxy; vargs...)
# NOTE: step 2 assume same axes sampling
bead(T::Type{<:Real}, d, Δxy::Length; vargs...) = bead(T, d, (Δxy, Δxy); vargs...)

struct Beads
    positions::AbstractVector{Tuple{Real,Real}}
    Δxy::Tuple{Length,Length} # pixel size
    wh::Tuple{Int,Int}
    diameter::Length
    α_evanescent::PerLength
end

export bead
