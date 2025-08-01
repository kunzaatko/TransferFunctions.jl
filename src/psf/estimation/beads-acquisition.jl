using Unitful
using Unitful: Length, Quantity, 𝐋
using OffsetArrays: OffsetMatrix, OffsetArray
using Statistics
using TransferFunctions: PixelSize, SpatialArray
const PerLength = Quantity{<:Any,inv(𝐋)}

"""
    bead([T=Float64], d, Δ; <kwargs>)
Generate a model of a fluorescent microsphere (AKA calibration bead) with diameter `d` and pixel-size `Δ`.

# Arguments
- `α=0u"nm^-1`: evanescent wave attenuation constant in the [TIR-FM microscopy modulation](@cite 2025). An acquisition that is not TIR-FM is equivalent to setting `α=0u"nm^-1"`.
- `pixel_grid=10`: aliasing of the pixels is done by averaging values at a larger grid. `pixel_grid` sets the dimensions equivalent to each pixel in the grid
- `intesity=1.0`: peak intensity value, i.e. the theoretical value at the center of the bead
- `position=(0.0, 0.0)`: set the subpixel position of the peak (center) of the bead within the generated array 

# Examples
```jldoctest; setup = :(using TransferFunctions: Estimation)
julia> Estimation.bead(100u"nm", 30.5u"nm");

julia> Estimation.bead(95u"nm", (30.5u"nm", 25u"nm"));

julia> Estimation.bead(Float32, 0.1u"μm", 30.5u"nm");

julia> Estimation.bead(0.1u"μm", 30.5u"nm"; α=0.01u"nm^-1")
5×5 SpatialArray{Float64, 2, Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}, OffsetArrays.OffsetMatrix{Float64, Matrix{Float64}}} with indices -2:2×-2:2:
 0.0        0.0       0.0705075  0.0       0.0
 0.0        0.606779  0.892914   0.606779  0.0
 0.0705075  0.892914  1.0        0.892914  0.0705075
 0.0        0.606779  0.892914   0.606779  0.0
 0.0        0.0       0.0705075  0.0       0.0

julia> Estimation.bead(0.1u"μm", 50.5u"nm"; position=(-0.3,-0.2))
3×3 SpatialArray{Float64, 2, Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}, OffsetArrays.OffsetMatrix{Float64, Matrix{Float64}}} with indices -1:1×-1:1:
 0.333333  0.737374  0.0808081
 0.59596   1.0       0.20202
 0.020202  0.141414  0.0

julia> Estimation.bead(0.1u"μm", 50.5u"nm"; intensity=0.5)
3×3 SpatialArray{Float64, 2, Quantity{Float64, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}, OffsetArrays.OffsetMatrix{Float64, Matrix{Float64}}} with indices -1:1×-1:1:
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
    pixel_grid=10,
    intensity=one(T),
    position::Tuple{Real,Real}=(0.0, 0.0)
)
    @assert all(abs.(position) .< one(eltype(position))) "|position| must be < 1"
    subpx_pad = map(position) do x
        (x < zero(x) ? x : 0, x > zero(x) ? x : 0)  # padding to apply to the buffer to fit the sub-pixel shift
    end
    buf_axes_px = map(Tuple((-x, x) for x in (d ./ (2 .* Δxy)) .- 1 / 2), subpx_pad) do ax_range, pad # center pixel size => -½
        round.(Int, (ax_range .+ pad), RoundFromZero) # range of the axes in px including the sub-pixel padding
    end

    grid_size = (((b - a + 1) * pixel_grid) for (a, b) in buf_axes_px)
    @assert all(>(0), grid_size) LazyString("Invalid grid size: $(collect(grid_size)). The problem may be in the combination of `d` and `Δxy` arguments such that the bead does not span a single pixel in some dimension.")
    grid = Matrix{T}(undef, grid_size...)
    grid_subpixel_shift = position .* pixel_grid
    grid_center = ((1 / 2 .- first.(buf_axes_px)) .* pixel_grid) .+ grid_subpixel_shift

    r_grid = map(CartesianIndices(grid)) do xy
        (Tuple(xy) .- grid_center .- 1 / 2) ./ pixel_grid .* Δxy |> splat(hypot)
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
            ((i-1)*pixel_grid+1):(i*pixel_grid)
        end
        mean(intensity_grid[indsx, indsy])
    end
    buf = OffsetArray(buf, (first.(buf_axes_px) .- 1)...)
    buf .*= intensity / buf[0, 0]
    return SpatialArray(buf, Δxy)
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
