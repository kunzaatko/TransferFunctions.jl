module TestExt
using TransferFunctions
using TransferFunctions: TransferFunctions as TF
using Distributions
using Unitful: Length, Quantity, 𝐋
using LinearAlgebra
import TransferFunctions.SyntheticData

const Size = Tuple{Int, Int}
const PerLength = Quantity{<:Any,inv(𝐋)}
const PixelSize = Tuple{Length, Length}

# TODO: Synthetic Data for N-dims <21-08-25> 
struct Beads{T, D} <: SyntheticData.SyntheticModel{2}
    positions::Vector{Tuple{T,T}}
    Δxy::Tuple{Length,Length} # pixel size
    wh::Size
    diameter::D
    α_evanescent::PerLength # FIX: Instead of this, I would like something like params to the beads that will be stored <21-08-25> 
end

# TODO: Randomized or vector of diameters <21-08-25> 
# FIX: This could be done by passing a Sampleable <21-08-25> 
function SyntheticData.beads(N::Int, d::Length, Δxy::PixelSize, wh::Size; 
        α::PerLength=0u"nm^-1", 
        spacing::Length = d * 2.5,  
        xdist=Uniform(0,wh[1]), 
        ydist=Uniform(0,wh[2]),
        maxiters::Int=30
    )
    positions = NTuple{2,Float64}[]
    iters = 0
    while length(positions) < N
        n = length(positions)
        positions = append!(positions, zip(rand(xdist, N-n), rand(ydist, N-n)))
        # PERF: Inefficient, because we are filtering the whole array and removing all the close points even though, one
        # of them may remain. It may be faster to iterate over the vector and remove the invalid points one by one. <21-08-25> 
        filter!(positions) do p
            !any(norm((p .- other) .* Δxy) < spacing for other in positions if other != p)
        end
        iters += 1
        if iters >= maxiters
            throw(error("Failed to generate $N bead positions in $maxiters attempts. Try to decrease `N` ($N) or the `spacing` ($spacing)"))
        end
    end
    return Beads(positions, Δxy, wh, d, α)
end

function SyntheticData.groundtruth(bs::Beads)
    buf = zeros(bs.wh)
    for p in bs.positions
		px_position = round.(Int, p)
		subpx_position = p .- px_position
		b = TF.Estimation.bead(bs.diameter,bs.Δxy;α=bs.α_evanescent, position=subpx_position)
        b_window = CartesianIndices(b) .+ CartesianIndex(px_position)
        buf_overlap = intersect(CartesianIndices(buf), CartesianIndices(b) .+ CartesianIndex(px_position))
        b_overlap = intersect(buf_overlap .- CartesianIndex(px_position), CartesianIndices(b))
        buf[buf_overlap] .+= b[b_overlap]
    end
    return SpatialArray(buf, bs.Δxy)
end

end
