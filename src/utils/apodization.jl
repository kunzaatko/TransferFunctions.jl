# TODO: Use reinterpret instead of `T.` for the changes of type in the places that it is used. <09-09-24> 
using ImageFiltering: AbstractBorder, borderinstance, BorderSpecAny


# TODO: This should really be a separate package <26-08-24> 
# FIX: Type stability <13-12-23> 
# FIX: Domain error, when ∉ (-1,1) <13-12-23> 
# https://mathworld.wolfram.com/ApodizationFunction.html
# TODO: Add documentation about what is an apodization function and how it is used. "An apodization function is ... zero-phase function ... Instrument function ... Links" <18-11-24> 
@doc raw"""
    Apodization

Abstract type for apodization functions.
"""
abstract type Apodization end
Broadcast.broadcastable(a::Apodization) = Ref(a)
apodization(apo::Apodization, x::Real, halfwidth::Int) = apodization(apo, x / halfwidth)

include("apodization-types.jl")

# FIX: The call stack must be rewritten to construct the arguments of the function from the beginning <02-09-24> 

# TODO: There should be a mutating method and a non-mutating method that allocates the output array <26-08-24> 
# TODO: Should allow padding with some scheme from `ImageFiltering.jl` before the tapering for an unobscured data
# behaviour  <26-08-24> 
# FIX: This API should be slightly different. Right now it requires to specify the apodization. Usually the user can be
# satisfied with the default `TransferFunctions.Cosine`. The `width` parameter should be part of the `apo` instance.
# When only the width is set, the default apodization should be instantiated. The tendency should be for lower
# parameter methods be easier and just call the general methods, which are the more flexible ones. <17-07-24> 
# TODO: Document: how the widths arguments work <02-09-24> 
@doc raw"""
    taperedges([apo], A, width::Int, [border="replicate"]; dims=1:N)
    taperedges([apo], A, (w1,...,wM), [border="replicate"]; dims =1:M)
    taperedges([apo], A, ((w1_start,...,wM_start), (w1_end,...,wM_end)), [border="replicate"]; dims=1:M)

Taper a given width of the array `A`'s edges using [`apo::Apodization`](@ref Apodization)

- `apo`: [Apodization](@ref) function to use for the tapering. Default is `Cosine`.
- `border="replicate"` : Similarly as with `imfilter`, the array `A` may be extended to avoid loss of the data at the edges
- `width::NTuple{M,Int}` or `Int` or `Tuple{NTuple{M,Int},NTuple{M,Int}}`: The width of tapering at the edges.
- `dims::NTuple{M,Int}`: The dimensions along which to taper the edges.
"""
taperedges

function taperedges( # STEP 1A: Fill the apodization type
    A::AbstractArray{<:Number,N},
    args...;
    kwargs...
) where {N}
    return taperedges(Cosine(), A, args...; kwargs...)
end
function taperedges( # STEP 2A: Fill from `width` single width
    apo::Apodization,
    A::AbstractArray{<:Number,N},
    width::Int,
    args...;
    # FIX: Instead of this, should be default `Colon()` <02-09-24> 
    dims::NTuple{M,Int}=Tuple(1:N), # TODO: Document by default taper along all the dimensions <03-09-24> 
    kwargs...
) where {N,M}
    # TODO: Test whether this works for permuted order of kwargs... I.e. whether dims must be supplied as the first
    # argument or not. Otherwise, it must be done by testing if kwargs has dims in it... <02-09-24> 
    return taperedges(apo, A, (Tuple(fill(width, M)), Tuple(fill(width, M))), args...; dims, kwargs...)
end
function taperedges( # STEP 2B: Fill from `width`s for each dimension
    apo::Apodization,
    A::AbstractArray{<:Number,N},
    edge_widths::NTuple{M,Int},
    args...;
    kwargs...
) where {N,M}
    return taperedges(apo, A, (edge_widths, edge_widths), args...; kwargs...)
end
function taperedges( # STEP 3: Fill in the default `border`
    apo::Apodization,
    A::AbstractArray{<:Number,N},
    widths::Tuple{NTuple{M,Int},NTuple{M,Int}},
    args...;
    kwargs...
) where {N,M}
    return taperedges(apo, A, widths, "replicate", args...; kwargs...)
end
function taperedges( # STEP 4: Create a border instance
    apo::Apodization,
    A::AbstractArray{<:Number,N},
    widths::Tuple{NTuple{M,Int},NTuple{M,Int}},
    border::AbstractString,
    args...;
    kwargs...
) where {N,M}
    return taperedges(apo, A, widths, borderinstance(border), args...; kwargs...)
end
function taperedges( # STEP 5: Set the border sizes
    apo::Apodization,
    A::AbstractArray{<:Number,N},
    widths::Tuple{NTuple{M,Int},NTuple{M,Int}},
    border::BorderSpecAny,
    args...;
    dims=Tuple(1:M),
    kwargs...
) where {N,M}
    # FIX: Must be tested here and in the final because otherwise we would index out of bounds <09-09-24> 
    (M > N || maximum(dims) > N) && throw(ArgumentError("The number of dimensions must be less than or equal to the number of axes."))

    concrete_border = full_padding_border(border, widths, dims, N)
    return taperedges(apo, A, widths, concrete_border, args...; kwargs...)
end

# TODO: Test this <09-09-24> 
function full_padding_border(border::BorderSpecAny, widths::Tuple{NTuple{M,Int},NTuple{M,Int}}, dims::NTuple{M,Int}, ndims::Int) where {M}
    all_left_widths, all_right_widths = zeros(Int, ndims), zeros(Int, ndims)
    foreach(dims, widths[1], widths[2]) do dim, w_left, w_right
        all_left_widths[dim] = w_left
        all_right_widths[dim] = w_right
    end
    all_left_widths, all_right_widths = Tuple(all_left_widths), Tuple(all_right_widths)
    if border isa Pad
        return Pad(border.style, all_left_widths, all_right_widths)
    elseif border isa Fill
        return Fill(border.value, all_left_widths, all_right_widths)
    elseif border isa Inner
        return Inner(all_left_widths, all_right_widths)
    else
        throw(ErrorException("`NA` and `NoPad` borders should not occur here. Type is $(typeof(border))."))
    end
end

# TODO: Perhaps there could be an argument to make the array odd sized for the Fourier transform. Since we do not have
# to have 0 at both edges. For the signal to be periodic, only one edge to be 0 is sufficient. An odd size is beneficial
# for a Fourier transform. <10-09-24> 
function taperedges( # FINAL # TODO: Instead of this should be something like `_taperedges` function <02-09-24> 
    apo::Apodization,
    A::AbstractArray{<:Number,N},
    width::Tuple{NTuple{M,Int},NTuple{M,Int}},
    border::AbstractBorder;
    dims::NTuple{M,Int}=Tuple(1:M) # TODO: Document that the number of dimensions is by default taken as the first M 
    # dimensions where M is the number of widths supplied <kunzaatko> 
) where {N,M}
    # TODO: How can one handle the dimensions and the various types that can define them (such as Colon())... Ask on
    # discourse and implement. <02-09-24> 
    (M > N || maximum(dims) > N) && throw(ArgumentError("The number of dimensions must be less than or equal to the number of axes."))

    A = padarray(A, border)

    for (dim, low, high) in zip(dims, width[1], width[2])
        lowedge = range(-1, 0; length=low + 1)[begin:(end-1)]
        for (e, i) in enumerate(firstindex(axes(A, dim)):(firstindex(axes(A, dim))+(low-1)))
            selectdim(A, dim, i) .*= apodization(apo, lowedge[e])
        end
        highedge = range(0, 1; length=high + 1)[(begin+1):end]
        @assert length(highedge) == length((lastindex(axes(A, dim))-(high-1)):lastindex(axes(A, dim)))
        for (e, i) in enumerate((lastindex(axes(A, dim))-(high-1)):lastindex(axes(A, dim)))
            selectdim(A, dim, i) .*= apodization(apo, highedge[e])
        end
    end
    return A
end

# TODO: There should be a mutating method and a non-mutating method that allocates the output array <26-08-24> 
# TODO: Add documentation <21-12-23> 
@doc raw"""
    apodize(apo::Apodization, A::AbstractArray{<:Number,N}, cutoff::Real, width::Real) where {N}

Apply apodization to the input array `A`.

# Arguments
- `apo::Apodization`: The type of apodization to apply.
- `A::AbstractArray{<:Number,N}`: The input array to be apodized.
- `cutoff::Real`: The cutoff frequency for the apodization.
- `width::Real`: The width of the transition region for the apodization.
"""
function apodize(
    apo::Apodization,
    A::AbstractArray{<:Number,N},
    cutoff::Real, # Tuple{NTuple{M,Int},NTuple{M,Int}}, # FIX: Generalize to non-symmetric cut-offs <21-12-23> 
    width::Real, # FIX: Generalize to non-symmetric widths <21-12-23> 
    dims::NTuple{M,Int}=Tuple(1:N) # FIX: Abstract over dimensions <22-12-23> 
) where {N,M}
    # PERF: Should be done with no allocation... This is the KISS solution  
    # FIX: Work for all dimensions  
    rs = [hypot(abs(x), abs(y)) for x in fftfreq(size(A, 1), size(A, 1)), y in fftfreq(size(A, 2), size(A, 2))]
    coefs = ones(eltype(A), size(rs)...)
    coefs[rs.>=cutoff] .= 0
    coefs[rs.<=(cutoff-width)] .= 1
    coefs[cutoff.>rs.>(cutoff-width)] .= map(r -> apodization(apo, r - (cutoff - width), width), rs[cutoff.>rs.>(cutoff-width)])
    A .* coefs
end
