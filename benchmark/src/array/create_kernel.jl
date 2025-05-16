function centered_monotone_kernel(s...)
    @assert all(isodd, s)
    K = OAs.centered(Array{Float64}(undef, s))
    max_hypot = hypot(maximum.(map(x -> abs.(x), extrema.(axes(K))))...)
    K .= [cos(hypot(Tuple(i)...) ./ max_hypot) for i in CartesianIndices(K)]
    K .+= samerand(s)
    K ./= sum(K)
    return K
end
