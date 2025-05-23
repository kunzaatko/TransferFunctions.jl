function Base.showarg(io::IO, A::Flattened, toplevel)
    print(io, "flatten(")
    Base.showarg(io, parent(A), false)
    if A.outermap != _default_outer(parent(A))
        print(io, "; outer=")
        print(io, A.outermap)
        print(io, ", inner=")
        print(io, A.innermap)
    end
    print(io, ")")
    toplevel && print(io, " with eltype ", eltype(A))
end

function Base.showarg(io::IO, A::CirculantTensor, toplevel)
    print(io, "circulant(")
    Base.showarg(io, basearray(A), false)
    print(io, ", ")
    print(io, A.kern)
    print(io, ")")
    toplevel && print(io, " with eltype ", eltype(A))
end
