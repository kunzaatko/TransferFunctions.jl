using TransferFunctions
using Documenter, DocumenterCitations, DocumenterInterLinks, MakieMaestro

MakieMaestro.Themes.width!(25u"cm")

makie_doc_blocks = MakieMaestro.MakieDocBlocks(;
    formats=[:png]
)

links = InterLinks(
    "Julia" => "https://docs.julialang.org/en/v1/",
    "Unitful" => "https://juliaphysics.github.io/Unitful.jl/stable/",
    "ImageFiltering" => "https://juliaimages.org/ImageFiltering.jl/stable/",
    # TODO: When `objects.inv` are generated for the stable branch, change to it <19-08-25> 
    "Rotations" => "https://juliageometry.github.io/Rotations.jl/dev/",
    "ComponentArrays" => "https://docs.sciml.ai/ComponentArrays/stable/"
)

DocMeta.setdocmeta!(TransferFunctions, :DocTestSetup, :(
        include(joinpath(@__DIR__, "../test", "doctestsetup.jl"))
    ); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
)

makedocs(;
    modules=[TransferFunctions],
    authors="Martin Kunz <martinkunz@email.cz> and contributors",
    repo="https://github.com/kunzaatko/TransferFunctions.jl/blob/{commit}{path}#{line}",
    sitename="TransferFunctions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kunzaatko.github.io/TransferFunctions.jl",
        edit_link="trunk",
        assets=String[]
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "pages/manual/01-interface.md",
            "pages/manual/02-linear-transfer-functions.md",
            "pages/manual/03-nonlinear-transfer-functions.md",
            "pages/manual/04-restoration.md",
        ],
        "Transfer Functions" => [
            "Models" => [
                "Airy Disc" => "pages/transfer-function-types/airy-disc.md",
                "Gaussian" => "pages/transfer-function-types/gaussian.md",
                "Born & Wolf" => "pages/transfer-function-types/born-wolf.md",
                "Gibson & Lanni" => "pages/transfer-function-types/gibson-lanni.md",
                "Circular Pupil" => "pages/transfer-function-types/circular-pupil.md",
            ],
        ],
        "Array Types" => [
            "pages/arrays/sampled-arrays.md"
            "pages/arrays/border-arrays.md"
            "pages/arrays/tapered-arrays.md"
            "pages/arrays/circulant-tensors.md"
            "pages/arrays/filtering-matrices.md"
            "pages/arrays/reflected-arrays.md"
        ],
        "Apodization" => "pages/apodization.md",
        "Filtering" => "pages/filter.md",
        "Reference" => [
            "API Index" => "pages/apireference.md",
            "Bibliography" => "pages/bibliography.md"
        ]],
    plugins=[bib, links],
    warnonly=[:missing_docs],
    doctest=false # tests run in `test/runtests.jl`
)

deploydocs(;
    repo="github.com/kunzaatko/TransferFunctions.jl",
    devbranch="trunk"
)

if !haskey(ENV, "GITHUB_ACTIONS")
    build_path = joinpath(@__DIR__, "build")
    cached_path = joinpath(@__DIR__, "cached")
    @info "Making cached docs at $cached_path"
    cp(build_path, cached_path; force=true)
end
