using TransferFunctions
using Documenter, DocumenterCitations, DocumenterInterLinks

links = InterLinks(
    "Julia" => "https://docs.julialang.org/en/v1/",
    "Unitful" => "https://painterqubits.github.io/Unitful.jl/stable/",
    "ImageFiltering" => "https://juliaimages.org/ImageFiltering.jl/stable/"
)

DocMeta.setdocmeta!(TransferFunctions, :DocTestSetup, :(
        using TransferFunctions;
        using TestImages
    ); recursive=true)

bib = CitationBibliography(
    joinpath(@__DIR__, "src", "refs.bib");
    # style=:authoryear
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
            "pages/01-interface.md",
            "pages/02-transfer-functions.md",
            "pages/03-estimation.md",
            "pages/04-apodization.md"
        ],
        "Reference" => [
            "Internals" => [
                "Array types" => "pages/internals/01-arrays.md"
            ],
            "API Index" => "pages/05-apireference.md",
            "Bibliography" => "pages/06-bibliography.md"
        ]],
    plugins=[bib, links],
    # NOTE: doctesting is done in the `runtests.jl` so it is not necessary to do here
    doctest=false
)

deploydocs(;
    repo="github.com/kunzaatko/TransferFunctions.jl",
    devbranch="trunk"
)
