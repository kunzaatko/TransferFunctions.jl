using TransferFunctions
using Documenter, DocumenterCitations, DocumenterInterLinks

links = InterLinks(
    "Unitful" => "https://painterqubits.github.io/Unitful.jl/stable/",
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
        "Theory" => "pages/01_theory.md",
        "General Interface" => "pages/02_interface.md",
        "Transfer Functions" => [
            "Transfer Function Models" => "pages/03_model_tfs.md",
            "Transfer Functions Measurements" => "pages/04_measured_tfs.md",
            "Sampled Transfer Functions" => "pages/05_sampled_tfs.md",],
        "References" => [
            "API" => "pages/06_apireference.md",
            "Bibliography" => "pages/07_bibliography.md"
        ]],
    plugins=[bib, links],
    # NOTE: doctesting is done in the `runtests.jl` so it is not necessary to do here
    doctest=false
)

deploydocs(;
    repo="github.com/kunzaatko/TransferFunctions.jl",
    devbranch="trunk"
)
