using TransferFunctions
using Documenter

DocMeta.setdocmeta!(TransferFunctions, :DocTestSetup, :(using TransferFunctions); recursive=true)

makedocs(;
    modules=[TransferFunctions],
    authors="Martin Kunz <martinkunz@email.cz> and contributors",
    repo="https://github.com/kunzaatko/TransferFunctions.jl/blob/{commit}{path}#{line}",
    sitename="TransferFunctions.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://kunzaatko.github.io/TransferFunctions.jl",
        edit_link="trunk",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/kunzaatko/TransferFunctions.jl",
    devbranch="trunk",
)
