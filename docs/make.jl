using GinzburgLandau
using Documenter

DocMeta.setdocmeta!(GinzburgLandau, :DocTestSetup, :(using GinzburgLandau); recursive=true)

makedocs(;
    modules=[GinzburgLandau],
    authors="Maximilian Gelbrecht <maximilian.gelbrecht@posteo.de> and contributors",
    repo="https://github.com/maximilian-gelbrecht/GinzburgLandau.jl/blob/{commit}{path}#{line}",
    sitename="GinzburgLandau.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://maximilian-gelbrecht.github.io/GinzburgLandau.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Reference" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/maximilian-gelbrecht/GinzburgLandau.jl",
    devbranch="main",
)
