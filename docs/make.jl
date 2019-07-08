using Documenter, Oiler

makedocs(;
    modules=[Oiler],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://gitlab.com/cossatot/Oiler.jl/blob/{commit}{path}#L{line}",
    sitename="Oiler.jl",
    authors="Richard Styron",
    assets=String[],
)
