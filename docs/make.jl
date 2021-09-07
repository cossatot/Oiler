using Documenter, Oiler

makedocs(;
    modules=[Oiler],
    format=Documenter.HTML(assets=String[]),
    pages=[
        "Home" => "index.md",
        "Examples" => "examples.md",
        "Module" => map(
            s -> "docstrings/$(s)",
            sort(readdir(joinpath(@__DIR__, "src/docstrings")))
        ),
    ],
    # repo="https://gitlab.com/cossatot/Oiler.jl/blob/{commit}{path}#L{line}",
    sitename="Oiler.jl",
    authors="Richard Styron",
)
