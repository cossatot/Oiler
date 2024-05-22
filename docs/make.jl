using Documenter, Oiler

include("makedocs.jl")

deploydocs(
    repo="github.com/cossatot/Oiler.git",
    branch="gh-pages",
    target="build",
)