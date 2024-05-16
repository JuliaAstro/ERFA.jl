using Documenter, ERFA

include("pages.jl")

makedocs(
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    sitename = "ERFA.jl",
    authors = "The JuliaAstro Contributors",
    pages = pages,
    doctest = true,
)

deploydocs(
    repo = "github.com/JuliaAstro/ERFA.jl.git",
    target = "build",
)