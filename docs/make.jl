using Documenter, ERFA

makedocs(
    format = :html,
    sitename = "ERFA.jl",
    authors = "The ERFA.jl Maintainers",
    pages = [
        "Home" => "index.md",
        "API" => "api.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaAstro/ERFA.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    julia = "0.6",
)
