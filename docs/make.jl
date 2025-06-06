using Documenter, ERFA

include("pages.jl")

makedocs(;
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold_warn = 400 * 1024, # 400 KiB
        size_threshold = 800 * 1024, # 800 KiB
        canonical = "https://juliaastro.org/ERFA/stable/",
    ),
    sitename = "ERFA.jl",
    authors = "The JuliaAstro Contributors",
    pages = pages,
    doctest = true,
)

deploydocs(;
    repo = "github.com/JuliaAstro/ERFA.jl.git",
    target = "build",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)
