using Downloads
using JuliaFormatter

const VERSION = v"2.0.0"
const URL = "https://github.com/liberfa/erfa/releases/download/v$VERSION/erfa-$VERSION.tar.gz"
const ARCHIVE = "erfa.tar.gz"
const SRC = joinpath("erfa-$VERSION", "src")
const OUTPUT = joinpath(@__DIR__(), "..", "src", "exports.jl")

workdir = tempname()
mkdir(workdir)
cd(workdir)

Downloads.download(URL, ARCHIVE)
run(`tar zxvf $ARCHIVE`)
cd(SRC)
sources = first.(splitext.(filter(x->!startswith(x, "t_") && !startswith(x, "erfa") && endswith(x, ".c"), readdir("."))))
exports = "export " * join(sources, ", ")

open(OUTPUT, "w") do f
	println(f, "# Auto-generated from $(basename(@__FILE__)). Do not edit!\n")
	println(f, exports)
end

format(OUTPUT, BlueStyle())
