# ERFA.jl

*Julia wrapper for [liberfa](https://github.com/liberfa/erfa)*

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.github.io/ERFA/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.github.io/ERFA.jl/dev)

[![CI](https://github.com/JuliaAstro/ERFA.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaAstro/ERFA.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/juliaastro/ERFA.jl/graph/badge.svg?token=q5WlA4TVIZ)](https://codecov.io/gh/juliaastro/ERFA.jl)

## Installation

```julia
julia> import Pkg; Pkg.add("ERFA")
```

## Example

```julia
julia> using ERFA

julia> u1, u2 = dtf2d("UTC", 2010, 7, 24, 11, 18, 7.318)
(2.4554015e6, 0.47091803240740737)

julia> a1, a2 = utctai(u1, u2)
(2.4554015e6, 0.4713115509259259)

julia> t1, t2 = taitt(a1, a2)
(2.4554015e6, 0.4716840509259259)

julia> d2dtf("TT", 3, t1, t2)
(2010, 7, 24, 11, 19, 13, 502)
```

## Documentation

Please refer to the [documentation](https://juliaastro.github.io/ERFA.jl/stable) for additional
information.

