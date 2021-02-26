# ERFA.jl

*Julia wrapper for [liberfa](https://github.com/liberfa/erfa)*

[![Build Status](https://github.com/juliaastro/ERFA.jl/workflows/CI/badge.svg)](https://github.com/juliaastro/ERFA.jl/actions)
[![Coverage](https://codecov.io/gh/juliaastro/ERFA.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/juliaastro/ERFA.jl)

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaastro.github.io/ERFA.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaastro.github.io/ERFA.jl/dev)


## Installation

```julia
julia> Pkg.add("ERFA")
```

## Example

```julia
julia> using ERFA

julia> u1,u2 = ERFA.dtf2d("UTC", 2010, 7, 24, 11, 18, 7.318)
(2.4554015e6,0.47091803240740737)

julia> a1,a2 = ERFA.utctai(u1, u2)
(2.4554015e6,0.4713115509259259)

julia> t1,t2 = ERFA.taitt(a1, a2)
(2.4554015e6,0.4716840509259259)

julia> ERFA.d2dtf("tt", 3, t1, t2)
(2010,7,24,11,19,13,502)
```

## Documentation

Please refer to the [documentation](https://juliaastro.github.io/ERFA.jl/stable) for additional
information.

