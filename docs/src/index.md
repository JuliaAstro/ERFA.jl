ERFA.jl
=======

Julia wrapper for [liberfa](https://github.com/liberfa/erfa).

Installation
------------

```julia-repl
julia> import Pkg; Pkg.add("ERFA")
```

Example
-------

```jldoctest
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
