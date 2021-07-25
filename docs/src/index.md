ERFA.jl
=======

Julia wrapper for [liberfa](https://github.com/liberfa/erfa).

Installation
------------

```julia
julia> Pkg.add("ERFA")
```

Example
-------

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
