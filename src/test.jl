include("ERFA.jl")
using ERFA
using .Test

u1,u2 = eraDtf2d("UTC", 2010, 7, 24, 11, 18, 7.318)
a1,a2 = eraUtctai(u1, u2)
t1,t2 = eraTaitt(a1, a2)
@test eraD2dtf("tt", 3, t1, t2) == (2010,7,24,11,19,13,502)
