include("ERFA.jl")
using ERFA
using .Test

u1,u2 = eraDtf2d("UTC", 2010, 7, 24, 11, 18, 7.318)
a1,a2 = eraUtctai(u1, u2)
t1,t2 = eraTaitt(a1, a2)
@test eraD2dtf("tt", 3, t1, t2) == (2010,7,24,11,19,13,502)

iy = 2008; imo = 2; id = 29
ihour = 23; imin = 59; sec = 59.9
d1,d2 = eraCal2jd(iy, imo, id)
d = eraTf2d('+', ihour, imin, sec)
d2 += d
@test d1 == 2400000.5
@test_approx_eq_eps d2 54525.999999 5e-7
iy,imo,id,fd = eraJd2cal(d1, d2)
@test (iy,imo,id) == (2008,2,29)
@test_approx_eq_eps fd 0.999999 5e-7
@test eraJdcalf(3, d1, d2) == (2008,3,1,0)

d = 2457073.05631
e = eraEpb(0., d)
@test_approx_eq_eps e 2015.1365941021 5e-11
d1,d2 = eraEpb2jd(e)
d = d1 + d2
@test_approx_eq_eps d 2457073.056310000 5e-10
e = eraEpj(0., d)
@test_approx_eq_eps e 2015.1349933196 5e-11
d1,d2 = eraEpj2jd(e)
d = d1 + d2
@test_approx_eq_eps d 2457073.056310000 5e-10
