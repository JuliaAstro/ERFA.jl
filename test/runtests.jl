include("../src/ERFA.jl")
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

## test from t_erfa_c.c
# eraCal2jd
dmj0, dmj = eraCal2jd(2003, 6, 1)
@test_approx_eq_eps dmj0  2400000.5 1e-9
@test_approx_eq_eps dmj  52791.0 1e-9

# eraDat
d = eraDat(2003, 6, 1, 0.0)
@test_approx_eq_eps d  32.0 1e-9
d = eraDat(2008, 1, 17, 0.0)
@test_approx_eq_eps d  33.0 1e-9

# eraD2dtf
y,m,d,H,M,S,F = eraD2dtf("UTC", 5, 2400000.5, 49533.99999)
@test (y,m,d,H,M,S,F) == (1994,6,30,23,59,60,13599)

# eraDtf2d
jd1, jd2 = eraDtf2d("UTC", 1994, 6, 30, 23, 59, 60.13599)
@test_approx_eq_eps jd1+jd2  2449534.49999 1e-6

# eraEe00a
ee = eraEe00a(2400000.5, 53736.0)
@test_approx_eq_eps ee  -0.8834192459222588227e-5 1e-18

# eraEe00b
ee = eraEe00b(2400000.5, 53736.0)
@test_approx_eq_eps ee  -0.8835700060003032831e-5 1e-18

# eraEe06a 
ee = eraEe06a(2400000.5, 53736.0)
@test_approx_eq_eps ee  -0.8834195072043790156e-5 1e-15

# eraEpb
b = eraEpb(2415019.8135, 30103.18648)
@test_approx_eq_eps b  1982.418424159278580 1e-12

# eraEpb2jd
dj0,dj1 = eraEpb2jd(1957.3)
@test_approx_eq_eps dj0  2400000.5 1e-9
@test_approx_eq_eps dj1  35948.1915101513 1e-9

# eraEpj
j = eraEpj(2451545, -7392.5)
@test_approx_eq_eps j  1979.760438056125941 1e-12

# eraEpj2jd
dj0,dj1 = eraEpj2jd(1996.8)
@test_approx_eq_eps dj0  2400000.5 1e-9
@test_approx_eq_eps dj1  50375.7 1e-9

# eraEqeq94
ee = eraEqeq94(2400000.5, 41234.0)
@test_approx_eq_eps ee  0.5357758254609256894e-4 1e-17

# eraEra00
era = eraEra00(2400000.5, 54388.0)
@test_approx_eq_eps era  0.4022837240028158102 1e-12

#eraFad03
d = eraFad03(0.80)
@test_approx_eq_eps d  1.946709205396925672 1e-12

# eraFae03
e = eraFae03(0.80)
@test_approx_eq_eps e  1.744713738913081846 1e-12

# eraFaf03
f = eraFaf03(0.80)
@test_approx_eq_eps f  0.2597711366745499518 1e-12

# eraFaju03
l = eraFaju03(0.80)
@test_approx_eq_eps l  5.275711665202481138 1e-12

# eraFal03
l = eraFal03(0.80)
@test_approx_eq_eps l  5.132369751108684150 1e-12

# eraFalp03
lp = eraFalp03(0.80)
@test_approx_eq_eps lp  6.226797973505507345 1e-12

# eraFama03
l = eraFama03(0.80)
@test_approx_eq_eps l  3.275506840277781492 1e-12

# eraFame03
l = eraFame03(0.80)
@test_approx_eq_eps l  5.417338184297289661 1e-12

# eraFane03
l = eraFane03(0.80)
@test_approx_eq_eps l  2.079343830860413523 1e-12

# eraFaom03
l = eraFaom03(0.80)
@test_approx_eq_eps l  -5.973618440951302183 1e-12

# eraFapa03
l = eraFapa03(0.80)
@test_approx_eq_eps l  0.1950884762240000000e-1 1e-12

# eraFasa03
l = eraFasa03(0.80)
@test_approx_eq_eps l  5.371574539440827046 1e-12

# eraFaur03
l = eraFaur03(0.80)
@test_approx_eq_eps l  5.180636450180413523 1e-12

# eraFave03
l = eraFave03(0.80)
@test_approx_eq_eps l  3.424900460533758000 1e-12

# eraGmst82
g = eraGmst82(2400000.5, 53736.0)
@test_approx_eq_eps g  1.754174981860675096 1e-14

# eraGst00b
g = eraGst00b(2400000.5, 53736.0)
@test_approx_eq_eps g  1.754166136510680589 1e-14

# eraGst94
g = eraGst94(2400000.5, 53736.0)
@test_approx_eq_eps g  1.754166136020645203 1e-14

# eraJd2cal
y, m, d, fd = eraJd2cal(2400000.5, 50123.9999)
@test (y, m, d) == (1996, 2, 10)
@test_approx_eq_eps fd  0.9999 1e-7

# eraJdcalf
y, m, d, fd = eraJdcalf(4, 2400000.5, 50123.9999)
@test (y, m, d, fd) == (1996, 2, 10, 9999)

# eraNum00a
rmatn = eraNum00a(2400000.5, 53736.0)
@test_approx_eq_eps rmatn[1]  0.9999999999536227949 1e-12
@test_approx_eq_eps rmatn[2]  0.8836238544090873336e-5 1e-12
@test_approx_eq_eps rmatn[3]  0.3830835237722400669e-5 1e-12
@test_approx_eq_eps rmatn[4]  -0.8836082880798569274e-5 1e-12
@test_approx_eq_eps rmatn[5]  0.9999999991354655028 1e-12
@test_approx_eq_eps rmatn[6]  -0.4063240865362499850e-4 1e-12
@test_approx_eq_eps rmatn[7]  -0.3831194272065995866e-5 1e-12
@test_approx_eq_eps rmatn[8]  0.4063237480216291775e-4 1e-12
@test_approx_eq_eps rmatn[9]  0.9999999991671660338 1e-12

# eraNum00b
rmatn = eraNum00b(2400000.5, 53736.0)
@test_approx_eq_eps rmatn[1]  0.9999999999536069682 1e-12
@test_approx_eq_eps rmatn[2]  0.8837746144871248011e-5 1e-12
@test_approx_eq_eps rmatn[3]  0.3831488838252202945e-5 1e-12
@test_approx_eq_eps rmatn[4]  -0.8837590456632304720e-5 1e-12
@test_approx_eq_eps rmatn[5]  0.9999999991354692733 1e-12
@test_approx_eq_eps rmatn[6]  -0.4063198798559591654e-4 1e-12
@test_approx_eq_eps rmatn[7]  -0.3831847930134941271e-5 1e-12
@test_approx_eq_eps rmatn[8]  0.4063195412258168380e-4 1e-12
@test_approx_eq_eps rmatn[9]  0.9999999991671806225 1e-12

# eraNum06a
rmatn = eraNum06a(2400000.5, 53736.)
@test_approx_eq_eps rmatn[1]  0.9999999999536227668 1e-12
@test_approx_eq_eps rmatn[2]  0.8836241998111535233e-5 1e-12
@test_approx_eq_eps rmatn[3]  0.3830834608415287707e-5 1e-12
@test_approx_eq_eps rmatn[4]  -0.8836086334870740138e-5 1e-12
@test_approx_eq_eps rmatn[5]  0.9999999991354657474 1e-12
@test_approx_eq_eps rmatn[6]  -0.4063240188248455065e-4 1e-12
@test_approx_eq_eps rmatn[7]  -0.3831193642839398128e-5 1e-12
@test_approx_eq_eps rmatn[8]  0.4063236803101479770e-4 1e-12
@test_approx_eq_eps rmatn[9]  0.9999999991671663114 1e-12

# eraNumat
epsa =  0.4090789763356509900
dpsi = -0.9630909107115582393e-5
deps =  0.4063239174001678826e-4
rmatn = eraNumat(epsa, dpsi, deps)
@test_approx_eq_eps rmatn[1]  0.9999999999536227949  1e-12
@test_approx_eq_eps rmatn[2]  0.8836239320236250577e-5  1e-12
@test_approx_eq_eps rmatn[3]  0.3830833447458251908e-5  1e-12
@test_approx_eq_eps rmatn[4]  -0.8836083657016688588e-5  1e-12
@test_approx_eq_eps rmatn[5]  0.9999999991354654959  1e-12
@test_approx_eq_eps rmatn[6]  -0.4063240865361857698e-4  1e-12
@test_approx_eq_eps rmatn[7]  -0.3831192481833385226e-5  1e-12
@test_approx_eq_eps rmatn[8]  0.4063237480216934159e-4  1e-12
@test_approx_eq_eps rmatn[9]  0.9999999991671660407  1e-12

# eraNut00a
dpsi, deps = eraNut00a(2400000.5, 53736.0)
@test_approx_eq_eps dpsi  -0.9630909107115518431e-5 1e-13
@test_approx_eq_eps deps  0.4063239174001678710e-4 1e-13

# eraNut00b
dpsi, deps = eraNut00b(2400000.5, 53736.0)
@test_approx_eq_eps dpsi  -0.9632552291148362783e-5 1e-13
@test_approx_eq_eps deps  0.4063197106621159367e-4 1e-13

# eraNut06a
dpsi, deps = eraNut06a(2400000.5, 53736.0)
@test_approx_eq_eps dpsi  -0.9630912025820308797e-5 1e-13
@test_approx_eq_eps deps  0.4063238496887249798e-4 1e-13

# eraNut80
dpsi, deps = eraNut80(2400000.5, 53736.0)
@test_approx_eq_eps dpsi  -0.9643658353226563966e-5 1e-13
@test_approx_eq_eps deps  0.4060051006879713322e-4 1e-13

# eraNutm80
rmatn = eraNutm80(2400000.5, 53736.)
@test_approx_eq_eps rmatn[1]  0.9999999999534999268 1e-12
@test_approx_eq_eps rmatn[2]  0.8847935789636432161e-5 1e-12
@test_approx_eq_eps rmatn[3]  0.3835906502164019142e-5 1e-12
@test_approx_eq_eps rmatn[4]  -0.8847780042583435924e-5 1e-12
@test_approx_eq_eps rmatn[5]  0.9999999991366569963 1e-12
@test_approx_eq_eps rmatn[6]  -0.4060052702727130809e-4 1e-12
@test_approx_eq_eps rmatn[7]  -0.3836265729708478796e-5 1e-12
@test_approx_eq_eps rmatn[8]  0.4060049308612638555e-4 1e-12
@test_approx_eq_eps rmatn[9]  0.9999999991684415129 1e-12

# eraObl06
obl = eraObl06(2400000.5, 54388.0)
@test_approx_eq_eps obl 0.4090749229387258204 1e-16

# eraObl80
obl = eraObl80(2400000.5, 54388.0)
@test_approx_eq_eps obl 0.409075134764381621 1e-16

# eraPmat00
rbp = eraPmat00(2400000.5, 50123.9999)
@test_approx_eq_eps rbp[1]  0.9999995505175087260 1e-12
@test_approx_eq_eps rbp[2]  0.8695405883617884705e-3 1e-14
@test_approx_eq_eps rbp[3]  0.3779734722239007105e-3 1e-14
@test_approx_eq_eps rbp[4]  -0.8695405990410863719e-3 1e-14
@test_approx_eq_eps rbp[5]  0.9999996219494925900 1e-12
@test_approx_eq_eps rbp[6]  -0.1360775820404982209e-6 1e-14
@test_approx_eq_eps rbp[7]  -0.3779734476558184991e-3 1e-14
@test_approx_eq_eps rbp[8]  -0.1925857585832024058e-6 1e-14
@test_approx_eq_eps rbp[9]  0.9999999285680153377 1e-12

# eraPmat06
rbp = eraPmat06(2400000.5, 50123.9999)
@test_approx_eq_eps rbp[1]  0.9999995505176007047 1e-12
@test_approx_eq_eps rbp[2]  0.8695404617348208406e-3 1e-14
@test_approx_eq_eps rbp[3]  0.3779735201865589104e-3 1e-14
@test_approx_eq_eps rbp[4]  -0.8695404723772031414e-3 1e-14
@test_approx_eq_eps rbp[5]  0.9999996219496027161 1e-12
@test_approx_eq_eps rbp[6]  -0.1361752497080270143e-6 1e-14
@test_approx_eq_eps rbp[7]  -0.3779734957034089490e-3 1e-14
@test_approx_eq_eps rbp[8]  -0.1924880847894457113e-6 1e-14
@test_approx_eq_eps rbp[9]  0.9999999285679971958 1e-12

# eraPmat76
rmatp = eraPmat76(2400000.5, 50123.9999)
@test_approx_eq_eps rmatp[1]  0.9999995504328350733 1e-12
@test_approx_eq_eps rmatp[2]  0.8696632209480960785e-3 1e-14
@test_approx_eq_eps rmatp[3]  0.3779153474959888345e-3 1e-14
@test_approx_eq_eps rmatp[4]  -0.8696632209485112192e-3 1e-14
@test_approx_eq_eps rmatp[5]  0.9999996218428560614 1e-12
@test_approx_eq_eps rmatp[6]  -0.1643284776111886407e-6 1e-14
@test_approx_eq_eps rmatp[7]  -0.3779153474950335077e-3 1e-14
@test_approx_eq_eps rmatp[8]  -0.1643306746147366896e-6 1e-14
@test_approx_eq_eps rmatp[9]  0.9999999285899790119 1e-12

# eraTaitt
t1, t2 = eraTaitt(2453750.5, 0.892482639)
@test_approx_eq_eps t1  2453750.5 1e-6
@test_approx_eq_eps t2  0.892855139 1e-12

# eraTaiut1
u1, u2 = eraTaiut1(2453750.5, 0.892482639, -32.6659)
@test_approx_eq_eps u1  2453750.5 1e-6
@test_approx_eq_eps u2  0.8921045614537037037 1e-12

# eraTaiutc
u1, u2 = eraTaiutc(2453750.5, 0.892482639)
@test_approx_eq_eps u1  2453750.5 1e-6
@test_approx_eq_eps u2  0.8921006945555555556 1e-12

# eraTcbtdb
b1, b2 = eraTcbtdb(2453750.5, 0.893019599)
@test_approx_eq_eps b1  2453750.5 1e-6
@test_approx_eq_eps b2  0.8928551362746343397 1e-12

# eraTcgtt
t1, t2 = eraTcgtt(2453750.5,  0.892862531)
@test_approx_eq_eps t1  2453750.5 1e-6
@test_approx_eq_eps t2  0.8928551387488816828 1e-12

# eraTdbtcb
b1, b2 = eraTdbtcb(2453750.5, 0.892855137)
@test_approx_eq_eps b1  2453750.5 1e-6
@test_approx_eq_eps b2  0.8930195997253656716 1e-12

# eraTdbtt
t1, t2 = eraTdbtt(2453750.5,  0.892855137, -0.000201)
@test_approx_eq_eps t1  2453750.5 1e-6
@test_approx_eq_eps t2  0.8928551393263888889 1e-12

# eraTf2d
d = eraTf2d('+',23,55,10.9)
@test_approx_eq_eps d  0.9966539351851851852 1e-12

# eraTttai
t1, t2 = eraTttai(2453750.5, 0.892482639)
@test_approx_eq_eps t1  2453750.5 1e-6
@test_approx_eq_eps t2  0.892110139 1e-12

# eraTttcg
t1, t2 = eraTttcg(2453750.5, 0.892482639)
@test_approx_eq_eps t1  2453750.5 1e-6
@test_approx_eq_eps t2  0.8924900312508587113 1e-12

# eraTttdb
t1, t2 = eraTttdb(2453750.5, 0.892855139, -0.000201)
@test_approx_eq_eps t1  2453750.5 1e-6
@test_approx_eq_eps t2  0.8928551366736111111 1e-12

# eraTtut1
t1, t2 = eraTtut1(2453750.5, 0.892855139, 64.8499)
@test_approx_eq_eps t1  2453750.5 1e-6
@test_approx_eq_eps t2  0.8921045614537037037 1e-12

# eraUt1tai
a1, a2 = eraUt1tai(2453750.5, 0.892104561, -32.6659)
@test_approx_eq_eps a1  2453750.5 1e-6
@test_approx_eq_eps a2  0.8924826385462962963 1e-12

# eraUt1tt
a1, a2 = eraUt1tt(2453750.5, 0.892104561, 64.8499)
@test_approx_eq_eps a1  2453750.5 1e-6
@test_approx_eq_eps a2  0.8928551385462962963 1e-15

# eraUt1utc
a1, a2 = eraUt1utc(2453750.5, 0.892104561, 0.3341)
@test_approx_eq_eps a1  2453750.5 1e-6
@test_approx_eq_eps a2  0.8921006941018518519 1e-13

# eraUtctai
u1, u2 = eraUtctai(2453750.5, 0.892100694)
@test_approx_eq_eps u1  2453750.5 1e-6
@test_approx_eq_eps u2  0.8924826384444444444 1e-13

# eraUtcut1
u1, u2 = eraUtcut1(2453750.5, 0.892100694, 0.3341)
@test_approx_eq_eps u1  2453750.5 1e-6
@test_approx_eq_eps u2  0.8921045608981481481 1e-13
