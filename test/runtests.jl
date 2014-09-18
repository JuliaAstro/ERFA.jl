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
# eraA2af
@test eraA2af(4, 2.345) == ('+', 134, 21, 30, 9706)

# eraA2tf
@test eraA2tf(4, -3.01234) == ('-', 11, 30, 22, 6484)

# eraAf2a
r = eraAf2a('-', 45, 13, 27.2)
@test_approx_eq_eps r  -0.7893115794313644842 1e-15
r = eraAf2a('+', 45, 13, 27.2)
@test_approx_eq_eps r  0.7893115794313644842 1e-15

# eraAnp
r = eraAnp(-0.1)
@test_approx_eq_eps r  6.183185307179586477 1e-15

# eraAnpm
r = eraAnpm(-4.0)
@test_approx_eq_eps r  2.283185307179586477 1e-15

# eraBp00
rb, rp,rbp = eraBp00(2400000.5, 50123.9999)
@test_approx_eq_eps rb[1]  0.9999999999999942498 1e-12
@test_approx_eq_eps rb[2]  -0.7078279744199196626e-7 1e-16
@test_approx_eq_eps rb[3]  0.8056217146976134152e-7 1e-16
@test_approx_eq_eps rb[4]  0.7078279477857337206e-7 1e-16
@test_approx_eq_eps rb[5]  0.9999999999999969484 1e-12
@test_approx_eq_eps rb[6]  0.3306041454222136517e-7 1e-16
@test_approx_eq_eps rb[7]  -0.8056217380986972157e-7 1e-16
@test_approx_eq_eps rb[8]  -0.3306040883980552500e-7 1e-16
@test_approx_eq_eps rb[9]  0.9999999999999962084 1e-12
@test_approx_eq_eps rp[1]  0.9999995504864048241 1e-12
@test_approx_eq_eps rp[2]  0.8696113836207084411e-3 1e-14
@test_approx_eq_eps rp[3]  0.3778928813389333402e-3 1e-14
@test_approx_eq_eps rp[4]  -0.8696113818227265968e-3 1e-14
@test_approx_eq_eps rp[5]  0.9999996218879365258 1e-12
@test_approx_eq_eps rp[6]  -0.1690679263009242066e-6 1e-14
@test_approx_eq_eps rp[7]  -0.3778928854764695214e-3 1e-14
@test_approx_eq_eps rp[8]  -0.1595521004195286491e-6 1e-14
@test_approx_eq_eps rp[9]  0.9999999285984682756 1e-12
@test_approx_eq_eps rbp[1]  0.9999995505175087260 1e-12
@test_approx_eq_eps rbp[2]  0.8695405883617884705e-3 1e-14
@test_approx_eq_eps rbp[3]  0.3779734722239007105e-3 1e-14
@test_approx_eq_eps rbp[4]  -0.8695405990410863719e-3 1e-14
@test_approx_eq_eps rbp[5]  0.9999996219494925900 1e-12
@test_approx_eq_eps rbp[6]  -0.1360775820404982209e-6 1e-14
@test_approx_eq_eps rbp[7]  -0.3779734476558184991e-3 1e-14
@test_approx_eq_eps rbp[8]  -0.1925857585832024058e-6 1e-14
@test_approx_eq_eps rbp[9]  0.9999999285680153377 1e-12

# eraBp06
rb, rp,rbp = eraBp06(2400000.5, 50123.9999)
@test_approx_eq_eps rb[1]  0.9999999999999942497 1e-12
@test_approx_eq_eps rb[2]  -0.7078368960971557145e-7 1e-14
@test_approx_eq_eps rb[3]  0.8056213977613185606e-7 1e-14
@test_approx_eq_eps rb[4]  0.7078368694637674333e-7 1e-14
@test_approx_eq_eps rb[5]  0.9999999999999969484 1e-12
@test_approx_eq_eps rb[6]  0.3305943742989134124e-7 1e-14
@test_approx_eq_eps rb[7]  -0.8056214211620056792e-7 1e-14
@test_approx_eq_eps rb[8]  -0.3305943172740586950e-7 1e-14
@test_approx_eq_eps rb[9]  0.9999999999999962084 1e-12
@test_approx_eq_eps rp[1]  0.9999995504864960278 1e-12
@test_approx_eq_eps rp[2]  0.8696112578855404832e-3 1e-14
@test_approx_eq_eps rp[3]  0.3778929293341390127e-3 1e-14
@test_approx_eq_eps rp[4]  -0.8696112560510186244e-3 1e-14
@test_approx_eq_eps rp[5]  0.9999996218880458820 1e-12
@test_approx_eq_eps rp[6]  -0.1691646168941896285e-6 1e-14
@test_approx_eq_eps rp[7]  -0.3778929335557603418e-3 1e-14
@test_approx_eq_eps rp[8]  -0.1594554040786495076e-6 1e-14
@test_approx_eq_eps rp[9]  0.9999999285984501222 1e-12
@test_approx_eq_eps rbp[1]  0.9999995505176007047 1e-12
@test_approx_eq_eps rbp[2]  0.8695404617348208406e-3 1e-14
@test_approx_eq_eps rbp[3]  0.3779735201865589104e-3 1e-14
@test_approx_eq_eps rbp[4]  -0.8695404723772031414e-3 1e-14
@test_approx_eq_eps rbp[5]  0.9999996219496027161 1e-12
@test_approx_eq_eps rbp[6]  -0.1361752497080270143e-6 1e-14
@test_approx_eq_eps rbp[7]  -0.3779734957034089490e-3 1e-14
@test_approx_eq_eps rbp[8]  -0.1924880847894457113e-6 1e-14
@test_approx_eq_eps rbp[9]  0.9999999285679971958 1e-12

# eraC2i00a
rc2i = eraC2i00a(2400000.5, 53736.0)
@test_approx_eq_eps rc2i[1]  0.9999998323037165557 1e-12
@test_approx_eq_eps rc2i[2]  0.5581526348992140183e-9 1e-12
@test_approx_eq_eps rc2i[3]  -0.5791308477073443415e-3 1e-12
@test_approx_eq_eps rc2i[4]  -0.2384266227870752452e-7 1e-12
@test_approx_eq_eps rc2i[5]  0.9999999991917405258 1e-12
@test_approx_eq_eps rc2i[6]  -0.4020594955028209745e-4 1e-12
@test_approx_eq_eps rc2i[7]  0.5791308472168152904e-3 1e-12
@test_approx_eq_eps rc2i[8]  0.4020595661591500259e-4 1e-12
@test_approx_eq_eps rc2i[9]  0.9999998314954572304 1e-12

# eraC2i00b
rc2i = eraC2i00b(2400000.5, 53736.0)
@test_approx_eq_eps rc2i[1]  0.9999998323040954356 1e-12
@test_approx_eq_eps rc2i[2]  0.5581526349131823372e-9 1e-12
@test_approx_eq_eps rc2i[3]  -0.5791301934855394005e-3 1e-12
@test_approx_eq_eps rc2i[4]  -0.2384239285499175543e-7 1e-12
@test_approx_eq_eps rc2i[5]  0.9999999991917574043 1e-12
@test_approx_eq_eps rc2i[6]  -0.4020552974819030066e-4 1e-12
@test_approx_eq_eps rc2i[7]  0.5791301929950208873e-3 1e-12
@test_approx_eq_eps rc2i[8]  0.4020553681373720832e-4 1e-12
@test_approx_eq_eps rc2i[9]  0.9999998314958529887 1e-12

# eraC2i06a
rc2i = eraC2i06a(2400000.5, 53736.0)
@test_approx_eq_eps rc2i[1]  0.9999998323037159379 1e-12
@test_approx_eq_eps rc2i[2]  0.5581121329587613787e-9 1e-12
@test_approx_eq_eps rc2i[3]  -0.5791308487740529749e-3 1e-12
@test_approx_eq_eps rc2i[4]  -0.2384253169452306581e-7 1e-12
@test_approx_eq_eps rc2i[5]  0.9999999991917467827 1e-12
@test_approx_eq_eps rc2i[6]  -0.4020579392895682558e-4 1e-12
@test_approx_eq_eps rc2i[7]  0.5791308482835292617e-3 1e-12
@test_approx_eq_eps rc2i[8]  0.4020580099454020310e-4 1e-12
@test_approx_eq_eps rc2i[9]  0.9999998314954628695 1e-12

# eraCal2jd
dmj0, dmj = eraCal2jd(2003, 6, 1)
@test_approx_eq_eps dmj0  2400000.5 1e-9
@test_approx_eq_eps dmj  52791.0 1e-9

# eraC2ixy
x = 0.5791308486706011000e-3
y = 0.4020579816732961219e-4
rc2i = eraC2ixy(2400000.5, 53736., x, y)
@test_approx_eq_eps rc2i[1]  0.9999998323037157138 1e-12
@test_approx_eq_eps rc2i[2]  0.5581526349032241205e-9 1e-12
@test_approx_eq_eps rc2i[3]  -0.5791308491611263745e-3 1e-12
@test_approx_eq_eps rc2i[4]  -0.2384257057469842953e-7 1e-12
@test_approx_eq_eps rc2i[5]  0.9999999991917468964 1e-12
@test_approx_eq_eps rc2i[6]  -0.4020579110172324363e-4 1e-12
@test_approx_eq_eps rc2i[7]  0.5791308486706011000e-3 1e-12
@test_approx_eq_eps rc2i[8]  0.4020579816732961219e-4 1e-12
@test_approx_eq_eps rc2i[9]  0.9999998314954627590 1e-12

# eraC2ixys
x =  0.5791308486706011000e-3
y =  0.4020579816732961219e-4
s = -0.1220040848472271978e-7
rc2i = eraC2ixys(x, y, s)
@test_approx_eq_eps rc2i[1]  0.9999998323037157138 1e-12
@test_approx_eq_eps rc2i[2]  0.5581984869168499149e-9 1e-12
@test_approx_eq_eps rc2i[3]  -0.5791308491611282180e-3 1e-12
@test_approx_eq_eps rc2i[4]  -0.2384261642670440317e-7 1e-12
@test_approx_eq_eps rc2i[5]  0.9999999991917468964 1e-12
@test_approx_eq_eps rc2i[6]  -0.4020579110169668931e-4 1e-12
@test_approx_eq_eps rc2i[7]  0.5791308486706011000e-3 1e-12
@test_approx_eq_eps rc2i[8]  0.4020579816732961219e-4 1e-12
@test_approx_eq_eps rc2i[9]  0.9999998314954627590 1e-12

# eraDat
d = eraDat(2003, 6, 1, 0.0)
@test_approx_eq_eps d  32.0 1e-9
d = eraDat(2008, 1, 17, 0.0)
@test_approx_eq_eps d  33.0 1e-9

# eraD2dtf
y,m,d,H,M,S,F = eraD2dtf("UTC", 5, 2400000.5, 49533.99999)
@test (y,m,d,H,M,S,F) == (1994,6,30,23,59,60,13599)

# eraD2tf
@test eraD2tf(4, -0.987654321) == ('-', 23, 42, 13, 3333)

# eraDtf2d
jd1, jd2 = eraDtf2d("UTC", 1994, 6, 30, 23, 59, 60.13599)
@test_approx_eq_eps jd1+jd2  2449534.49999 1e-6

# eraEe00
epsa =  0.4090789763356509900
dpsi = -0.9630909107115582393e-5
ee = eraEe00(2400000.5, 53736.0, epsa, dpsi)
@test_approx_eq_eps ee  -0.8834193235367965479e-5 1e-18

# eraEe00a
ee = eraEe00a(2400000.5, 53736.0)
@test_approx_eq_eps ee  -0.8834192459222588227e-5 1e-18

# eraEe00b
ee = eraEe00b(2400000.5, 53736.0)
@test_approx_eq_eps ee  -0.8835700060003032831e-5 1e-18

# eraEe06a 
ee = eraEe06a(2400000.5, 53736.0)
@test_approx_eq_eps ee  -0.8834195072043790156e-5 1e-15

# eraEform
a, f = eraEform(1)
@test_approx_eq_eps a  6378137.0 1e-10
@test_approx_eq_eps f  0.0033528106647474807 1e-18
a, f = eraEform(2)
@test_approx_eq_eps a  6378137.0 1e-10
@test_approx_eq_eps f  0.0033528106811823189 1e-18
a, f = eraEform(3)
@test_approx_eq_eps a  6378135.0 1e-10
@test_approx_eq_eps f  0.0033527794541675049 1e-18

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

# eraEpv00
pvh,pvb = eraEpv00(2400000.5, 53411.52501161)
@test_approx_eq_eps pvh[1]  -0.7757238809297706813 1e-14
@test_approx_eq_eps pvh[2]  0.5598052241363340596 1e-14
@test_approx_eq_eps pvh[3]  0.2426998466481686993 1e-14
@test_approx_eq_eps pvh[4]  -0.1091891824147313846e-1 1e-15
@test_approx_eq_eps pvh[5]  -0.1247187268440845008e-1 1e-15
@test_approx_eq_eps pvh[6]  -0.5407569418065039061e-2 1e-15
@test_approx_eq_eps pvb[1]  -0.7714104440491111971 1e-14
@test_approx_eq_eps pvb[2]  0.5598412061824171323 1e-14
@test_approx_eq_eps pvb[3]  0.2425996277722452400 1e-14
@test_approx_eq_eps pvb[4]  -0.1091874268116823295e-1 1e-15
@test_approx_eq_eps pvb[5]  -0.1246525461732861538e-1 1e-15
@test_approx_eq_eps pvb[6]  -0.5404773180966231279e-2 1e-15

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

# eraFw2m
gamb = -0.2243387670997992368e-5
phib =  0.4091014602391312982
psi  = -0.9501954178013015092e-3
eps  =  0.4091014316587367472
r = eraFw2m(gamb, phib, psi, eps)
@test_approx_eq_eps r[1]  0.9999995505176007047 1e-12
@test_approx_eq_eps r[2]  0.8695404617348192957e-3 1e-12
@test_approx_eq_eps r[3]  0.3779735201865582571e-3 1e-12
@test_approx_eq_eps r[4]  -0.8695404723772016038e-3 1e-12
@test_approx_eq_eps r[5]  0.9999996219496027161 1e-12
@test_approx_eq_eps r[6]  -0.1361752496887100026e-6 1e-12
@test_approx_eq_eps r[7]  -0.3779734957034082790e-3 1e-12
@test_approx_eq_eps r[8]  -0.1924880848087615651e-6 1e-12
@test_approx_eq_eps r[9]  0.9999999285679971958 1e-12

# eraGmst00
g = eraGmst00(2400000.5, 53736.0, 2400000.5, 53736.0)
@test_approx_eq_eps g  1.754174972210740592 1e-14

# eraGmst06
g = eraGmst06(2400000.5, 53736.0, 2400000.5, 53736.0)
@test_approx_eq_eps g  1.754174971870091203 1e-14

# eraGmst82
g = eraGmst82(2400000.5, 53736.0)
@test_approx_eq_eps g  1.754174981860675096 1e-14

# eraGst00a
g = eraGst00a(2400000.5, 53736.0, 2400000.5, 53736.0)
@test_approx_eq_eps g  1.754166138018281369 1e-14

# eraGst00b
g = eraGst00b(2400000.5, 53736.0)
@test_approx_eq_eps g  1.754166136510680589 1e-14

# eraGst06a
g = eraGst06a(2400000.5, 53736.0, 2400000.5, 53736.0)
@test_approx_eq_eps g  1.754166137675019159 1e-14

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

# eraLdn
sc = [-0.763276255, -0.608633767, -0.216735543]
ob = [-0.974170437, -0.2115201, -0.0917583114]
pv1 = [-7.81014427,-5.60956681,-1.98079819,
       0.0030723249,-0.00406995477,-0.00181335842]
pv2 = [0.738098796, 4.63658692,1.9693136,
       -0.00755816922, 0.00126913722, 0.000727999001]
pv3 = [-0.000712174377, -0.00230478303, -0.00105865966,
       6.29235213e-6, -3.30888387e-7, -2.96486623e-7]
b1 = eraLDBODY(0.00028574, 3e-10, pv1)
b2 = eraLDBODY(0.00095435, 3e-9, pv2)
b3 = eraLDBODY(1.0, 6e-6, pv3)
l = [b1, b2, b3]
sn = eraLdn(l, ob, sc)
@test_approx_eq_eps sn[1]  -0.7632762579693333866 1e-12
@test_approx_eq_eps sn[2]  -0.6086337636093002660 1e-12
@test_approx_eq_eps sn[3]  -0.2167355420646328159 1e-12

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

# eraPap
a = [1.,0.1,0.2]
b = [-3.,1e-3,0.2]
theta = eraPap(a,b)
@test_approx_eq_eps theta  0.3671514267841113674 1e-12

# eraPas
p = eraPas(1.0,0.1,0.2,-1.0)
@test_approx_eq_eps p  -2.724544922932270424 1e-12

# eraPdp
ab =eraPdp([2.,2.,3.],[1.,3.,4.])
@test_approx_eq_eps ab  20 1e-12

# eraPlan94
pv = eraPlan94(2400000.5, -320000., 3)
@test_approx_eq_eps pv[1]  0.9308038666832975759 1e-11
@test_approx_eq_eps pv[2]  0.3258319040261346000 1e-11
@test_approx_eq_eps pv[3]  0.1422794544481140560 1e-11
@test_approx_eq_eps pv[4]  -0.6429458958255170006e-2 1e-11
@test_approx_eq_eps pv[5]  0.1468570657704237764e-1 1e-11
@test_approx_eq_eps pv[6]  0.6406996426270981189e-2 1e-11

pv = eraPlan94(2400000.5, 43999.9, 1)
@test_approx_eq_eps pv[1]  0.2945293959257430832 1e-11
@test_approx_eq_eps pv[2]  -0.2452204176601049596 1e-11
@test_approx_eq_eps pv[3]  -0.1615427700571978153 1e-11
@test_approx_eq_eps pv[4]  0.1413867871404614441e-1 1e-11
@test_approx_eq_eps pv[5]  0.1946548301104706582e-1 1e-11
@test_approx_eq_eps pv[6]  0.8929809783898904786e-2 1e-11

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

# eraPn00
dpsi = -0.9632552291149335877e-5
deps =  0.4063197106621141414e-4
epsa, rb, rp, rbp, rn, rbpn = eraPn00(2400000.5, 53736.0, dpsi, deps)
@test_approx_eq_eps epsa  0.4090791789404229916 1e-12
@test_approx_eq_eps rb[1]  0.9999999999999942498 1e-12
@test_approx_eq_eps rb[2]  -0.7078279744199196626e-7 1e-18
@test_approx_eq_eps rb[3]  0.8056217146976134152e-7 1e-18
@test_approx_eq_eps rb[4]  0.7078279477857337206e-7 1e-18
@test_approx_eq_eps rb[5]  0.9999999999999969484 1e-12
@test_approx_eq_eps rb[6]  0.3306041454222136517e-7 1e-18
@test_approx_eq_eps rb[7]  -0.8056217380986972157e-7 1e-18
@test_approx_eq_eps rb[8]  -0.3306040883980552500e-7 1e-18
@test_approx_eq_eps rb[9]  0.9999999999999962084 1e-12
@test_approx_eq_eps rp[1]  0.9999989300532289018 1e-12
@test_approx_eq_eps rp[2]  -0.1341647226791824349e-2 1e-14
@test_approx_eq_eps rp[3]  -0.5829880927190296547e-3 1e-14
@test_approx_eq_eps rp[4]  0.1341647231069759008e-2 1e-14
@test_approx_eq_eps rp[5]  0.9999990999908750433 1e-12
@test_approx_eq_eps rp[6]  -0.3837444441583715468e-6 1e-14
@test_approx_eq_eps rp[7]  0.5829880828740957684e-3 1e-14
@test_approx_eq_eps rp[8]  -0.3984203267708834759e-6 1e-14
@test_approx_eq_eps rp[9]  0.9999998300623538046 1e-12
@test_approx_eq_eps rbp[1]  0.9999989300052243993 1e-12
@test_approx_eq_eps rbp[2]  -0.1341717990239703727e-2 1e-14
@test_approx_eq_eps rbp[3]  -0.5829075749891684053e-3 1e-14
@test_approx_eq_eps rbp[4]  0.1341718013831739992e-2 1e-14
@test_approx_eq_eps rbp[5]  0.9999990998959191343 1e-12
@test_approx_eq_eps rbp[6]  -0.3505759733565421170e-6 1e-14
@test_approx_eq_eps rbp[7]  0.5829075206857717883e-3 1e-14
@test_approx_eq_eps rbp[8]  -0.4315219955198608970e-6 1e-14
@test_approx_eq_eps rbp[9]  0.9999998301093036269 1e-12
@test_approx_eq_eps rn[1]  0.9999999999536069682 1e-12
@test_approx_eq_eps rn[2]  0.8837746144872140812e-5 1e-16
@test_approx_eq_eps rn[3]  0.3831488838252590008e-5 1e-16
@test_approx_eq_eps rn[4]  -0.8837590456633197506e-5 1e-16
@test_approx_eq_eps rn[5]  0.9999999991354692733 1e-12
@test_approx_eq_eps rn[6]  -0.4063198798559573702e-4 1e-16
@test_approx_eq_eps rn[7]  -0.3831847930135328368e-5 1e-16
@test_approx_eq_eps rn[8]  0.4063195412258150427e-4 1e-16
@test_approx_eq_eps rn[9]  0.9999999991671806225 1e-12
@test_approx_eq_eps rbpn[1]  0.9999989440499982806 1e-12
@test_approx_eq_eps rbpn[2]  -0.1332880253640848301e-2 1e-14
@test_approx_eq_eps rbpn[3]  -0.5790760898731087295e-3 1e-14
@test_approx_eq_eps rbpn[4]  0.1332856746979948745e-2 1e-14
@test_approx_eq_eps rbpn[5]  0.9999991109064768883 1e-12
@test_approx_eq_eps rbpn[6]  -0.4097740555723063806e-4 1e-14
@test_approx_eq_eps rbpn[7]  0.5791301929950205000e-3 1e-14
@test_approx_eq_eps rbpn[8]  0.4020553681373702931e-4 1e-14
@test_approx_eq_eps rbpn[9]  0.9999998314958529887 1e-12

# eraPn06
dpsi = -0.9632552291149335877e-5
deps =  0.4063197106621141414e-4
epsa, rb, rp, rbp, rn, rbpn = eraPn06(2400000.5, 53736.0, dpsi, deps)
@test_approx_eq_eps epsa  0.4090789763356509926 1e-12
@test_approx_eq_eps rb[1]  0.9999999999999942497 1e-12
@test_approx_eq_eps rb[2]  -0.7078368960971557145e-7 1e-14
@test_approx_eq_eps rb[3]  0.8056213977613185606e-7 1e-14
@test_approx_eq_eps rb[4]  0.7078368694637674333e-7 1e-14
@test_approx_eq_eps rb[5]  0.9999999999999969484 1e-12
@test_approx_eq_eps rb[6]  0.3305943742989134124e-7 1e-14
@test_approx_eq_eps rb[7]  -0.8056214211620056792e-7 1e-14
@test_approx_eq_eps rb[8]  -0.3305943172740586950e-7 1e-14
@test_approx_eq_eps rb[9]  0.9999999999999962084 1e-12
@test_approx_eq_eps rp[1]  0.9999989300536854831 1e-12
@test_approx_eq_eps rp[2]  -0.1341646886204443795e-2 1e-14
@test_approx_eq_eps rp[3]  -0.5829880933488627759e-3 1e-14
@test_approx_eq_eps rp[4]  0.1341646890569782183e-2 1e-14
@test_approx_eq_eps rp[5]  0.9999990999913319321 1e-12
@test_approx_eq_eps rp[6]  -0.3835944216374477457e-6 1e-14
@test_approx_eq_eps rp[7]  0.5829880833027867368e-3 1e-14
@test_approx_eq_eps rp[8]  -0.3985701514686976112e-6 1e-14
@test_approx_eq_eps rp[9]  0.9999998300623534950 1e-12
@test_approx_eq_eps rbp[1]  0.9999989300056797893 1e-12
@test_approx_eq_eps rbp[2]  -0.1341717650545059598e-2 1e-14
@test_approx_eq_eps rbp[3]  -0.5829075756493728856e-3 1e-14
@test_approx_eq_eps rbp[4]  0.1341717674223918101e-2 1e-14
@test_approx_eq_eps rbp[5]  0.9999990998963748448 1e-12
@test_approx_eq_eps rbp[6]  -0.3504269280170069029e-6 1e-14
@test_approx_eq_eps rbp[7]  0.5829075211461454599e-3 1e-14
@test_approx_eq_eps rbp[8]  -0.4316708436255949093e-6 1e-14
@test_approx_eq_eps rbp[9]  0.9999998301093032943 1e-12
@test_approx_eq_eps rn[1]  0.9999999999536069682 1e-12
@test_approx_eq_eps rn[2]  0.8837746921149881914e-5 1e-14
@test_approx_eq_eps rn[3]  0.3831487047682968703e-5 1e-14
@test_approx_eq_eps rn[4]  -0.8837591232983692340e-5 1e-14
@test_approx_eq_eps rn[5]  0.9999999991354692664 1e-12
@test_approx_eq_eps rn[6]  -0.4063198798558931215e-4 1e-14
@test_approx_eq_eps rn[7]  -0.3831846139597250235e-5 1e-14
@test_approx_eq_eps rn[8]  0.4063195412258792914e-4 1e-14
@test_approx_eq_eps rn[9]  0.9999999991671806293 1e-12
@test_approx_eq_eps rbpn[1]  0.9999989440504506688 1e-12
@test_approx_eq_eps rbpn[2]  -0.1332879913170492655e-2 1e-14
@test_approx_eq_eps rbpn[3]  -0.5790760923225655753e-3 1e-14
@test_approx_eq_eps rbpn[4]  0.1332856406595754748e-2 1e-14
@test_approx_eq_eps rbpn[5]  0.9999991109069366795 1e-12
@test_approx_eq_eps rbpn[6]  -0.4097725651142641812e-4 1e-14
@test_approx_eq_eps rbpn[7]  0.5791301952321296716e-3 1e-14
@test_approx_eq_eps rbpn[8]  0.4020538796195230577e-4 1e-14
@test_approx_eq_eps rbpn[9]  0.9999998314958576778 1e-12

# eraPnm00a
rbpn = eraPnm00a(2400000.5, 50123.9999)
@test_approx_eq_eps rbpn[1]  0.9999995832793134257 1e-12
@test_approx_eq_eps rbpn[2]  0.8372384254137809439e-3 1e-14
@test_approx_eq_eps rbpn[3]  0.3639684306407150645e-3 1e-14
@test_approx_eq_eps rbpn[4]  -0.8372535226570394543e-3 1e-14
@test_approx_eq_eps rbpn[5]  0.9999996486491582471 1e-12
@test_approx_eq_eps rbpn[6]  0.4132915262664072381e-4 1e-14
@test_approx_eq_eps rbpn[7]  -0.3639337004054317729e-3 1e-14
@test_approx_eq_eps rbpn[8]  -0.4163386925461775873e-4 1e-14
@test_approx_eq_eps rbpn[9]  0.9999999329094390695 1e-12

# eraPnm00b
rbpn = eraPnm00b(2400000.5, 50123.9999)
@test_approx_eq_eps rbpn[1]  0.9999995832776208280 1e-12
@test_approx_eq_eps rbpn[2]  0.8372401264429654837e-3 1e-14
@test_approx_eq_eps rbpn[3]  0.3639691681450271771e-3 1e-14
@test_approx_eq_eps rbpn[4]  -0.8372552234147137424e-3 1e-14
@test_approx_eq_eps rbpn[5]  0.9999996486477686123 1e-12
@test_approx_eq_eps rbpn[6]  0.4132832190946052890e-4 1e-14
@test_approx_eq_eps rbpn[7]  -0.3639344385341866407e-3 1e-14
@test_approx_eq_eps rbpn[8]  -0.4163303977421522785e-4 1e-14
@test_approx_eq_eps rbpn[9]  0.9999999329092049734 1e-12

# eraPnm06a
rbpn = eraPnm06a(2400000.5, 50123.9999)
@test_approx_eq_eps rbpn[1]  0.9999995832794205484 1e-12
@test_approx_eq_eps rbpn[2]  0.8372382772630962111e-3 1e-14
@test_approx_eq_eps rbpn[3]  0.3639684771140623099e-3 1e-14
@test_approx_eq_eps rbpn[4]  -0.8372533744743683605e-3 1e-14
@test_approx_eq_eps rbpn[5]  0.9999996486492861646 1e-12
@test_approx_eq_eps rbpn[6]  0.4132905944611019498e-4 1e-14
@test_approx_eq_eps rbpn[7]  -0.3639337469629464969e-3 1e-14
@test_approx_eq_eps rbpn[8]  -0.4163377605910663999e-4 1e-14
@test_approx_eq_eps rbpn[9]  0.9999999329094260057 1e-12

# eraPnm80
rmatpn = eraPnm80(2400000.5, 50123.9999)
@test_approx_eq_eps rmatpn[1]  0.9999995831934611169 1e-12
@test_approx_eq_eps rmatpn[2]  0.8373654045728124011e-3 1e-14
@test_approx_eq_eps rmatpn[3]  0.3639121916933106191e-3 1e-14
@test_approx_eq_eps rmatpn[4]  -0.8373804896118301316e-3 1e-14
@test_approx_eq_eps rmatpn[5]  0.9999996485439674092 1e-12
@test_approx_eq_eps rmatpn[6]  0.4130202510421549752e-4 1e-14
@test_approx_eq_eps rmatpn[7]  -0.3638774789072144473e-3 1e-14
@test_approx_eq_eps rmatpn[8]  -0.4160674085851722359e-4 1e-14
@test_approx_eq_eps rmatpn[9]  0.9999999329310274805 1e-12

# eraPom00
xp =  2.55060238e-7
yp =  1.860359247e-6
sp = -0.1367174580728891460e-10
rpom = eraPom00(xp, yp, sp)
@test_approx_eq_eps rpom[1]  0.9999999999999674721 1e-12
@test_approx_eq_eps rpom[2]  -0.1367174580728846989e-10 1e-16
@test_approx_eq_eps rpom[3]  0.2550602379999972345e-6 1e-16
@test_approx_eq_eps rpom[4]  0.1414624947957029801e-10 1e-16
@test_approx_eq_eps rpom[5]  0.9999999999982695317 1e-12
@test_approx_eq_eps rpom[6]  -0.1860359246998866389e-5 1e-16
@test_approx_eq_eps rpom[7]  -0.2550602379741215021e-6 1e-16
@test_approx_eq_eps rpom[8]  0.1860359247002414021e-5 1e-16
@test_approx_eq_eps rpom[9]  0.9999999999982370039 1e-12

# eraRx
phi = 0.3456789
r = [[2.0,3.0,2.0],
     [3.0,2.0,3.0],
     [3.0,4.0,5.0]]
r = eraRx(phi, r)
@test_approx_eq_eps r[1]  2.0 0.0
@test_approx_eq_eps r[2]  3.0 0.0
@test_approx_eq_eps r[3]  2.0 0.0
@test_approx_eq_eps r[4]  3.839043388235612460 1e-12
@test_approx_eq_eps r[5]  3.237033249594111899 1e-12
@test_approx_eq_eps r[6]  4.516714379005982719 1e-12
@test_approx_eq_eps r[7]  1.806030415924501684 1e-12
@test_approx_eq_eps r[8]  3.085711545336372503 1e-12
@test_approx_eq_eps r[9]  3.687721683977873065 1e-12

# eraRxp
r = [[2.0,3.0,2.0],
     [3.0,2.0,3.0],
     [3.0,4.0,5.0]]
p = [0.2,1.5,0.1]
rp = eraRxp(r, p)
@test_approx_eq_eps rp[1]  5.1 1e-12
@test_approx_eq_eps rp[2]  3.9 1e-12
@test_approx_eq_eps rp[3]  7.1 1e-12

# eraRxpv
r = [[2.0,3.0,2.0],
     [3.0,2.0,3.0],
     [3.0,4.0,5.0]]
pv = [[0.2,1.5,0.1],
      [1.5,0.2,0.1]]
rpv = eraRxpv(r, pv)
@test_approx_eq_eps rpv[1]  5.1 1e-12
@test_approx_eq_eps rpv[4]  3.8 1e-12
@test_approx_eq_eps rpv[2]  3.9 1e-12
@test_approx_eq_eps rpv[5]  5.2 1e-12
@test_approx_eq_eps rpv[3]  7.1 1e-12
@test_approx_eq_eps rpv[6]  5.8 1e-12

# eraRxr
a = [[2.0,3.0,2.0],
     [3.0,2.0,3.0],
     [3.0,4.0,5.0]]
b = [[1.0,2.0,2.0],
     [4.0,1.0,1.0],
     [3.0,0.0,1.0]]
atb = eraRxr(a, b)
@test_approx_eq_eps atb[1]  20.0 1e-12
@test_approx_eq_eps atb[2]   7.0 1e-12
@test_approx_eq_eps atb[3]   9.0 1e-12
@test_approx_eq_eps atb[4]  20.0 1e-12
@test_approx_eq_eps atb[5]   8.0 1e-12
@test_approx_eq_eps atb[6]  11.0 1e-12
@test_approx_eq_eps atb[7]  34.0 1e-12
@test_approx_eq_eps atb[8]  10.0 1e-12
@test_approx_eq_eps atb[9]  15.0 1e-12

# eraRy
theta = 0.3456789
r = [[2.0,3.0,2.0],
     [3.0,2.0,3.0],
     [3.0,4.0,5.0]]
r = eraRy(theta, r)
@test_approx_eq_eps r[1]  0.8651847818978159930 1e-12
@test_approx_eq_eps r[2]  1.467194920539316554 1e-12
@test_approx_eq_eps r[3]  0.1875137911274457342 1e-12
@test_approx_eq_eps r[4]  3 1e-12
@test_approx_eq_eps r[5]  2 1e-12
@test_approx_eq_eps r[6]  3 1e-12
@test_approx_eq_eps r[7]  3.500207892850427330 1e-12
@test_approx_eq_eps r[8]  4.779889022262298150 1e-12
@test_approx_eq_eps r[9]  5.381899160903798712 1e-12

# eraRz
psi = 0.3456789
r = [[2.0,3.0,2.0],
     [3.0,2.0,3.0],
     [3.0,4.0,5.0]]
r = eraRz(psi, r)
@test_approx_eq_eps r[1]  2.898197754208926769 1e-12
@test_approx_eq_eps r[2]  3.500207892850427330 1e-12
@test_approx_eq_eps r[3]  2.898197754208926769 1e-12
@test_approx_eq_eps r[4]  2.144865911309686813 1e-12
@test_approx_eq_eps r[5]  0.865184781897815993 1e-12
@test_approx_eq_eps r[6]  2.144865911309686813 1e-12
@test_approx_eq_eps r[7]  3.0 1e-12
@test_approx_eq_eps r[8]  4.0 1e-12
@test_approx_eq_eps r[9]  5.0 1e-12

# eraS00
x = 0.5791308486706011000e-3
y = 0.4020579816732961219e-4
s = eraS00(2400000.5, 53736.0, x, y)
@test_approx_eq_eps s  -0.1220036263270905693e-7 1e-18

# eraS00a
s = eraS00a(2400000.5, 52541.0)
@test_approx_eq_eps s  -0.1340684448919163584e-7 1e-18

# eraS00b
s = eraS00b(2400000.5, 52541.0)
@test_approx_eq_eps s  -0.1340695782951026584e-7 1e-18

# eraS06
x = 0.5791308486706011000e-3
y = 0.4020579816732961219e-4
s = eraS06(2400000.5, 53736.0, x, y)
@test_approx_eq_eps s  -0.1220032213076463117e-7 1e-18

# eraS06a
s = eraS06a(2400000.5, 52541.0)
@test_approx_eq_eps s  -0.1340680437291812383e-7 1e-18

# eraSepp
a = [1.,0.1,0.2]
b = [-3.,1e-3,0.2]
s = eraSepp(a,b)
@test_approx_eq_eps s  2.860391919024660768 1e-12

# eraSeps
s = eraSeps(1.,.1,.2,-3.)
@test_approx_eq_eps s  2.346722016996998842 1e-14

# eraSp00
s = eraSp00(2400000.5, 52541.0)
@test_approx_eq_eps s  -0.6216698469981019309e-11 1e-12

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

# eraTr
r = [[2.0,3.0,2.0],
     [3.0,2.0,3.0],
     [3.0,4.0,5.0]]
rt = eraTr(r)
@test_approx_eq_eps rt[1]  2.0 0.0
@test_approx_eq_eps rt[2]  3.0 0.0
@test_approx_eq_eps rt[3]  3.0 0.0
@test_approx_eq_eps rt[4]  3.0 0.0
@test_approx_eq_eps rt[5]  2.0 0.0
@test_approx_eq_eps rt[6]  4.0 0.0
@test_approx_eq_eps rt[7]  2.0 0.0
@test_approx_eq_eps rt[8]  3.0 0.0
@test_approx_eq_eps rt[9]  5.0 0.0

# eraTrxp
r = [[2.0,3.0,2.0],
     [3.0,2.0,3.0],
     [3.0,4.0,5.0]]
p = [0.2,1.5,0.1]
trp = eraTrxp(r, p)
@test_approx_eq_eps trp[1]  5.2 1e-12
@test_approx_eq_eps trp[2]  4.0 1e-12
@test_approx_eq_eps trp[3]  5.4 1e-12

# eraTrxpv
r = [[2.0,3.0,2.0],
     [3.0,2.0,3.0],
     [3.0,4.0,5.0]]
pv = [[0.2,1.5,0.1],
      [1.5,0.2,0.1]]
trpv = eraTrxpv(r, pv)
@test_approx_eq_eps trpv[1]  5.2 1e-12
@test_approx_eq_eps trpv[2]  4.0 1e-12
@test_approx_eq_eps trpv[3]  5.4 1e-12
@test_approx_eq_eps trpv[4]  3.9 1e-12
@test_approx_eq_eps trpv[5]  5.3 1e-12
@test_approx_eq_eps trpv[6]  4.1 1e-12

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

# eraXy06
x, y = eraXy06(2400000.5, 53736.0)
@test_approx_eq_eps x  0.5791308486706010975e-3 1e-16
@test_approx_eq_eps y  0.4020579816732958141e-4 1e-17

# eraXys00a
x, y, s = eraXys00a(2400000.5, 53736.0)
@test_approx_eq_eps x  0.5791308472168152904e-3 1e-16
@test_approx_eq_eps y  0.4020595661591500259e-4 1e-17
@test_approx_eq_eps s  -0.1220040848471549623e-7 1e-20

# eraXys00b
x, y, s = eraXys00b(2400000.5, 53736.0)
@test_approx_eq_eps x  0.5791301929950208873e-3 1e-16
@test_approx_eq_eps y  0.4020553681373720832e-4 1e-16
@test_approx_eq_eps s  -0.1220027377285083189e-7 1e-19

# eraXys06a
x, y, s = eraXys06a(2400000.5, 53736.0)
@test_approx_eq_eps x  0.5791308482835292617e-3 1e-16
@test_approx_eq_eps y  0.4020580099454020310e-4 1e-15
@test_approx_eq_eps s  -0.1220032294164579896e-7 1e-19
