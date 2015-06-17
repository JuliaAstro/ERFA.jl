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

# eraAb
pnat = [-0.76321968546737951,-0.60869453983060384,-0.21676408580639883]
v = [2.1044018893653786e-5,-8.9108923304429319e-5,-3.8633714797716569e-5]
s = 0.99980921395708788
bm1 = 0.99999999506209258
ppr = eraAb(pnat, v, s, bm1)
@test_approx_eq_eps ppr[1]  -0.7631631094219556269 1e-12
@test_approx_eq_eps ppr[2]  -0.6087553082505590832 1e-12
@test_approx_eq_eps ppr[3]  -0.2167926269368471279 1e-12

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

# eraApcg
date1 = 2456165.5
date2 = 0.401182685
ebpv=[[0.901310875,-0.417402664,-0.180982288];
      [ 0.00742727954,0.0140507459,0.00609045792]]
ehp=[0.903358544,-0.415395237,-0.180084014]
astrom = eraApcg(date1,date2,ebpv,ehp)
@test_approx_eq_eps astrom.pmt  12.65133794027378508 1e-11
@test_approx_eq_eps astrom.eb[1]  0.901310875 1e-12
@test_approx_eq_eps astrom.eb[2]  -0.417402664 1e-12
@test_approx_eq_eps astrom.eb[3]  -0.180982288 1e-12
@test_approx_eq_eps astrom.eh[1]  0.8940025429324143045 1e-12
@test_approx_eq_eps astrom.eh[2]  -0.4110930268679817955 1e-12
@test_approx_eq_eps astrom.eh[3]  -0.1782189004872870264 1e-12
@test_approx_eq_eps astrom.em  1.010465295811013146 1e-12
@test_approx_eq_eps astrom.v[1]  0.4289638897813379954e-4 1e-16
@test_approx_eq_eps astrom.v[2]  0.8115034021720941898e-4 1e-16
@test_approx_eq_eps astrom.v[3]  0.3517555123437237778e-4 1e-16
@test_approx_eq_eps astrom.bm1  0.9999999951686013336 1e-15
@test_approx_eq_eps astrom.bpn[1]  1.0 1e-10
@test_approx_eq_eps astrom.bpn[4]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[7]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[2]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[5]  1.0 1e-10
@test_approx_eq_eps astrom.bpn[8]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[3]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[6]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[9]  1.0 1e-10

# eraApcg13
date1 = 2456165.5
date2 = 0.401182685
astrom = eraApcg13(date1,date2)
@test_approx_eq_eps astrom.pmt  12.65133794027378508 1e-12
@test_approx_eq_eps astrom.eb[1]  0.9013108747340644755 1e-12
@test_approx_eq_eps astrom.eb[2]  -0.4174026640406119957 1e-12
@test_approx_eq_eps astrom.eb[3]  -0.1809822877867817771 1e-12
@test_approx_eq_eps astrom.eh[1]  0.8940025429255499549 1e-12
@test_approx_eq_eps astrom.eh[2]  -0.4110930268331896318 1e-12
@test_approx_eq_eps astrom.eh[3]  -0.1782189006019749850 1e-12
@test_approx_eq_eps astrom.em  1.010465295964664178 1e-12
@test_approx_eq_eps astrom.v[1]  0.4289638897157027528e-4 1e-16
@test_approx_eq_eps astrom.v[2]  0.8115034002544663526e-4 1e-16
@test_approx_eq_eps astrom.v[3]  0.3517555122593144633e-4 1e-16
@test_approx_eq_eps astrom.bm1  0.9999999951686013498 1e-12
@test_approx_eq_eps astrom.bpn[1]  1.0 1e-10
@test_approx_eq_eps astrom.bpn[4]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[7]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[2]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[5]  1.0 1e-10
@test_approx_eq_eps astrom.bpn[8]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[3]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[6]  0.0 1e-10
@test_approx_eq_eps astrom.bpn[9]  1.0 1e-10

# eraApci
date1 = 2456165.5
date2 = 0.401182685
ebpv = [[0.901310875,-0.417402664,-0.180982288];
        [0.00742727954,0.0140507459,0.00609045792]]
ehp = [0.903358544,-0.415395237,-0.180084014]
x =  0.0013122272
y = -2.92808623e-5
s =  3.05749468e-8
astrom = eraApci(date1, date2, ebpv, ehp, x, y, s)
@test_approx_eq_eps astrom.pmt  12.65133794027378508 1e-11
@test_approx_eq_eps astrom.eb[1]  0.901310875 1e-12
@test_approx_eq_eps astrom.eb[2]  -0.417402664 1e-12
@test_approx_eq_eps astrom.eb[3]  -0.180982288 1e-12
@test_approx_eq_eps astrom.eh[1]  0.8940025429324143045 1e-12
@test_approx_eq_eps astrom.eh[2]  -0.4110930268679817955 1e-12
@test_approx_eq_eps astrom.eh[3]  -0.1782189004872870264 1e-12
@test_approx_eq_eps astrom.em  1.010465295811013146 1e-12
@test_approx_eq_eps astrom.v[1]  0.4289638897813379954e-4 1e-16
@test_approx_eq_eps astrom.v[2]  0.8115034021720941898e-4 1e-16
@test_approx_eq_eps astrom.v[3]  0.3517555123437237778e-4 1e-16
@test_approx_eq_eps astrom.bm1  0.9999999951686013336 1e-12
@test_approx_eq_eps astrom.bpn[1]  0.9999991390295159156 1e-12
@test_approx_eq_eps astrom.bpn[4]  0.4978650072505016932e-7 1e-12
@test_approx_eq_eps astrom.bpn[7]  0.1312227200000000000e-2 1e-12
@test_approx_eq_eps astrom.bpn[2]  -0.1136336653771609630e-7 1e-12
@test_approx_eq_eps astrom.bpn[5]  0.9999999995713154868 1e-12
@test_approx_eq_eps astrom.bpn[8]  -0.2928086230000000000e-4 1e-12
@test_approx_eq_eps astrom.bpn[3]  -0.1312227200895260194e-2 1e-12
@test_approx_eq_eps astrom.bpn[6]  0.2928082217872315680e-4 1e-12
@test_approx_eq_eps astrom.bpn[9]  0.9999991386008323373 1e-12

# eraApci13
date1 = 2456165.5
date2 = 0.401182685
astrom, eo = eraApci13(date1, date2)
@test_approx_eq_eps astrom.pmt  12.65133794027378508 1e-11
@test_approx_eq_eps astrom.eb[1]  0.9013108747340644755 1e-12
@test_approx_eq_eps astrom.eb[2]  -0.4174026640406119957 1e-12
@test_approx_eq_eps astrom.eb[3]  -0.1809822877867817771 1e-12
@test_approx_eq_eps astrom.eh[1]  0.8940025429255499549 1e-12
@test_approx_eq_eps astrom.eh[2]  -0.4110930268331896318 1e-12
@test_approx_eq_eps astrom.eh[3]  -0.1782189006019749850 1e-12
@test_approx_eq_eps astrom.em  1.010465295964664178 1e-12
@test_approx_eq_eps astrom.v[1]  0.4289638897157027528e-4 1e-16
@test_approx_eq_eps astrom.v[2]  0.8115034002544663526e-4 1e-16
@test_approx_eq_eps astrom.v[3]  0.3517555122593144633e-4 1e-16
@test_approx_eq_eps astrom.bm1  0.9999999951686013498 1e-12
@test_approx_eq_eps astrom.bpn[1]  0.9999992060376761710 1e-12
@test_approx_eq_eps astrom.bpn[4]  0.4124244860106037157e-7 1e-12
@test_approx_eq_eps astrom.bpn[7]  0.1260128571051709670e-2 1e-12
@test_approx_eq_eps astrom.bpn[2]  -0.1282291987222130690e-7 1e-12
@test_approx_eq_eps astrom.bpn[5]  0.9999999997456835325 1e-12
@test_approx_eq_eps astrom.bpn[8]  -0.2255288829420524935e-4 1e-12
@test_approx_eq_eps astrom.bpn[3]  -0.1260128571661374559e-2 1e-12
@test_approx_eq_eps astrom.bpn[6]  0.2255285422953395494e-4 1e-12
@test_approx_eq_eps astrom.bpn[9]  0.9999992057833604343 1e-12
@test_approx_eq_eps eo  -0.2900618712657375647e-2 1e-12

# eraApco
date1 = 2456384.5
date2 = 0.970031644
ebpv = [[-0.974170438,-0.211520082,-0.0917583024];
        [0.00364365824,-0.0154287319,-0.00668922024]]
ehp = [-0.973458265,-0.209215307,-0.0906996477]
x = 0.0013122272
y = -2.92808623e-5
s = 3.05749468e-8
theta = 3.14540971
elong = -0.527800806
phi = -1.2345856
hm = 2738.0
xp = 2.47230737e-7
yp = 1.82640464e-6
sp = -3.01974337e-11
refa = 0.000201418779
refb = -2.36140831e-7
astrom = eraApco(date1, date2, ebpv, ehp, x, y, s,theta, elong, phi, hm, xp, yp, sp,refa, refb)
@test_approx_eq_eps astrom.pmt  13.25248468622587269 1e-11
@test_approx_eq_eps astrom.eb[1]  -0.9741827110630897003 1e-12
@test_approx_eq_eps astrom.eb[2]  -0.2115130190135014340 1e-12
@test_approx_eq_eps astrom.eb[3]  -0.09179840186968295686 1e-12
@test_approx_eq_eps astrom.eh[1]  -0.9736425571689670428 1e-12
@test_approx_eq_eps astrom.eh[2]  -0.2092452125848862201 1e-12
@test_approx_eq_eps astrom.eh[3]  -0.09075578152261439954 1e-12
@test_approx_eq_eps astrom.em  0.9998233241710617934 1e-12
@test_approx_eq_eps astrom.v[1]  0.2078704985147609823e-4 1e-16
@test_approx_eq_eps astrom.v[2]  -0.8955360074407552709e-4 1e-16
@test_approx_eq_eps astrom.v[3]  -0.3863338980073114703e-4 1e-16
@test_approx_eq_eps astrom.bm1  0.9999999950277561600 1e-12
@test_approx_eq_eps astrom.bpn[1]  0.9999991390295159156 1e-12
@test_approx_eq_eps astrom.bpn[4]  0.4978650072505016932e-7 1e-12
@test_approx_eq_eps astrom.bpn[7]  0.1312227200000000000e-2 1e-12
@test_approx_eq_eps astrom.bpn[2]  -0.1136336653771609630e-7 1e-12
@test_approx_eq_eps astrom.bpn[5]  0.9999999995713154868 1e-12
@test_approx_eq_eps astrom.bpn[8]  -0.2928086230000000000e-4 1e-12
@test_approx_eq_eps astrom.bpn[3]  -0.1312227200895260194e-2 1e-12
@test_approx_eq_eps astrom.bpn[6]  0.2928082217872315680e-4 1e-12
@test_approx_eq_eps astrom.bpn[9]  0.9999991386008323373 1e-12
@test_approx_eq_eps astrom.along  -0.5278008060301974337 1e-12
@test_approx_eq_eps astrom.xpl  0.1133427418174939329e-5 1e-17
@test_approx_eq_eps astrom.ypl  0.1453347595745898629e-5 1e-17
@test_approx_eq_eps astrom.sphi  -0.9440115679003211329 1e-12
@test_approx_eq_eps astrom.cphi  0.3299123514971474711 1e-12
@test_approx_eq_eps astrom.diurab  0 1e-10
@test_approx_eq_eps astrom.eral  2.617608903969802566 1e-12
@test_approx_eq_eps astrom.refa  0.2014187790000000000e-3 1e-15
@test_approx_eq_eps astrom.refb  -0.2361408310000000000e-6 1e-18

# eraApco13
utc1 = 2456384.5
utc2 = 0.969254051
dut1 = 0.1550675
elong = -0.527800806
phi = -1.2345856
hm = 2738.0
p = 2.47230737e-7
xp = 2.47230737e-7
yp = 1.82640464e-6
phpa = 731.0
tc = 12.8
rh = 0.59
wl = 0.55
astrom, eo = eraApco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
@test_approx_eq_eps astrom.pmt  13.25248468622475727 1e-11
@test_approx_eq_eps astrom.eb[1]  -0.9741827107321449445 1e-12
@test_approx_eq_eps astrom.eb[2]  -0.2115130190489386190 1e-12
@test_approx_eq_eps astrom.eb[3]  -0.09179840189515518726 1e-12
@test_approx_eq_eps astrom.eh[1]  -0.9736425572586866640 1e-12
@test_approx_eq_eps astrom.eh[2]  -0.2092452121602867431 1e-12
@test_approx_eq_eps astrom.eh[3]  -0.09075578153903832650 1e-12
@test_approx_eq_eps astrom.em  0.9998233240914558422 1e-12
@test_approx_eq_eps astrom.v[1]  0.2078704986751370303e-4 1e-16
@test_approx_eq_eps astrom.v[2]  -0.8955360100494469232e-4 1e-16
@test_approx_eq_eps astrom.v[3]  -0.3863338978840051024e-4 1e-16
@test_approx_eq_eps astrom.bm1  0.9999999950277561368 1e-12
@test_approx_eq_eps astrom.bpn[1]  0.9999991390295147999 1e-12
@test_approx_eq_eps astrom.bpn[4]  0.4978650075315529277e-7 1e-12
@test_approx_eq_eps astrom.bpn[7]  0.001312227200850293372 1e-12
@test_approx_eq_eps astrom.bpn[2]  -0.1136336652812486604e-7 1e-12
@test_approx_eq_eps astrom.bpn[5]  0.9999999995713154865 1e-12
@test_approx_eq_eps astrom.bpn[8]  -0.2928086230975367296e-4 1e-12
@test_approx_eq_eps astrom.bpn[3]  -0.001312227201745553566 1e-12
@test_approx_eq_eps astrom.bpn[6]  0.2928082218847679162e-4 1e-12
@test_approx_eq_eps astrom.bpn[9]  0.9999991386008312212 1e-12
@test_approx_eq_eps astrom.along  -0.5278008060301974337 1e-12
@test_approx_eq_eps astrom.xpl  0.1133427418174939329e-5 1e-17
@test_approx_eq_eps astrom.ypl  0.1453347595745898629e-5 1e-17
@test_approx_eq_eps astrom.sphi  -0.9440115679003211329 1e-12
@test_approx_eq_eps astrom.cphi  0.3299123514971474711 1e-12
@test_approx_eq_eps astrom.diurab  0 1e-10
@test_approx_eq_eps astrom.eral  2.617608909189066140 1e-12
@test_approx_eq_eps astrom.refa  0.2014187785940396921e-3 1e-15
@test_approx_eq_eps astrom.refb  -0.2361408314943696227e-6 1e-18
@test_approx_eq_eps eo  -0.003020548354802412839 1e-14

# eraApcs
date1 = 2456384.5
date2 = 0.970031644
pv = [[-1836024.09,1056607.72,-5998795.26];
      [-77.0361767,-133.310856,0.0971855934]]
ebpv = [[-0.974170438,-0.211520082,-0.0917583024];
        [0.00364365824,-0.0154287319,-0.00668922024]]
ehp = [-0.973458265,-0.209215307,-0.0906996477]
astrom = eraApcs(date1, date2, pv, ebpv, ehp)
@test_approx_eq_eps astrom.pmt  13.25248468622587269 1e-11
@test_approx_eq_eps astrom.eb[1]  -0.9741827110630456169 1e-12
@test_approx_eq_eps astrom.eb[2]  -0.2115130190136085494 1e-12
@test_approx_eq_eps astrom.eb[3]  -0.09179840186973175487 1e-12
@test_approx_eq_eps astrom.eh[1]  -0.9736425571689386099 1e-12
@test_approx_eq_eps astrom.eh[2]  -0.2092452125849967195 1e-12
@test_approx_eq_eps astrom.eh[3]  -0.09075578152266466572 1e-12
@test_approx_eq_eps astrom.em  0.9998233241710457140 1e-12
@test_approx_eq_eps astrom.v[1]  0.2078704985513566571e-4 1e-16
@test_approx_eq_eps astrom.v[2]  -0.8955360074245006073e-4 1e-16
@test_approx_eq_eps astrom.v[3]  -0.3863338980073572719e-4 1e-16
@test_approx_eq_eps astrom.bm1  0.9999999950277561601 1e-12
@test_approx_eq_eps astrom.bpn[1]  1 1e-10
@test_approx_eq_eps astrom.bpn[4]  0 1e-10
@test_approx_eq_eps astrom.bpn[7]  0 1e-10
@test_approx_eq_eps astrom.bpn[2]  0 1e-10
@test_approx_eq_eps astrom.bpn[5]  1 1e-10
@test_approx_eq_eps astrom.bpn[8]  0 1e-10
@test_approx_eq_eps astrom.bpn[3]  0 1e-10
@test_approx_eq_eps astrom.bpn[6]  0 1e-10
@test_approx_eq_eps astrom.bpn[9]  1 1e-10

# eraApcs13
date1 = 2456165.5
date2 = 0.401182685
pv = [[-6241497.16,401346.896,-1251136.04];
      [-29.264597,-455.021831,0.0266151194]]
astrom = eraApcs13(date1, date2, pv)
@test_approx_eq_eps astrom.pmt  12.65133794027378508 1e-11
@test_approx_eq_eps astrom.eb[1]  0.9012691529023298391 1e-12
@test_approx_eq_eps astrom.eb[2]  -0.4173999812023068781 1e-12
@test_approx_eq_eps astrom.eb[3]  -0.1809906511146821008 1e-12
@test_approx_eq_eps astrom.eh[1]  0.8939939101759726824 1e-12
@test_approx_eq_eps astrom.eh[2]  -0.4111053891734599955 1e-12
@test_approx_eq_eps astrom.eh[3]  -0.1782336880637689334 1e-12
@test_approx_eq_eps astrom.em  1.010428384373318379 1e-12
@test_approx_eq_eps astrom.v[1]  0.4279877278327626511e-4 1e-16
@test_approx_eq_eps astrom.v[2]  0.7963255057040027770e-4 1e-16
@test_approx_eq_eps astrom.v[3]  0.3517564000441374759e-4 1e-16
@test_approx_eq_eps astrom.bm1  0.9999999952947981330 1e-12
@test_approx_eq_eps astrom.bpn[1]  1 1e-10
@test_approx_eq_eps astrom.bpn[4]  0 1e-10
@test_approx_eq_eps astrom.bpn[7]  0 1e-10
@test_approx_eq_eps astrom.bpn[2]  0 1e-10
@test_approx_eq_eps astrom.bpn[5]  1 1e-10
@test_approx_eq_eps astrom.bpn[8]  0 1e-10
@test_approx_eq_eps astrom.bpn[3]  0 1e-10
@test_approx_eq_eps astrom.bpn[6]  0 1e-10
@test_approx_eq_eps astrom.bpn[9]  1 1e-10

# eraAper
theta = 5.678
pmt = 0.
eb = zeros(3)
eh = zeros(3)
em = 0.
v = zeros(3)
bm1 = 0.
bpn = zeros(9)
along = 1.234
phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb = 0., 0., 0., 0., 0., 0., 0., 0., 0.
astrom = eraASTROM(pmt, eb, eh, em, v, bm1, bpn, along,
                   phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb)
astrom = eraAper(theta, astrom)
@test_approx_eq_eps astrom.eral  6.912000000000000000 1e-12

# eraAper13
ut11 = 2456165.5
ut12 = 0.401182685
pmt = 0.
eb = zeros(3)
eh = zeros(3)
em = 0.
v = zeros(3)
bm1 = 0.
bpn = zeros(9)
along = 1.234
phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb = 0., 0., 0., 0., 0., 0., 0., 0., 0.
astrom = eraASTROM(pmt, eb, eh, em, v, bm1, bpn, along,
                   phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb)
astrom = eraAper13(ut11, ut12, astrom)
@test_approx_eq_eps astrom.eral  3.316236661789694933 1e-12

# eraApio
sp = -3.01974337e-11
theta = 3.14540971
elong = -0.527800806
phi = -1.2345856
hm = 2738.0
xp = 2.47230737e-7
yp = 1.82640464e-6
refa = 0.000201418779
refb = -2.36140831e-7
astrom = eraApio(sp, theta, elong, phi, hm, xp, yp, refa, refb)
@test_approx_eq_eps astrom.along  -0.5278008060301974337 1e-12
@test_approx_eq_eps astrom.xpl  0.1133427418174939329e-5 1e-17
@test_approx_eq_eps astrom.ypl  0.1453347595745898629e-5 1e-17
@test_approx_eq_eps astrom.sphi  -0.9440115679003211329 1e-12
@test_approx_eq_eps astrom.cphi  0.3299123514971474711 1e-12
@test_approx_eq_eps astrom.diurab  0.5135843661699913529e-6 1e-12
@test_approx_eq_eps astrom.eral  2.617608903969802566 1e-12
@test_approx_eq_eps astrom.refa  0.2014187790000000000e-3 1e-15
@test_approx_eq_eps astrom.refb  -0.2361408310000000000e-6 1e-18

# eraApio13
utc1 = 2456384.5
utc2 = 0.969254051
dut1 = 0.1550675
elong = -0.527800806
phi = -1.2345856
hm = 2738.0
xp = 2.47230737e-7
yp = 1.82640464e-6
phpa = 731.0
tc = 12.8
rh = 0.59
wl = 0.55
astrom = eraApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
@test_approx_eq_eps astrom.along  -0.5278008060301974337 1e-12
@test_approx_eq_eps astrom.xpl  0.1133427418174939329e-5 1e-17
@test_approx_eq_eps astrom.ypl  0.1453347595745898629e-5 1e-17
@test_approx_eq_eps astrom.sphi  -0.9440115679003211329 1e-12
@test_approx_eq_eps astrom.cphi  0.3299123514971474711 1e-12
@test_approx_eq_eps astrom.diurab  0.5135843661699913529e-6 1e-12
@test_approx_eq_eps astrom.eral  2.617608909189066140 1e-12
@test_approx_eq_eps astrom.refa  0.2014187785940396921e-3 1e-15
@test_approx_eq_eps astrom.refb  -0.2361408314943696227e-6 1e-18

# eraAtci13
rc = 2.71
dc = 0.174
pr = 1e-5
pd = 5e-6
px = 0.1
rv = 55.0
date1 = 2456165.5
date2 = 0.401182685
ri, di, eo = eraAtci13(rc, dc, pr, pd, px, rv, date1, date2)
@test_approx_eq_eps ri  2.710121572969038991 1e-12
@test_approx_eq_eps di  0.1729371367218230438 1e-12
@test_approx_eq_eps eo  -0.002900618712657375647 1e-14

# eraAtciq
date1 = 2456165.5
date2 = 0.401182685
astrom, eo = eraApci13(date1, date2)
rc = 2.71
dc = 0.174
pr = 1e-5
pd = 5e-6
px = 0.1
rv = 55.0
ri, di = eraAtciq(rc, dc, pr, pd, px, rv, astrom)
@test_approx_eq_eps ri  2.710121572969038991 1e-12
@test_approx_eq_eps di  0.1729371367218230438 1e-12

# eraAtciqn
date1 = 2456165.5
date2 = 0.401182685
astrom, eo = eraApci13(date1, date2)
rc = 2.71
dc = 0.174
pr = 1e-5
pd = 5e-6
px = 0.1
rv = 55.0
b1 = eraLDBODY(0.00028574,3e-10,
               [[-7.81014427,-5.60956681,-1.98079819];
                [0.0030723249,-0.00406995477,-0.00181335842]])
b2 = eraLDBODY(0.00095435,3e-9,
                 [[0.738098796, 4.63658692,1.9693136];
                  [-0.00755816922, 0.00126913722, 0.000727999001]])
b3 = eraLDBODY(1.0,6e-6,
                 [[-0.000712174377, -0.00230478303, -0.00105865966];
                  [6.29235213e-6, -3.30888387e-7, -2.96486623e-7]])
b = [b1; b2; b3]
ri, di = eraAtciqn(rc, dc, pr, pd, px, rv, astrom, b)
@test_approx_eq_eps ri  2.710122008105325582 1e-12
@test_approx_eq_eps di  0.1729371916491459122 1e-12

# eraAtciqz
date1 = 2456165.5
date2 = 0.401182685
astrom, eo = eraApci13(date1, date2)
rc = 2.71
dc = 0.174
ri, di = eraAtciqz(rc, dc, astrom)
@test_approx_eq_eps ri  2.709994899247599271 1e-12
@test_approx_eq_eps di  0.1728740720983623469 1e-12

# eraAtco13
rc = 2.71
dc = 0.174
pr = 1e-5
pd = 5e-6
px = 0.1
rv = 55.0
utc1 = 2456384.5
utc2 = 0.969254051
dut1 = 0.1550675
elong = -0.527800806
phi = -1.2345856
hm = 2738.0
xp = 2.47230737e-7
yp = 1.82640464e-6
phpa = 731.0
tc = 12.8
rh = 0.59
wl = 0.55
aob, zob, hob, dob, rob, eo = eraAtco13(rc, dc, pr, pd, px, rv,
                                        utc1, utc2, dut1, elong, phi, hm, xp, yp,
                                        phpa, tc, rh, wl)
@test_approx_eq_eps aob  0.09251774485358230653 1e-12
@test_approx_eq_eps zob  1.407661405256767021 1e-12
@test_approx_eq_eps hob  -0.09265154431403157925 1e-12
@test_approx_eq_eps dob  0.1716626560075591655 1e-12
@test_approx_eq_eps rob  2.710260453503097719 1e-12
@test_approx_eq_eps eo  -0.003020548354802412839 1e-14

# eraAtic13
ri = 2.710121572969038991
di = 0.1729371367218230438
date1 = 2456165.5
date2 = 0.401182685
rc, dc, eo = eraAtic13(ri, di, date1, date2)
@test_approx_eq_eps rc  2.710126504531374930 1e-12
@test_approx_eq_eps dc  0.1740632537628342320 1e-12
@test_approx_eq_eps eo  -0.002900618712657375647 1e-14

# eraAticq
date1 = 2456165.5
date2 = 0.401182685
rc, dc = eraAticq(ri, di, astrom)
astrom, eo = eraApci13(date1, date2)
@test_approx_eq_eps rc  2.710126504531374930 1e-12
@test_approx_eq_eps dc  0.1740632537628342320 1e-12

# eraAticqn
date1 = 2456165.5
date2 = 0.401182685
astrom, eo = eraApci13(date1, date2)
ri = 2.709994899247599271
di = 0.1728740720983623469
b1 = eraLDBODY(0.00028574,3e-10,
               [[-7.81014427,-5.60956681,-1.98079819];
                [0.0030723249,-0.00406995477,-0.00181335842]])
b2 = eraLDBODY(0.00095435,3e-9,
               [[0.738098796, 4.63658692,1.9693136];
                [-0.00755816922, 0.00126913722, 0.000727999001]])
b3 = eraLDBODY(1.0,6e-6,
               [[-0.000712174377, -0.00230478303, -0.00105865966];
                [6.29235213e-6, -3.30888387e-7, -2.96486623e-7]])
b = [b1; b2; b3]
rc, dc = eraAticqn(ri, di, astrom, b)
@test_approx_eq_eps rc  2.709999575032685412 1e-12
@test_approx_eq_eps dc  0.1739999656317778034 1e-12

# eraAtio13
ri = 2.710121572969038991
di = 0.1729371367218230438
utc1 = 2456384.5
utc2 = 0.969254051
dut1 = 0.1550675
elong = -0.527800806
phi = -1.2345856
hm = 2738.0
xp = 2.47230737e-7
yp = 1.82640464e-6
phpa = 731.0
tc = 12.8
rh = 0.59
wl = 0.55
aob, zob, hob, dob, rob = eraAtio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
@test_approx_eq_eps aob  0.09233952224794989993 1e-12
@test_approx_eq_eps zob  1.407758704513722461 1e-12
@test_approx_eq_eps hob  -0.09247619879782006106 1e-12
@test_approx_eq_eps dob  0.1717653435758265198 1e-12
@test_approx_eq_eps rob  2.710085107986886201 1e-12

# eraAtioq
utc1 = 2456384.5
utc2 = 0.969254051
dut1 = 0.1550675
elong = -0.527800806
phi = -1.2345856
hm = 2738.0
xp = 2.47230737e-7
yp = 1.82640464e-6
phpa = 731.0
tc = 12.8
rh = 0.59
wl = 0.55
astrom = eraApio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                   phpa, tc, rh, wl)
ri = 2.710121572969038991
di = 0.1729371367218230438
aob, zob, hob, dob, rob = eraAtioq(ri, di, astrom)
@test_approx_eq_eps aob  0.09233952224794989993 1e-12
@test_approx_eq_eps zob  1.407758704513722461 1e-12
@test_approx_eq_eps hob  -0.09247619879782006106 1e-12
@test_approx_eq_eps dob  0.1717653435758265198 1e-12
@test_approx_eq_eps rob  2.710085107986886201 1e-12

# eraAtoc13
utc1 = 2456384.5
utc2 = 0.969254051
dut1 = 0.1550675
elong = -0.527800806
phi = -1.2345856
hm = 2738.0
xp = 2.47230737e-7
yp = 1.82640464e-6
phpa = 731.0
tc = 12.8
rh = 0.59
wl = 0.55
ob1 = 2.710085107986886201
ob2 = 0.1717653435758265198
rc, dc = eraAtoc13('r', ob1, ob2, utc1, utc2, dut1,
                   elong, phi, hm, xp, yp, phpa, tc, rh, wl)
@test_approx_eq_eps rc  2.709956744661000609 1e-12
@test_approx_eq_eps dc  0.1741696500895398562 1e-12
ob1 = -0.09247619879782006106
ob2 = 0.1717653435758265198
rc, dc = eraAtoc13('h', ob1, ob2, utc1, utc2, dut1,
                   elong, phi, hm, xp, yp, phpa, tc, rh, wl)
@test_approx_eq_eps rc  2.709956744661000609 1e-12
@test_approx_eq_eps dc  0.1741696500895398562 1e-12
ob1 = 0.09233952224794989993
ob2 = 1.407758704513722461
rc, dc = eraAtoc13('a', ob1, ob2, utc1, utc2, dut1,
                   elong, phi, hm, xp, yp, phpa, tc, rh, wl)
@test_approx_eq_eps rc  2.709956744661000609 1e-12
@test_approx_eq_eps dc  0.1741696500895398562 1e-12

# eraAtoi13
utc1 = 2456384.5
utc2 = 0.969254051
dut1 = 0.1550675
elong = -0.527800806
phi = -1.2345856
hm = 2738.0
xp = 2.47230737e-7
yp = 1.82640464e-6
phpa = 731.0
tc = 12.8
rh = 0.59
wl = 0.55
ob1 = 2.710085107986886201
ob2 = 0.1717653435758265198
ri, di = eraAtoi13('r', ob1, ob2, utc1, utc2, dut1,
                   elong, phi, hm, xp, yp, phpa, tc, rh, wl) 
@test_approx_eq_eps ri  2.710121574449135955 1e-12
@test_approx_eq_eps di  0.1729371839114567725 1e-12
ob1 = -0.09247619879782006106
ob2 = 0.1717653435758265198
ri, di = eraAtoi13('h', ob1, ob2, utc1, utc2, dut1,
                   elong, phi, hm, xp, yp, phpa, tc, rh, wl) 
@test_approx_eq_eps ri  2.710121574449135955 1e-12
@test_approx_eq_eps di  0.1729371839114567725 1e-12
ob1 = 0.09233952224794989993
ob2 = 1.407758704513722461
ri, di = eraAtoi13('a', ob1, ob2, utc1, utc2, dut1,
                   elong, phi, hm, xp, yp, phpa, tc, rh, wl) 
@test_approx_eq_eps ri  2.710121574449135955 1e-12
@test_approx_eq_eps di  0.1729371839114567728 1e-12

# eraAtoiq
utc1 = 2456384.5
utc2 = 0.969254051
dut1 = 0.1550675
elong = -0.527800806
phi = -1.2345856
hm = 2738.0
xp = 2.47230737e-7
yp = 1.82640464e-6
phpa = 731.0
tc = 12.8
rh = 0.59
wl = 0.55
astrom = eraApio13(utc1, utc2, dut1,
                   elong, phi, hm, xp, yp, phpa, tc, rh, wl)
ob1 = 2.710085107986886201
ob2 = 0.1717653435758265198
ri, di = eraAtoiq('r', ob1, ob2, astrom)
@test_approx_eq_eps ri  2.710121574449135955 1e-12
@test_approx_eq_eps di  0.1729371839114567725 1e-12
ob1 = -0.09247619879782006106
ob2 = 0.1717653435758265198
ri, di = eraAtoiq('h', ob1, ob2, astrom)
@test_approx_eq_eps ri  2.710121574449135955 1e-12
@test_approx_eq_eps di  0.1729371839114567725 1e-12
ob1 = 0.09233952224794989993
ob2 = 1.407758704513722461
ri, di = eraAtoiq('a', ob1, ob2, astrom)
@test_approx_eq_eps ri  2.710121574449135955 1e-12
@test_approx_eq_eps di  0.1729371839114567728 1e-12

# eraBi00
dp, de, dr = eraBi00()
@test_approx_eq_eps dp  -0.2025309152835086613e-6 1e-15
@test_approx_eq_eps de  -0.3306041454222147847e-7 1e-15
@test_approx_eq_eps dr  -0.7078279744199225506e-7 1e-15

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

# eraBpn2xy
rbpn = [[9.999962358680738e-1,-2.516417057665452e-3,-1.093569785342370e-3];
        [2.516462370370876e-3,9.999968329010883e-1,4.006159587358310e-5];
        [1.093465510215479e-3,-4.281337229063151e-5,9.999994012499173e-1]]
x, y = eraBpn2xy(rbpn)
@test_approx_eq_eps x  1.093465510215479e-3 1e-12
@test_approx_eq_eps y  -4.281337229063151e-5 1e-12

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

# eraC2ibpn
rbpn = [[9.999962358680738e-1,-2.516417057665452e-3,-1.093569785342370e-3];
        [2.516462370370876e-3,9.999968329010883e-1,4.006159587358310e-5];
        [1.093465510215479e-3,-4.281337229063151e-5,9.999994012499173e-1]]
rc2i = eraC2ibpn(2400000.5, 50123.9999, rbpn)
@test_approx_eq_eps rc2i[1]  0.9999994021664089977 1e-12
@test_approx_eq_eps rc2i[2]  -0.3869195948017503664e-8 1e-12
@test_approx_eq_eps rc2i[3]  -0.1093465511383285076e-2 1e-12
@test_approx_eq_eps rc2i[4]  0.5068413965715446111e-7 1e-12
@test_approx_eq_eps rc2i[5]  0.9999999990835075686 1e-12
@test_approx_eq_eps rc2i[6]  0.4281334246452708915e-4 1e-12
@test_approx_eq_eps rc2i[7]  0.1093465510215479000e-2 1e-12
@test_approx_eq_eps rc2i[8]  -0.4281337229063151000e-4 1e-12
@test_approx_eq_eps rc2i[9]  0.9999994012499173103 1e-12

# eraC2s
t,p = eraC2s([100.,-50.,25.])
@test_approx_eq_eps t  -0.4636476090008061162 1e-15
@test_approx_eq_eps p  0.2199879773954594463 1e-15

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

# eraC2t00a
tta = 2400000.5
uta = 2400000.5
ttb = 53736.0
utb = 53736.0
xp = 2.55060238e-7
yp = 1.860359247e-6
rc2t = eraC2t00a(tta, ttb, uta, utb, xp, yp)
@test_approx_eq_eps rc2t[1]  -0.1810332128307182668 1e-12
@test_approx_eq_eps rc2t[2]  0.9834769806938457836 1e-12
@test_approx_eq_eps rc2t[3]  0.6555535638688341725e-4 1e-12
@test_approx_eq_eps rc2t[4]  -0.9834768134135984552 1e-12
@test_approx_eq_eps rc2t[5]  -0.1810332203649520727 1e-12
@test_approx_eq_eps rc2t[6]  0.5749801116141056317e-3 1e-12
@test_approx_eq_eps rc2t[7]  0.5773474014081406921e-3 1e-12
@test_approx_eq_eps rc2t[8]  0.3961832391770163647e-4 1e-12
@test_approx_eq_eps rc2t[9]  0.9999998325501692289 1e-12

# eraC2t00b
rc2t = eraC2t00b(tta, ttb, uta, utb, xp, yp)
@test_approx_eq_eps rc2t[1]  -0.1810332128439678965 1e-12
@test_approx_eq_eps rc2t[2]  0.9834769806913872359 1e-12
@test_approx_eq_eps rc2t[3]  0.6555565082458415611e-4 1e-12
@test_approx_eq_eps rc2t[4]  -0.9834768134115435923 1e-12
@test_approx_eq_eps rc2t[5]  -0.1810332203784001946 1e-12
@test_approx_eq_eps rc2t[6]  0.5749793922030017230e-3 1e-12
@test_approx_eq_eps rc2t[7]  0.5773467471863534901e-3 1e-12
@test_approx_eq_eps rc2t[8]  0.3961790411549945020e-4 1e-12
@test_approx_eq_eps rc2t[9]  0.9999998325505635738 1e-12

# eraC2t06a
rc2t = eraC2t06a(tta, ttb, uta, utb, xp, yp)
@test_approx_eq_eps rc2t[1]  -0.1810332128305897282 1e-12
@test_approx_eq_eps rc2t[2]  0.9834769806938592296 1e-12
@test_approx_eq_eps rc2t[3]  0.6555550962998436505e-4 1e-12
@test_approx_eq_eps rc2t[4]  -0.9834768134136214897 1e-12
@test_approx_eq_eps rc2t[5]  -0.1810332203649130832 1e-12
@test_approx_eq_eps rc2t[6]  0.5749800844905594110e-3 1e-12
@test_approx_eq_eps rc2t[7]  0.5773474024748545878e-3 1e-12
@test_approx_eq_eps rc2t[8]  0.3961816829632690581e-4 1e-12
@test_approx_eq_eps rc2t[9]  0.9999998325501747785 1e-12

# eraC2tcio
c = [[0.9999998323037164738,0.5581526271714303683e-9,-0.5791308477073443903e-3];
     [-0.2384266227524722273e-7,0.9999999991917404296,-0.4020594955030704125e-4];
     [0.5791308472168153320e-3,.4020595661593994396e-4,0.9999998314954572365]]
era = 1.75283325530307
p = [[0.9999999999999674705,-0.1367174580728847031e-10,0.2550602379999972723e-6];
     [0.1414624947957029721e-10,0.9999999999982694954,-0.1860359246998866338e-5];
     [-0.2550602379741215275e-6,0.1860359247002413923e-5,0.9999999999982369658]]
rc2t = eraC2tcio(c,era,p)
@test_approx_eq_eps rc2t[1]  -0.1810332128307110439 1e-12
@test_approx_eq_eps rc2t[2]  0.9834769806938470149 1e-12
@test_approx_eq_eps rc2t[3]  0.6555535638685466874e-4 1e-12
@test_approx_eq_eps rc2t[4]  -0.9834768134135996657 1e-12
@test_approx_eq_eps rc2t[5]  -0.1810332203649448367 1e-12
@test_approx_eq_eps rc2t[6]  0.5749801116141106528e-3 1e-12
@test_approx_eq_eps rc2t[7]  0.5773474014081407076e-3 1e-12
@test_approx_eq_eps rc2t[8]  0.3961832391772658944e-4 1e-12
@test_approx_eq_eps rc2t[9]  0.9999998325501691969 1e-12

# eraC2teqx
c = [[0.9999989440476103608,-0.1332881761240011518e-2,-0.5790767434730085097e-3];
     [0.1332858254308954453e-2,0.9999991109044505944,-0.4097782710401555759e-4];
     [0.5791308472168153320e-3,0.4020595661593994396e-4,0.9999998314954572365]]
gst = 1.754166138040730516
p = [[0.9999999999999674705,-0.1367174580728847031e-10,0.2550602379999972723e-6];
     [0.1414624947957029721e-10,0.9999999999982694954,-0.1860359246998866338e-5];
     [-0.2550602379741215275e-6,0.1860359247002413923e-5,0.9999999999982369658]]
rc2t = eraC2teqx(c,gst,p)
@test_approx_eq_eps rc2t[1]  -0.1810332128528685730 1e-12
@test_approx_eq_eps rc2t[2]  0.9834769806897685071 1e-12
@test_approx_eq_eps rc2t[3]  0.6555535639982634449e-4 1e-12
@test_approx_eq_eps rc2t[4]  -0.9834768134095211257 1e-12
@test_approx_eq_eps rc2t[5]  -0.1810332203871023800 1e-12
@test_approx_eq_eps rc2t[6]  0.5749801116126438962e-3 1e-12
@test_approx_eq_eps rc2t[7]  0.5773474014081539467e-3 1e-12
@test_approx_eq_eps rc2t[8]  0.3961832391768640871e-4 1e-12
@test_approx_eq_eps rc2t[9]  0.9999998325501691969 1e-12

# eraC2tpe
tta = 2400000.5
uta = 2400000.5
ttb = 53736.0
utb = 53736.0
deps =  0.4090789763356509900
dpsi = -0.9630909107115582393e-5
xp = 2.55060238e-7
yp = 1.860359247e-6
rc2t = eraC2tpe(tta, ttb, uta, utb, dpsi, deps, xp, yp)
@test_approx_eq_eps rc2t[1]  -0.1813677995763029394 1e-12
@test_approx_eq_eps rc2t[2]  0.9023482206891683275 1e-12
@test_approx_eq_eps rc2t[3]  -0.3909902938641085751 1e-12
@test_approx_eq_eps rc2t[4]  -0.9834147641476804807 1e-12
@test_approx_eq_eps rc2t[5]  -0.1659883635434995121 1e-12
@test_approx_eq_eps rc2t[6]  0.7309763898042819705e-1 1e-12
@test_approx_eq_eps rc2t[7]  0.1059685430673215247e-2 1e-12
@test_approx_eq_eps rc2t[8]  0.3977631855605078674 1e-12
@test_approx_eq_eps rc2t[9]  0.9174875068792735362 1e-12

# eraC2txy
tta = 2400000.5
uta = 2400000.5
ttb = 53736.0
utb = 53736.0
x = 0.5791308486706011000e-3
y = 0.4020579816732961219e-4
xp = 2.55060238e-7
yp = 1.860359247e-6
rc2t = eraC2txy(tta, ttb, uta, utb, x, y, xp, yp)
@test_approx_eq_eps rc2t[1]  -0.1810332128306279253 1e-12
@test_approx_eq_eps rc2t[2]  0.9834769806938520084 1e-12
@test_approx_eq_eps rc2t[3]  0.6555551248057665829e-4 1e-12
@test_approx_eq_eps rc2t[4]  -0.9834768134136142314 1e-12
@test_approx_eq_eps rc2t[5]  -0.1810332203649529312 1e-12
@test_approx_eq_eps rc2t[6]  0.5749800843594139912e-3 1e-12
@test_approx_eq_eps rc2t[7]  0.5773474028619264494e-3 1e-12
@test_approx_eq_eps rc2t[8]  0.3961816546911624260e-4 1e-12
@test_approx_eq_eps rc2t[9]  0.9999998325501746670 1e-12

# eraCal2jd
dmj0, dmj = eraCal2jd(2003, 6, 1)
@test_approx_eq_eps dmj0  2400000.5 1e-9
@test_approx_eq_eps dmj  52791.0 1e-9

# eraCp
p = [0.3,1.2,-2.5]
c = eraCp(p)
@test_approx_eq_eps c[1]  0.3 1e-12
@test_approx_eq_eps c[2]  1.2 1e-12
@test_approx_eq_eps c[3]  -2.5 1e-12

# eraCpv
pv = [[0.3,1.2,-2.5];
      [-0.5,3.1,0.9]]
c = eraCpv(pv)
@test_approx_eq_eps c[1]  0.3 1e-12
@test_approx_eq_eps c[2]  1.2 1e-12
@test_approx_eq_eps c[3]  -2.5 1e-12
@test_approx_eq_eps c[4]  -0.5 1e-12
@test_approx_eq_eps c[5]  3.1 1e-12
@test_approx_eq_eps c[6]  0.9 1e-12

# eraCr
a = [[2.0,3.0,2.0];
     [3.0,2.0,3.0];
     [3.0,4.0,5.0]]
c = eraCr(a)
@test_approx_eq_eps c[1]  2.0 1e-12
@test_approx_eq_eps c[2]  3.0 1e-12
@test_approx_eq_eps c[3]  2.0 1e-12
@test_approx_eq_eps c[4]  3.0 1e-12
@test_approx_eq_eps c[5]  2.0 1e-12
@test_approx_eq_eps c[6]  3.0 1e-12
@test_approx_eq_eps c[7]  3.0 1e-12
@test_approx_eq_eps c[8]  4.0 1e-12
@test_approx_eq_eps c[9]  5.0 1e-12

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

# eraDtdb
d = eraDtdb(2448939.5, 0.123, 0.76543, 5.0123, 5525.242, 3190.0)
@test_approx_eq_eps d  -0.1280368005936998991e-2 1e-17

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

# eraEect00
ct = eraEect00(2400000.5, 53736.0)
@test_approx_eq_eps ct  0.2046085004885125264e-8 1e-20

# eraEo06a
eo = eraEo06a(2400000.5, 53736.0)
@test_approx_eq_eps eo  -0.1332882371941833644e-2 1e-15

# eraEors
r = [[0.9999989440476103608,-0.1332881761240011518e-2,-0.5790767434730085097e-3];
     [0.1332858254308954453e-2,0.9999991109044505944,-0.4097782710401555759e-4];
     [0.5791308472168153320e-3,0.4020595661593994396e-4,0.9999998314954572365]]
s = -0.1220040848472271978e-7
eo = eraEors(r,s)
@test_approx_eq_eps eo  -0.1332882715130744606e-2 1e-15

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

# eraFad03
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

# eraFk52h
r5  =  1.76779433
d5  = -0.2917517103
dr5 = -1.91851572e-7
dd5 = -5.8468475e-6
px5 =  0.379210
rv5 = -7.6
rh, dh, drh, ddh, pxh, rvh = eraFk52h(r5, d5, dr5, dd5, px5, rv5)
@test_approx_eq_eps rh  1.767794226299947632 1e-14
@test_approx_eq_eps dh   -0.2917516070530391757 1e-14
@test_approx_eq_eps drh  -0.19618741256057224e-6 1e-19
@test_approx_eq_eps ddh  -0.58459905176693911e-5 1e-19
@test_approx_eq_eps pxh   0.37921 1e-14
@test_approx_eq_eps rvh  -7.6000000940000254 1e-11

# eraFk5hz
r5 =  1.76779433
d5 = -0.2917517103
rh,dh = eraFk5hz(r5,d5,2400000.5, 54479.0)
@test_approx_eq_eps rh  1.767794191464423978 1e-12
@test_approx_eq_eps dh  -0.2917516001679884419 1e-12

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

# eraFw2xy
gamb = -0.2243387670997992368e-5
phib =  0.4091014602391312982
psi  = -0.9501954178013015092e-3
eps  =  0.4091014316587367472
x, y = eraFw2xy(gamb, phib, psi, eps)
@test_approx_eq_eps x -0.3779734957034082790e-3 1e-14
@test_approx_eq_eps y -0.1924880848087615651e-6 1e-14

# eraGc2gd
xyz = [2.e6, 3.e6, 5.244e6]
e,p,h = eraGc2gd(1, xyz)
@test_approx_eq_eps e  0.98279372324732907 1e-14
@test_approx_eq_eps p  0.97160184819075459 1e-14
@test_approx_eq_eps h  331.41724614260599 1e-8
e,p,h = eraGc2gd(2, xyz)
@test_approx_eq_eps e  0.98279372324732907 1e-14
@test_approx_eq_eps p  0.97160184820607853 1e-14
@test_approx_eq_eps h  331.41731754844348 1e-8
e,p,h = eraGc2gd(3, xyz)
@test_approx_eq_eps e  0.98279372324732907 1e-14
@test_approx_eq_eps p  0.97160181811015119 1e-14
@test_approx_eq_eps h  333.27707261303181 1e-8

# eraGc2gde
a = 6378136.0
f = 0.0033528
xyz = [2e6, 3e6, 5.244e6]
e, p, h = eraGc2gde(a, f, xyz)
@test_approx_eq_eps e  0.98279372324732907 1e-14
@test_approx_eq_eps p  0.97160183775704115 1e-14
@test_approx_eq_eps h  332.36862495764397 1e-8

# eraGd2gc
e = 3.1
p = -0.5
h = 2500.0
xyz = eraGd2gc(1,e,p,h)
@test_approx_eq_eps xyz[1]  -5599000.5577049947 1e-7
@test_approx_eq_eps xyz[2]  233011.67223479203 1e-7
@test_approx_eq_eps xyz[3]  -3040909.4706983363 1e-7
xyz = eraGd2gc(2,e,p,h)
@test_approx_eq_eps xyz[1]  -5599000.5577260984 1e-7
@test_approx_eq_eps xyz[2]  233011.6722356703 1e-7
@test_approx_eq_eps xyz[3]  -3040909.4706095476 1e-7
xyz = eraGd2gc(3,e,p,h)
@test_approx_eq_eps xyz[1]  -5598998.7626301490 1e-7
@test_approx_eq_eps xyz[2]  233011.5975297822 1e-7
@test_approx_eq_eps xyz[3]  -3040908.6861467111 1e-7

# eraGd2gce
a = 6378136.0
f = 0.0033528
e = 3.1
p = -0.5
h = 2500.0
xyz = eraGd2gce(a, f, e, p, h)
@test_approx_eq_eps xyz[1]  -5598999.6665116328 1e-7
@test_approx_eq_eps xyz[2]  233011.63514630572 1e-7
@test_approx_eq_eps xyz[3]  -3040909.0517314132 1e-7

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

# eraGst06
rnpb = [[0.9999989440476103608,-0.1332881761240011518e-2,-0.5790767434730085097e-3];
        [0.1332858254308954453e-2,0.9999991109044505944,-0.4097782710401555759e-4];
        [0.5791308472168153320e-3,0.4020595661593994396e-4,0.9999998314954572365]]
g = eraGst06(2400000.5, 53736.0, 2400000.5, 53736.0, rnpb)
@test_approx_eq_eps g  1.754166138018167568 1e-14

# eraGst06a
g = eraGst06a(2400000.5, 53736.0, 2400000.5, 53736.0)
@test_approx_eq_eps g  1.754166137675019159 1e-14

# eraGst94
g = eraGst94(2400000.5, 53736.0)
@test_approx_eq_eps g  1.754166136020645203 1e-14

# eraH2fk5
rh  =  1.767794352
dh  = -0.2917512594
drh = -2.76413026e-6
ddh = -5.92994449e-6
pxh =  0.379210
rvh = -7.6
r5, d5, dr5, dd5, px5, rv5 = eraH2fk5(rh, dh, drh, ddh, pxh, rvh)
@test_approx_eq_eps r5  1.767794455700065506 1e-13
@test_approx_eq_eps d5  -0.2917513626469638890 1e-13
@test_approx_eq_eps dr5  -0.27597945024511204e-5 1e-18
@test_approx_eq_eps dd5  -0.59308014093262838e-5 1e-18
@test_approx_eq_eps px5  0.37921 1e-13
@test_approx_eq_eps rv5  -7.6000001309071126 1e-10

# eraFk5hip
r5h, s5h = eraFk5hip()
@test_approx_eq_eps r5h[1]  0.9999999999999928638 1e-14
@test_approx_eq_eps r5h[2]  0.1110223351022919694e-6 1e-17
@test_approx_eq_eps r5h[3]  0.4411803962536558154e-7 1e-17
@test_approx_eq_eps r5h[4]  -0.1110223308458746430e-6 1e-17
@test_approx_eq_eps r5h[5]  0.9999999999999891830 1e-14
@test_approx_eq_eps r5h[6]  -0.9647792498984142358e-7 1e-17
@test_approx_eq_eps r5h[7]  -0.4411805033656962252e-7 1e-17
@test_approx_eq_eps r5h[8]  0.9647792009175314354e-7 1e-17
@test_approx_eq_eps r5h[9]  0.9999999999999943728 1e-14
@test_approx_eq_eps s5h[1]  -0.1454441043328607981e-8 1e-17
@test_approx_eq_eps s5h[2]  0.2908882086657215962e-8 1e-17
@test_approx_eq_eps s5h[3]  0.3393695767766751955e-8 1e-17

# eraHfk5z
rh =  1.767794352
dh = -0.2917512594
r5, d5, dr5, dd5 = eraHfk5z(rh, dh, 2400000.5, 54479.0)
@test_approx_eq_eps r5  1.767794490535581026 1e-13
@test_approx_eq_eps d5  -0.2917513695320114258 1e-14
@test_approx_eq_eps dr5  0.4335890983539243029e-8 1e-22
@test_approx_eq_eps dd5  -0.8569648841237745902e-9 1e-23

# eraJd2cal
y, m, d, fd = eraJd2cal(2400000.5, 50123.9999)
@test (y, m, d) == (1996, 2, 10)
@test_approx_eq_eps fd  0.9999 1e-7

# eraJdcalf
y, m, d, fd = eraJdcalf(4, 2400000.5, 50123.9999)
@test (y, m, d, fd) == (1996, 2, 10, 9999)

# eraLd
bm = 0.00028574
p = [-0.763276255, -0.608633767, -0.216735543]
q = [-0.763276255, -0.608633767, -0.216735543]
e = [0.76700421, 0.605629598, 0.211937094]
em = 8.91276983
dlim = 3e-10
p1 = eraLd(bm, p, q, e, em, dlim)
@test_approx_eq_eps p1[1]  -0.7632762548968159627 1e-12
@test_approx_eq_eps p1[2]  -0.6086337670823762701 1e-12
@test_approx_eq_eps p1[3]  -0.2167355431320546947 1e-12

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

# eraLdsun
p = [-0.763276255, -0.608633767, -0.216735543]
e = [-0.973644023, -0.20925523, -0.0907169552]
em = 0.999809214
p1 = eraLdsun(p, e, em)
@test_approx_eq_eps p1[1]  -0.7632762580731413169 1e-12
@test_approx_eq_eps p1[2]  -0.6086337635262647900 1e-12
@test_approx_eq_eps p1[3]  -0.2167355419322321302 1e-12

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

# eraP06e
eps0,psia,oma,bpa,bqa,pia,bpia,epsa,chia,za,zetaa,thetaa,pa,gam,phi,psi = eraP06e(2400000.5, 52541.0)
@test_approx_eq_eps eps0  0.4090926006005828715 1e-14
@test_approx_eq_eps psia  0.6664369630191613431e-3 1e-14
@test_approx_eq_eps oma   0.4090925973783255982 1e-14
@test_approx_eq_eps bpa  0.5561149371265209445e-6 1e-14
@test_approx_eq_eps bqa  -0.6191517193290621270e-5 1e-14
@test_approx_eq_eps pia  0.6216441751884382923e-5 1e-14
@test_approx_eq_eps bpia  3.052014180023779882 1e-14
@test_approx_eq_eps epsa  0.4090864054922431688 1e-14
@test_approx_eq_eps chia  0.1387703379530915364e-5 1e-14
@test_approx_eq_eps za  0.2921789846651790546e-3 1e-14
@test_approx_eq_eps zetaa  0.3178773290332009310e-3 1e-14
@test_approx_eq_eps thetaa  0.2650932701657497181e-3 1e-14
@test_approx_eq_eps pa  0.6651637681381016344e-3 1e-14
@test_approx_eq_eps gam  0.1398077115963754987e-5 1e-14
@test_approx_eq_eps phi  0.4090864090837462602 1e-14
@test_approx_eq_eps psi  0.6664464807480920325e-3 1e-14

# eraP2pv
pv = eraP2pv([0.25,1.2,3.0])
@test_approx_eq_eps pv[1]  0.25 1e-12
@test_approx_eq_eps pv[2]  1.2 1e-12
@test_approx_eq_eps pv[3]  3.0 1e-12
@test_approx_eq_eps pv[4]  0.0 1e-12
@test_approx_eq_eps pv[5]  0.0 1e-12
@test_approx_eq_eps pv[6]  0.0 1e-12

# eraP2s
theta, phi, r = eraP2s([100.,-50.,25.])
@test_approx_eq_eps theta  -0.4636476090008061162 1e-12
@test_approx_eq_eps phi  0.2199879773954594463 1e-12
@test_approx_eq_eps r  114.5643923738960002 1e-12

# eraPap
a = [1.,0.1,0.2]
b = [-3.,1e-3,0.2]
theta = eraPap(a,b)
@test_approx_eq_eps theta  0.3671514267841113674 1e-12

# eraPas
p = eraPas(1.0,0.1,0.2,-1.0)
@test_approx_eq_eps p  -2.724544922932270424 1e-12

# eraPb06
bzeta, bz, btheta = eraPb06(2400000.5, 50123.9999)
@test_approx_eq_eps bzeta  -0.5092634016326478238e-3 1e-12
@test_approx_eq_eps bz  -0.3602772060566044413e-3 1e-12
@test_approx_eq_eps btheta  -0.3779735537167811177e-3 1e-12

# eraPfw06
gamb, phib, psib, epsa = eraPfw06(2400000.5, 50123.9999)
@test_approx_eq_eps gamb  -0.2243387670997995690e-5 1e-16
@test_approx_eq_eps phib   0.4091014602391312808 1e-12
@test_approx_eq_eps psib  -0.9501954178013031895e-3 1e-14
@test_approx_eq_eps epsa   0.4091014316587367491 1e-12

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

# eraPm
m = eraPm([0.3,1.2,-2.5])
@test_approx_eq_eps m 2.789265136196270604 1e-14

# eraPmp
a = [2.0,2.0,3.0]
b = [1.0,3.0,4.0]
amb= eraPmp(a,b)
@test_approx_eq_eps amb[1]  1.0 1e-12
@test_approx_eq_eps amb[2]  -1.0 1e-12
@test_approx_eq_eps amb[3]  -1.0 1e-12

# eraPmpx
rc = 1.234
dc = 0.789
pr = 1e-5
pd = -2e-5
px = 1e-2
rv = 10.0
pmt = 8.75
pob = [0.9, 0.4, 0.1]
pco = eraPmpx(rc, dc, pr, pd, px, rv, pmt, pob)
@test_approx_eq_eps pco[1]  0.2328137623960308440 1e-12
@test_approx_eq_eps pco[2]  0.6651097085397855317 1e-12
@test_approx_eq_eps pco[3]  0.7095257765896359847 1e-12

# eraPmsafe
ra1 = 1.234
dec1 = 0.789
pmr1 = 1e-5
pmd1 = -2e-5
px1 = 1e-2
rv1 = 10.0
ep1a = 2400000.5
ep1b = 48348.5625
ep2a = 2400000.5
ep2b = 51544.5
ra2, dec2, pmr2, pmd2, px2, rv2 = eraPmsafe(ra1, dec1, pmr1, pmd1, px1,rv1, ep1a, ep1b, ep2a, ep2b)
@test_approx_eq_eps ra2  1.234087484501017061 1e-12
@test_approx_eq_eps dec2  0.7888249982450468574 1e-12
@test_approx_eq_eps pmr2  0.9996457663586073988e-5 1e-12
@test_approx_eq_eps pmd2  -0.2000040085106737816e-4 1e-16
@test_approx_eq_eps px2  0.9999997295356765185e-2 1e-12
@test_approx_eq_eps rv2  10.38468380113917014 1e-10

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

# eraPn
r, u = eraPn([0.3,1.2,-2.5])
@test_approx_eq_eps r  2.789265136196270604 1e-12
@test_approx_eq_eps u[1]  0.1075552109073112058 1e-12
@test_approx_eq_eps u[2]  0.4302208436292448232 1e-12
@test_approx_eq_eps u[3]  -0.8962934242275933816 1e-12

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

# eraPn00a
dpsi, deps, epsa, rb, rp, rbp, rn, rbpn = eraPn00a(2400000.5, 53736.0)
@test_approx_eq_eps dpsi  -0.9630909107115518431e-5 1e-12
@test_approx_eq_eps deps   0.4063239174001678710e-4 1e-12
@test_approx_eq_eps epsa   0.4090791789404229916 1e-12
@test_approx_eq_eps rb[1]  0.9999999999999942498 1e-12
@test_approx_eq_eps rb[2]  -0.7078279744199196626e-7 1e-16
@test_approx_eq_eps rb[3]  0.8056217146976134152e-7 1e-16
@test_approx_eq_eps rb[4]  0.7078279477857337206e-7 1e-16
@test_approx_eq_eps rb[5]  0.9999999999999969484 1e-12
@test_approx_eq_eps rb[6]  0.3306041454222136517e-7 1e-16
@test_approx_eq_eps rb[7]  -0.8056217380986972157e-7 1e-16
@test_approx_eq_eps rb[8]  -0.3306040883980552500e-7 1e-16
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
@test_approx_eq_eps rn[1]  0.9999999999536227949 1e-12
@test_approx_eq_eps rn[2]  0.8836238544090873336e-5 1e-14
@test_approx_eq_eps rn[3]  0.3830835237722400669e-5 1e-14
@test_approx_eq_eps rn[4]  -0.8836082880798569274e-5 1e-14
@test_approx_eq_eps rn[5]  0.9999999991354655028 1e-12
@test_approx_eq_eps rn[6]  -0.4063240865362499850e-4 1e-14
@test_approx_eq_eps rn[7]  -0.3831194272065995866e-5 1e-14
@test_approx_eq_eps rn[8]  0.4063237480216291775e-4 1e-14
@test_approx_eq_eps rn[9]  0.9999999991671660338 1e-12
@test_approx_eq_eps rbpn[1]  0.9999989440476103435 1e-12
@test_approx_eq_eps rbpn[2]  -0.1332881761240011763e-2 1e-14
@test_approx_eq_eps rbpn[3]  -0.5790767434730085751e-3 1e-14
@test_approx_eq_eps rbpn[4]  0.1332858254308954658e-2 1e-14
@test_approx_eq_eps rbpn[5]  0.9999991109044505577 1e-12
@test_approx_eq_eps rbpn[6]  -0.4097782710396580452e-4 1e-14
@test_approx_eq_eps rbpn[7]  0.5791308472168152904e-3 1e-14
@test_approx_eq_eps rbpn[8]  0.4020595661591500259e-4 1e-14
@test_approx_eq_eps rbpn[9]  0.9999998314954572304 1e-12

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

# eraPn06a
dpsi, deps, epsa, rb, rp, rbp, rn, rbpn = eraPn06a(2400000.5, 53736.0)
@test_approx_eq_eps dpsi  -0.9630912025820308797e-5 1e-12
@test_approx_eq_eps deps   0.4063238496887249798e-4 1e-12
@test_approx_eq_eps epsa   0.4090789763356509926 1e-12
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
@test_approx_eq_eps rn[1]  0.9999999999536227668 1e-12
@test_approx_eq_eps rn[2]  0.8836241998111535233e-5 1e-14
@test_approx_eq_eps rn[3]  0.3830834608415287707e-5 1e-14
@test_approx_eq_eps rn[4]  -0.8836086334870740138e-5 1e-14
@test_approx_eq_eps rn[5]  0.9999999991354657474 1e-12
@test_approx_eq_eps rn[6]  -0.4063240188248455065e-4 1e-14
@test_approx_eq_eps rn[7]  -0.3831193642839398128e-5 1e-14
@test_approx_eq_eps rn[8]  0.4063236803101479770e-4 1e-14
@test_approx_eq_eps rn[9]  0.9999999991671663114 1e-12
@test_approx_eq_eps rbpn[1]  0.9999989440480669738 1e-12
@test_approx_eq_eps rbpn[2]  -0.1332881418091915973e-2 1e-14
@test_approx_eq_eps rbpn[3]  -0.5790767447612042565e-3 1e-14
@test_approx_eq_eps rbpn[4]  0.1332857911250989133e-2 1e-14
@test_approx_eq_eps rbpn[5]  0.9999991109049141908 1e-12
@test_approx_eq_eps rbpn[6]  -0.4097767128546784878e-4 1e-14
@test_approx_eq_eps rbpn[7]  0.5791308482835292617e-3 1e-14
@test_approx_eq_eps rbpn[8]  0.4020580099454020310e-4 1e-14
@test_approx_eq_eps rbpn[9]  0.9999998314954628695 1e-12

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

# eraPn00b
dpsi, deps, epsa, rb, rp, rbp, rn, rbpn = eraPn00b(2400000.5, 53736.0)
@test_approx_eq_eps dpsi  -0.9632552291148362783e-5 1e-12
@test_approx_eq_eps deps   0.4063197106621159367e-4 1e-12
@test_approx_eq_eps epsa   0.4090791789404229916 1e-12
@test_approx_eq_eps rb[1]  0.9999999999999942498 1e-12
@test_approx_eq_eps rb[2]  -0.7078279744199196626e-7 1e-16
@test_approx_eq_eps rb[3]  0.8056217146976134152e-7 1e-16
@test_approx_eq_eps rb[4]  0.7078279477857337206e-7 1e-16
@test_approx_eq_eps rb[5]  0.9999999999999969484 1e-12
@test_approx_eq_eps rb[6]  0.3306041454222136517e-7 1e-16
@test_approx_eq_eps rb[7]  -0.8056217380986972157e-7 1e-16
@test_approx_eq_eps rb[8]  -0.3306040883980552500e-7 1e-16
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
@test_approx_eq_eps rn[2]  0.8837746144871248011e-5 1e-14
@test_approx_eq_eps rn[3]  0.3831488838252202945e-5 1e-14
@test_approx_eq_eps rn[4]  -0.8837590456632304720e-5 1e-14
@test_approx_eq_eps rn[5]  0.9999999991354692733 1e-12
@test_approx_eq_eps rn[6]  -0.4063198798559591654e-4 1e-14
@test_approx_eq_eps rn[7]  -0.3831847930134941271e-5 1e-14
@test_approx_eq_eps rn[8]  0.4063195412258168380e-4 1e-14
@test_approx_eq_eps rn[9]  0.9999999991671806225 1e-12
@test_approx_eq_eps rbpn[1]  0.9999989440499982806 1e-12
@test_approx_eq_eps rbpn[2]  -0.1332880253640849194e-2 1e-14
@test_approx_eq_eps rbpn[3]  -0.5790760898731091166e-3 1e-14
@test_approx_eq_eps rbpn[4]  0.1332856746979949638e-2 1e-14
@test_approx_eq_eps rbpn[5]  0.9999991109064768883 1e-12
@test_approx_eq_eps rbpn[6]  -0.4097740555723081811e-4 1e-14
@test_approx_eq_eps rbpn[7]  0.5791301929950208873e-3 1e-14
@test_approx_eq_eps rbpn[8]  0.4020553681373720832e-4 1e-14
@test_approx_eq_eps rbpn[9]  0.9999998314958529887 1e-12

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

# eraPpp
apb = eraPpp([2.0,2.0,3.0],[1.0,3.0,4.0])
@test_approx_eq_eps apb[1]  3.0 1e-12
@test_approx_eq_eps apb[2]  5.0 1e-12
@test_approx_eq_eps apb[3]  7.0 1e-12

# eraPpsp
apsb = eraPpsp([2.0,2.0,3.0],5.0,[1.0,3.0,4.0])
@test_approx_eq_eps apsb[1]  7.0 1e-12
@test_approx_eq_eps apsb[2]  17.0 1e-12
@test_approx_eq_eps apsb[3]  23.0 1e-12

# eraPr00
dpsipr, depspr = eraPr00(2400000.5, 53736.)
@test_approx_eq_eps dpsipr  -0.8716465172668347629e-7 1e-22
@test_approx_eq_eps depspr  -0.7342018386722813087e-8 1e-22

# eraPrec76
ep01 = 2400000.5
ep02 = 33282.0
ep11 = 2400000.5
ep12 = 51544.0
zeta, z, theta = eraPrec76(ep01, ep02, ep11, ep12)
@test_approx_eq_eps zeta   0.5588961642000161243e-2 1e-12
@test_approx_eq_eps z      0.5589922365870680624e-2 1e-12
@test_approx_eq_eps theta  0.4858945471687296760e-2 1e-12

# eraPv2p
p = eraPv2p([[0.3,1.2,-2.5];[-0.5,3.1,0.9]])
@test_approx_eq_eps p[1]   0.3 0.0
@test_approx_eq_eps p[2]   1.2 0.0
@test_approx_eq_eps p[3]  -2.5 0.0

# eraPv2s
pv = [[-0.4514964673880165,0.03093394277342585,0.05594668105108779];
      [1.292270850663260e-5,2.652814182060692e-6,2.568431853930293e-6]]
theta, phi, r, td, pd, rd = eraPv2s(pv)
@test_approx_eq_eps theta  3.073185307179586515 1e-12
@test_approx_eq_eps phi  0.1229999999999999992 1e-12
@test_approx_eq_eps r  0.4559999999999999757 1e-12
@test_approx_eq_eps td  -0.7800000000000000364e-5 1e-16
@test_approx_eq_eps pd  0.9010000000000001639e-5 1e-16
@test_approx_eq_eps rd  -0.1229999999999999832e-4 1e-16

# eraPvdpv
a = [[2.,2.,3.];[6.,0.,4.]]
b = [[1.,3.,4.];[0.,2.,8.]]
adb = eraPvdpv(a, b)
@test_approx_eq_eps adb[1]  20.0 1e-12
@test_approx_eq_eps adb[2]  50.0 1e-12

# eraPvm
pv = [[0.3,1.2,-2.5];[0.45,-0.25,1.1]]
r, s = eraPvm(pv)
@test_approx_eq_eps r  2.789265136196270604 1e-12
@test_approx_eq_eps s  1.214495780149111922 1e-12

# eraPvmpv
a = [[2.0,2.0,3.0];[5.0,6.0,3.0]]
b = [[1.0,3.0,4.0];[3.0,2.0,1.0]]
amb = eraPvmpv(a,b)
@test_approx_eq_eps amb[1]   1.0 1e-12
@test_approx_eq_eps amb[2]  -1.0 1e-12
@test_approx_eq_eps amb[3]  -1.0 1e-12
@test_approx_eq_eps amb[4]   2.0 1e-12
@test_approx_eq_eps amb[5]   4.0 1e-12
@test_approx_eq_eps amb[6]   2.0 1e-12

# eraPvppv
a = [[2.0,2.0,3.0];[5.0,6.0,3.0]]
b = [[1.0,3.0,4.0];[3.0,2.0,1.0]]
apb = eraPvppv(a,b)
@test_approx_eq_eps apb[1]  3.0 1e-12
@test_approx_eq_eps apb[2]  5.0 1e-12
@test_approx_eq_eps apb[3]  7.0 1e-12
@test_approx_eq_eps apb[4]  8.0 1e-12
@test_approx_eq_eps apb[5]  8.0 1e-12
@test_approx_eq_eps apb[6]  4.0 1e-12

# eraPvxpv
a = [[2.0,2.0,3.0];[6.0,0.0,4.0]]
b = [[1.0,3.0,4.0];[0.0,2.0,8.0]]
axb =eraPvxpv(a, b)
@test_approx_eq_eps axb[1]   -1.0 1e-12
@test_approx_eq_eps axb[2]   -5.0 1e-12
@test_approx_eq_eps axb[3]   4.0 1e-12
@test_approx_eq_eps axb[4]   -2.0 1e-12
@test_approx_eq_eps axb[5]   -36.0 1e-12
@test_approx_eq_eps axb[6]   22.0 1e-12

# eraPvstar
pv = [[126668.5912743160601,2136.792716839935195,-245251.2339876830091];
      [-0.4051854035740712739e-2,-0.6253919754866173866e-2,0.1189353719774107189e-1]]
ra, dec, pmr, pmd, px, rv = eraPvstar(pv)
@test_approx_eq_eps ra  0.1686756e-1 1e-12
@test_approx_eq_eps dec  -1.093989828 1e-12
@test_approx_eq_eps pmr  -0.178323516e-4 1e-16
@test_approx_eq_eps pmd  0.2336024047e-5 1e-16
@test_approx_eq_eps px  0.74723 1e-12
@test_approx_eq_eps rv  -21.6 1e-11

# eraPvtob
elong = 2.0
phi = 0.5
hm = 3000.0
xp = 1e-6
yp = -0.5e-6
sp = 1e-8
theta = 5.0
pv = eraPvtob(elong, phi, hm, xp, yp, sp, theta)
@test_approx_eq_eps pv[1]  4225081.367071159207 1e-5
@test_approx_eq_eps pv[2]  3681943.215856198144 1e-5
@test_approx_eq_eps pv[3]  3041149.399241260785 1e-5
@test_approx_eq_eps pv[4]  -268.4915389365998787 1e-9
@test_approx_eq_eps pv[5]  308.0977983288903123 1e-9
@test_approx_eq_eps pv[6]  0 1e-0

# eraPvu
dt = 2920.0
pv = [[126668.5912743160734,2136.792716839935565,-245251.2339876830229];
      [-0.4051854035740713039e-2,-0.6253919754866175788e-2,0.1189353719774107615e-1]]
upv = eraPvu(dt, pv)
@test_approx_eq_eps upv[1]  126656.7598605317105 1e-12
@test_approx_eq_eps upv[2]  2118.531271155726332 1e-12
@test_approx_eq_eps upv[3]  -245216.5048590656190 1e-12
@test_approx_eq_eps upv[4]  -0.4051854035740713039e-2 1e-12
@test_approx_eq_eps upv[5]  -0.6253919754866175788e-2 1e-12
@test_approx_eq_eps upv[6]  0.1189353719774107615e-1 1e-12

# eraPvup
dt = 2920.0
pv = [[126668.5912743160734,2136.792716839935565,-245251.2339876830229];
      [-0.4051854035740713039e-2,-0.6253919754866175788e-2,0.1189353719774107615e-1]]
p = eraPvup(dt, pv)
@test_approx_eq_eps p[1]   126656.7598605317105 1e-12
@test_approx_eq_eps p[2]   2118.531271155726332 1e-12
@test_approx_eq_eps p[3]  -245216.5048590656190 1e-12

# eraPxp
axb = eraPxp([2.0,2.0,3.0],[1.0,3.0,4.0])
@test_approx_eq_eps axb[1]  -1.0 1e-12
@test_approx_eq_eps axb[2]  -5.0 1e-12
@test_approx_eq_eps axb[3]   4.0 1e-12

# eraRefco
phpa = 800.0
tc = 10.0
rh = 0.9
wl = 0.4
refa, refb = eraRefco(phpa, tc, rh, wl)
@test_approx_eq_eps refa  0.2264949956241415009e-3 1e-15
@test_approx_eq_eps refb  -0.2598658261729343970e-6 1e-18

# eraRm2v
w = eraRm2v([[0.0,-0.8,-0.6];
             [0.8,-0.36,0.48];
             [0.6,0.48,-0.64]])
@test_approx_eq_eps w[1]   0.0 1e-12
@test_approx_eq_eps w[2]   1.413716694115406957 1e-12
@test_approx_eq_eps w[3]  -1.884955592153875943 1e-12

# eraRv2m
r = eraRv2m([0.0, 1.41371669, -1.88495559])
@test_approx_eq_eps r[1]  -0.7071067782221119905 1e-14
@test_approx_eq_eps r[2]  -0.5656854276809129651 1e-14
@test_approx_eq_eps r[3]  -0.4242640700104211225 1e-14
@test_approx_eq_eps r[4]   0.5656854276809129651 1e-14
@test_approx_eq_eps r[5]  -0.0925483394532274246 1e-14
@test_approx_eq_eps r[6]  -0.8194112531408833269 1e-14
@test_approx_eq_eps r[7]   0.4242640700104211225 1e-14
@test_approx_eq_eps r[8]  -0.8194112531408833269 1e-14
@test_approx_eq_eps r[9]   0.3854415612311154341 1e-14

# eraRx
phi = 0.3456789
r = [[2.0,3.0,2.0];
     [3.0,2.0,3.0];
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
r = [[2.0,3.0,2.0];
     [3.0,2.0,3.0];
     [3.0,4.0,5.0]]
p = [0.2,1.5,0.1]
rp = eraRxp(r, p)
@test_approx_eq_eps rp[1]  5.1 1e-12
@test_approx_eq_eps rp[2]  3.9 1e-12
@test_approx_eq_eps rp[3]  7.1 1e-12

# eraRxpv
r = [[2.0,3.0,2.0];
     [3.0,2.0,3.0];
     [3.0,4.0,5.0]]
pv = [[0.2,1.5,0.1];
      [1.5,0.2,0.1]]
rpv = eraRxpv(r, pv)
@test_approx_eq_eps rpv[1]  5.1 1e-12
@test_approx_eq_eps rpv[4]  3.8 1e-12
@test_approx_eq_eps rpv[2]  3.9 1e-12
@test_approx_eq_eps rpv[5]  5.2 1e-12
@test_approx_eq_eps rpv[3]  7.1 1e-12
@test_approx_eq_eps rpv[6]  5.8 1e-12

# eraRxr
a = [[2.0,3.0,2.0];
     [3.0,2.0,3.0];
     [3.0,4.0,5.0]]
b = [[1.0,2.0,2.0];
     [4.0,1.0,1.0];
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
r = [[2.0,3.0,2.0];
     [3.0,2.0,3.0];
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
r = [[2.0,3.0,2.0];
     [3.0,2.0,3.0];
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

# eraS2c
c = eraS2c(3.0123, -0.999)
@test_approx_eq_eps c[1]  -0.5366267667260523906 1e-12
@test_approx_eq_eps c[2]   0.0697711109765145365 1e-12
@test_approx_eq_eps c[3]  -0.8409302618566214041 1e-12

# eraS2p
p = eraS2p(-3.21, 0.123, 0.456)
@test_approx_eq_eps p[1]  -0.4514964673880165228 1e-12
@test_approx_eq_eps p[2]   0.0309339427734258688 1e-12
@test_approx_eq_eps p[3]   0.0559466810510877933 1e-12

# eraS2pv
pv = eraS2pv(-3.21, 0.123, 0.456, -7.8e-6, 9.01e-6, -1.23e-5)
@test_approx_eq_eps pv[1]  -0.4514964673880165228 1e-12
@test_approx_eq_eps pv[2]   0.0309339427734258688 1e-12
@test_approx_eq_eps pv[3]   0.0559466810510877933 1e-12
@test_approx_eq_eps pv[4]   0.1292270850663260170e-4 1e-16
@test_approx_eq_eps pv[5]   0.2652814182060691422e-5 1e-16
@test_approx_eq_eps pv[6]   0.2568431853930292259e-5 1e-16

# eraS2xpv
s1 = 2.0
s2 = 3.0
pv = [[0.3,1.2,-2.5];
      [0.5,2.3,-0.4]]
spv = eraS2xpv(s1, s2, pv)
@test_approx_eq_eps spv[1]   0.6 1e-12
@test_approx_eq_eps spv[2]   2.4 1e-12
@test_approx_eq_eps spv[3]  -5.0 1e-12
@test_approx_eq_eps spv[4]   1.5 1e-12
@test_approx_eq_eps spv[5]   6.9 1e-12
@test_approx_eq_eps spv[6]  -1.2 1e-12

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

# eraStarpm
ra1 =   0.01686756
dec1 = -1.093989828
pmr1 = -1.78323516e-5
pmd1 =  2.336024047e-6
px1 =   0.74723
rv1 = -21.6
ra2, dec2, pmr2, pmd2, px2, rv2 = eraStarpm(ra1, dec1, pmr1, pmd1, px1, rv1,
                                            2400000.5, 50083.0, 2400000.5, 53736.0)
@test_approx_eq_eps ra2  0.01668919069414242368 1e-13
@test_approx_eq_eps dec2  -1.093966454217127879 1e-13
@test_approx_eq_eps pmr2  -0.1783662682155932702e-4 1e-17
@test_approx_eq_eps pmd2  0.2338092915987603664e-5 1e-17
@test_approx_eq_eps px2  0.7473533835323493644 1e-13
@test_approx_eq_eps rv2  -21.59905170476860786 1e-11

# eraStarpv
ra =   0.01686756
dec = -1.093989828
pmr = -1.78323516e-5
pmd =  2.336024047e-6
px =   0.74723
rv = -21.6
pv = eraStarpv(ra, dec, pmr, pmd, px, rv)
@test_approx_eq_eps pv[1]  126668.5912743160601 1e-10
@test_approx_eq_eps pv[2]  2136.792716839935195 1e-12
@test_approx_eq_eps pv[3]  -245251.2339876830091 1e-10
@test_approx_eq_eps pv[4]  -0.4051854035740712739e-2 1e-13
@test_approx_eq_eps pv[5]  -0.6253919754866173866e-2 1e-15
@test_approx_eq_eps pv[6]  0.1189353719774107189e-1 1e-13

# eraSxp
s = 2.0
p = [0.3,1.2,-2.5]
sp = eraSxp(s, p)
@test_approx_eq_eps sp[1]   0.6 0.0
@test_approx_eq_eps sp[2]   2.4 0.0
@test_approx_eq_eps sp[3]  -5.0 0.0

# eraSxpv
s = 2.0
pv = [[0.3,1.2,-2.5];[0.5,3.2,-0.7]]
spv = eraSxpv(s, pv)
@test_approx_eq_eps spv[1]   0.6 0.0
@test_approx_eq_eps spv[2]   2.4 0.0
@test_approx_eq_eps spv[3]  -5.0 0.0
@test_approx_eq_eps spv[4]   1.0 0.0
@test_approx_eq_eps spv[5]   6.4 0.0
@test_approx_eq_eps spv[6]  -1.4 0.0

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

# eraTf2a
a = eraTf2a('+',4,58,20.2)
@test_approx_eq_eps a  1.301739278189537429 1e-12

# eraTf2d
d = eraTf2d('+',23,55,10.9)
@test_approx_eq_eps d  0.9966539351851851852 1e-12

# eraTr
r = [[2.0,3.0,2.0];
     [3.0,2.0,3.0];
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
r = [[2.0,3.0,2.0];
     [3.0,2.0,3.0];
     [3.0,4.0,5.0]]
p = [0.2,1.5,0.1]
trp = eraTrxp(r, p)
@test_approx_eq_eps trp[1]  5.2 1e-12
@test_approx_eq_eps trp[2]  4.0 1e-12
@test_approx_eq_eps trp[3]  5.4 1e-12

# eraTrxpv
r = [[2.0,3.0,2.0];
     [3.0,2.0,3.0];
     [3.0,4.0,5.0]]
pv = [[0.2,1.5,0.1];
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

# eraIcrs2g
dl, db = eraIcrs2g(5.9338074302227188048671087,-1.1784870613579944551540570)
@test_approx_eq_eps dl 5.5850536063818546461558 1e-14
@test_approx_eq_eps db -0.7853981633974483096157 1e14

# eraG2icrs
dr, dd = eraG2icrs(5.5850536063818546461558105, -0.7853981633974483096156608)
@test_approx_eq_eps dr 5.9338074302227188048671 1e-14
@test_approx_eq_eps dd -1.1784870613579944551541 1e14
