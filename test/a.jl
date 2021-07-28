@testset "a2af" begin
    @test a2af(4, 2.345) == ('+', 134, 21, 30, 9706)
end

@testset "a2tf" begin
    @test a2tf(4, -3.01234) == ('-', 11, 30, 22, 6484)
end

@testset "ab" begin
    pnat = [-0.76321968546737951,-0.60869453983060384,-0.21676408580639883]
    v = [2.1044018893653786e-5,-8.9108923304429319e-5,-3.8633714797716569e-5]
    s = 0.99980921395708788
    bm1 = 0.99999999506209258
    ppr = ab(pnat, v, s, bm1)
    @test ppr[1] ≈ -0.7631631094219556269 atol=1e-12
    @test ppr[2] ≈ -0.6087553082505590832 atol=1e-12
    @test ppr[3] ≈ -0.2167926269368471279 atol=1e-12
    @test_throws ArgumentError ab(pnat[1:2], v, s, bm1)
    @test_throws ArgumentError ab(pnat, v[1:2], s, bm1)
end

@testset "ae2hd" begin
   a = 5.5
   e = 1.1
   p = 0.7

   h, d = ERFA.ae2hd(a, e, p)

   @test h ≈ 0.5933291115507309663 atol=1e-14
   @test d ≈ 0.9613934761647817620 atol=1e-14
end

@testset "af2a" begin
    r = af2a('-', 45, 13, 27.2)
    @test r ≈ -0.7893115794313644842 atol=1e-15
    r = af2a('+', 45, 13, 27.2)
    @test r ≈ 0.7893115794313644842 atol=1e-15
    @test_throws ERFAException af2a('+', 361, 13, 27.2)
    @test_throws ERFAException af2a('+', 45, 60, 27.2)
    @test_throws ERFAException af2a('+', 45, 13, 60.2)
end

@testset "anp" begin
    r = ERFA._anp(-0.1)
    @test r ≈ 6.183185307179586477 atol=1e-15
    @test r ≈ mod2pi(-0.1) atol=1-15
end

@testset "anpm" begin
    r = anpm(-4.0)
    @test r ≈ 2.283185307179586477 atol=1e-15
end

@testset "apcg" begin
    date1 = 2456165.5
    date2 = 0.401182685
    ebpv = [[0.901310875,-0.417402664,-0.180982288],
            [0.00742727954,0.0140507459,0.00609045792]]
    ehp = [0.903358544,-0.415395237,-0.180084014]
    astrom = apcg(date1, date2, ebpv, ehp)
    @test astrom.pmt ≈ 12.65133794027378508 atol=1e-11
    @test astrom.eb[1] ≈ 0.901310875 atol=1e-12
    @test astrom.eb[2] ≈ -0.417402664 atol=1e-12
    @test astrom.eb[3] ≈ -0.180982288 atol=1e-12
    @test astrom.eh[1] ≈ 0.8940025429324143045 atol=1e-12
    @test astrom.eh[2] ≈ -0.4110930268679817955 atol=1e-12
    @test astrom.eh[3] ≈ -0.1782189004872870264 atol=1e-12
    @test astrom.em ≈ 1.010465295811013146 atol=1e-12
    @test astrom.v[1] ≈ 0.4289638913597693554e-4 atol=1e-16
    @test astrom.v[2] ≈ 0.8115034051581320575e-4 atol=1e-16
    @test astrom.v[3] ≈ 0.3517555136380563427e-4 atol=1e-16
    @test astrom.bm1 ≈ 0.9999999951686012981 atol=1e-12
    @test astrom.bpn[1,1] ≈ 1.0 atol=1e-10
    @test astrom.bpn[2,1] ≈ 0.0 atol=1e-10
    @test astrom.bpn[3,1] ≈ 0.0 atol=1e-10
    @test astrom.bpn[1,2] ≈ 0.0 atol=1e-10
    @test astrom.bpn[2,2] ≈ 1.0 atol=1e-10
    @test astrom.bpn[3,2] ≈ 0.0 atol=1e-10
    @test astrom.bpn[1,3] ≈ 0.0 atol=1e-10
    @test astrom.bpn[2,3] ≈ 0.0 atol=1e-10
    @test astrom.bpn[3,3] ≈ 1.0 atol=1e-10
    ebpve = [[-0.417402664,-0.180982288],
            [0.00742727954,0.0140507459,0.00609045792]]
    @test_throws ArgumentError apcg(date1, date2, ebpve, ehp)
    @test_throws ArgumentError apcg(date1, date2, ebpv, ehp[1:2])
end

@testset "apcg13" begin
    date1 = 2456165.5
    date2 = 0.401182685
    astrom = apcg13(date1, date2)
    @test astrom.pmt ≈ 12.65133794027378508 atol=1e-12
    @test astrom.eb[1] ≈ 0.9013108747340644755 atol=1e-12
    @test astrom.eb[2] ≈ -0.4174026640406119957 atol=1e-12
    @test astrom.eb[3] ≈ -0.1809822877867817771 atol=1e-12
    @test astrom.eh[1] ≈ 0.8940025429255499549 atol=1e-12
    @test astrom.eh[2] ≈ -0.4110930268331896318 atol=1e-12
    @test astrom.eh[3] ≈ -0.1782189006019749850 atol=1e-12
    @test astrom.em ≈ 1.010465295964664178 atol=1e-12
    @test astrom.v[1] ≈ 0.4289638912941341125e-4 atol=1e-16
    @test astrom.v[2] ≈ 0.8115034032405042132e-4 atol=1e-16
    @test astrom.v[3] ≈ 0.3517555135536470279e-4 atol=1e-16
    @test astrom.bm1 ≈ 0.9999999951686013142 atol=1e-12
    @test astrom.bpn[1,1] ≈ 1.0 atol=1e-10
    @test astrom.bpn[2,1] ≈ 0.0 atol=1e-10
    @test astrom.bpn[3,1] ≈ 0.0 atol=1e-10
    @test astrom.bpn[1,2] ≈ 0.0 atol=1e-10
    @test astrom.bpn[2,2] ≈ 1.0 atol=1e-10
    @test astrom.bpn[3,2] ≈ 0.0 atol=1e-10
    @test astrom.bpn[1,3] ≈ 0.0 atol=1e-10
    @test astrom.bpn[2,3] ≈ 0.0 atol=1e-10
    @test astrom.bpn[3,3] ≈ 1.0 atol=1e-10
end

@testset "apci" begin
    date1 = 2456165.5
    date2 = 0.401182685
    ebpv = [[0.901310875,-0.417402664,-0.180982288],
            [0.00742727954,0.0140507459,0.00609045792]]
    ehp = [0.903358544,-0.415395237,-0.180084014]
    x =  0.0013122272
    y = -2.92808623e-5
    s =  3.05749468e-8
    astrom = apci(date1, date2, ebpv, ehp, x, y, s)
    @test astrom.pmt ≈ 12.65133794027378508 atol=1e-11
    @test astrom.eb[1] ≈ 0.901310875 atol=1e-12
    @test astrom.eb[2] ≈ -0.417402664 atol=1e-12
    @test astrom.eb[3] ≈ -0.180982288 atol=1e-12
    @test astrom.eh[1] ≈ 0.8940025429324143045 atol=1e-12
    @test astrom.eh[2] ≈ -0.4110930268679817955 atol=1e-12
    @test astrom.eh[3] ≈ -0.1782189004872870264 atol=1e-12
    @test astrom.em ≈ 1.010465295811013146 atol=1e-12
    @test astrom.v[1] ≈ 0.4289638913597693554e-4 atol=1e-16
    @test astrom.v[2] ≈ 0.8115034051581320575e-4 atol=1e-16
    @test astrom.v[3] ≈ 0.3517555136380563427e-4 atol=1e-16
    @test astrom.bm1 ≈ 0.9999999951686012981 atol=1e-12
    @test astrom.bpn[1,1] ≈ 0.9999991390295159156 atol=1e-12
    @test astrom.bpn[2,1] ≈ 0.4978650072505016932e-7 atol=1e-12
    @test astrom.bpn[3,1] ≈ 0.1312227200000000000e-2 atol=1e-12
    @test astrom.bpn[1,2] ≈ -0.1136336653771609630e-7 atol=1e-12
    @test astrom.bpn[2,2] ≈ 0.9999999995713154868 atol=1e-12
    @test astrom.bpn[3,2] ≈ -0.2928086230000000000e-4 atol=1e-12
    @test astrom.bpn[1,3] ≈ -0.1312227200895260194e-2 atol=1e-12
    @test astrom.bpn[2,3] ≈ 0.2928082217872315680e-4 atol=1e-12
    @test astrom.bpn[3,3] ≈ 0.9999991386008323373 atol=1e-12
    ebpve = [[-0.417402664,-0.180982288],
            [0.00742727954,0.0140507459,0.00609045792]]
    @test_throws ArgumentError apci(date1, date2, ebpve, ehp, x, y, s)
    @test_throws ArgumentError apci(date1, date2, ebpv, ehp[1:2], x, y, s)
end

@testset "apci13" begin
    date1 = 2456165.5
    date2 = 0.401182685
    astrom, eo = apci13(date1, date2)
    @test astrom.pmt ≈ 12.65133794027378508 atol=1e-11
    @test astrom.eb[1] ≈ 0.9013108747340644755 atol=1e-12
    @test astrom.eb[2] ≈ -0.4174026640406119957 atol=1e-12
    @test astrom.eb[3] ≈ -0.1809822877867817771 atol=1e-12
    @test astrom.eh[1] ≈ 0.8940025429255499549 atol=1e-12
    @test astrom.eh[2] ≈ -0.4110930268331896318 atol=1e-12
    @test astrom.eh[3] ≈ -0.1782189006019749850 atol=1e-12
    @test astrom.em ≈ 1.010465295964664178 atol=1e-12
    @test astrom.v[1] ≈ 0.4289638912941341125e-4 atol=1e-16
    @test astrom.v[2] ≈ 0.8115034032405042132e-4 atol=1e-16
    @test astrom.v[3] ≈ 0.3517555135536470279e-4 atol=1e-16
    @test astrom.bm1 ≈ 0.9999999951686013142 atol=1e-12
    @test astrom.bpn[1,1] ≈ 0.9999992060376761710 atol=1e-12
    @test astrom.bpn[2,1] ≈ 0.4124244860106037157e-7 atol=1e-12
    @test astrom.bpn[3,1] ≈ 0.1260128571051709670e-2 atol=1e-12
    @test astrom.bpn[1,2] ≈ -0.1282291987222130690e-7 atol=1e-12
    @test astrom.bpn[2,2] ≈ 0.9999999997456835325 atol=1e-12
    @test astrom.bpn[3,2] ≈ -0.2255288829420524935e-4 atol=1e-12
    @test astrom.bpn[1,3] ≈ -0.1260128571661374559e-2 atol=1e-12
    @test astrom.bpn[2,3] ≈ 0.2255285422953395494e-4 atol=1e-12
    @test astrom.bpn[3,3] ≈ 0.9999992057833604343 atol=1e-12
    @test eo ≈ -0.2900618712657375647e-2 atol=1e-12
end

@testset "apco" begin
    date1 = 2456384.5
    date2 = 0.970031644
    ebpv = [[-0.974170438,-0.211520082,-0.0917583024],
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
    astrom = apco(date1, date2, ebpv, ehp, x, y, s, theta, elong, phi, hm, xp, yp, sp, refa, refb)
    @test astrom.pmt ≈ 13.25248468622587269 atol=1e-11
    @test astrom.eb[1] ≈ -0.9741827110630897003 atol=1e-12
    @test astrom.eb[2] ≈ -0.2115130190135014340 atol=1e-12
    @test astrom.eb[3] ≈ -0.09179840186968295686 atol=1e-12
    @test astrom.eh[1] ≈ -0.9736425571689670428 atol=1e-12
    @test astrom.eh[2] ≈ -0.2092452125848862201 atol=1e-12
    @test astrom.eh[3] ≈ -0.09075578152261439954 atol=1e-12
    @test astrom.em ≈ 0.9998233241710617934 atol=1e-12
    @test astrom.v[1] ≈ 0.2078704992916728762e-4 atol=1e-16
    @test astrom.v[2] ≈ -0.8955360107151952319e-4 atol=1e-16
    @test astrom.v[3] ≈ -0.3863338994288951082e-4 atol=1e-16
    @test astrom.bm1 ≈ 0.9999999950277561236 atol=1e-12
    @test astrom.bpn[1,1] ≈ 0.9999991390295159156 atol=1e-12
    @test astrom.bpn[2,1] ≈ 0.4978650072505016932e-7 atol=1e-12
    @test astrom.bpn[3,1] ≈ 0.1312227200000000000e-2 atol=1e-12
    @test astrom.bpn[1,2] ≈ -0.1136336653771609630e-7 atol=1e-12
    @test astrom.bpn[2,2] ≈ 0.9999999995713154868 atol=1e-12
    @test astrom.bpn[3,2] ≈ -0.2928086230000000000e-4 atol=1e-12
    @test astrom.bpn[1,3] ≈ -0.1312227200895260194e-2 atol=1e-12
    @test astrom.bpn[2,3] ≈ 0.2928082217872315680e-4 atol=1e-12
    @test astrom.bpn[3,3] ≈ 0.9999991386008323373 atol=1e-12
    @test astrom.along ≈ -0.5278008060295995734 atol=1e-12
    @test astrom.xpl ≈ 0.1133427418130752958e-5 atol=1e-17
    @test astrom.ypl ≈ 0.1453347595780646207e-5 atol=1e-17
    @test astrom.sphi ≈ -0.9440115679003211329 atol=1e-12
    @test astrom.cphi ≈ 0.3299123514971474711 atol=1e-12
    @test astrom.diurab ≈ 0 atol=1e-10
    @test astrom.eral ≈ 2.617608903970400427 atol=1e-12
    @test astrom.refa ≈ 0.2014187790000000000e-3 atol=1e-15
    @test astrom.refb ≈ -0.2361408310000000000e-6 atol=1e-18
    ebpve = [[-0.211520082,-0.0917583024],
            [0.00364365824,-0.0154287319,-0.00668922024]]
    @test_throws ArgumentError apco(date1, date2, ebpve, ehp, x, y, s, theta, elong,
                                         phi, hm, xp, yp, sp, refa, refb)
    @test_throws ArgumentError apco(date1, date2, ebpve, ehp[1:2], x, y, s, theta,
                                         elong, phi, hm, xp, yp, sp, refa, refb)
end

@testset "apco13" begin
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
    astrom, eo = apco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test astrom.pmt ≈ 13.25248468622475727 atol=1e-11
    @test astrom.eb[1] ≈ -0.9741827107321449445 atol=1e-12
    @test astrom.eb[2] ≈ -0.2115130190489386190 atol=1e-12
    @test astrom.eb[3] ≈ -0.09179840189515518726 atol=1e-12
    @test astrom.eh[1] ≈ -0.9736425572586866640 atol=1e-12
    @test astrom.eh[2] ≈ -0.2092452121602867431 atol=1e-12
    @test astrom.eh[3] ≈ -0.09075578153903832650 atol=1e-12
    @test astrom.em ≈ 0.9998233240914558422 atol=1e-12
    @test astrom.v[1] ≈ 0.2078704994520489246e-4 atol=1e-16
    @test astrom.v[2] ≈ -0.8955360133238868938e-4 atol=1e-16
    @test astrom.v[3] ≈ -0.3863338993055887398e-4 atol=1e-16
    @test astrom.bm1 ≈ 0.9999999950277561004 atol=1e-12
    @test astrom.bpn[1,1] ≈ 0.9999991390295147999 atol=1e-12
    @test astrom.bpn[2,1] ≈ 0.4978650075315529277e-7 atol=1e-12
    @test astrom.bpn[3,1] ≈ 0.001312227200850293372 atol=1e-12
    @test astrom.bpn[1,2] ≈ -0.1136336652812486604e-7 atol=1e-12
    @test astrom.bpn[2,2] ≈ 0.9999999995713154865 atol=1e-12
    @test astrom.bpn[3,2] ≈ -0.2928086230975367296e-4 atol=1e-12
    @test astrom.bpn[1,3] ≈ -0.001312227201745553566 atol=1e-12
    @test astrom.bpn[2,3] ≈ 0.2928082218847679162e-4 atol=1e-12
    @test astrom.bpn[3,3] ≈ 0.9999991386008312212 atol=1e-12
    @test astrom.along ≈ -0.5278008060295995733 atol=1e-12
    @test astrom.xpl ≈ 0.1133427418130752958e-5 atol=1e-17
    @test astrom.ypl ≈ 0.1453347595780646207e-5 atol=1e-17
    @test astrom.sphi ≈ -0.9440115679003211329 atol=1e-12
    @test astrom.cphi ≈ 0.3299123514971474711 atol=1e-12
    @test astrom.diurab ≈ 0 atol=1e-10
    @test astrom.eral ≈ 2.617608909189664000 atol=1e-12
    @test astrom.refa ≈ 0.2014187785940396921e-3 atol=1e-15
    @test astrom.refb ≈ -0.2361408314943696227e-6 atol=1e-18
    @test eo ≈ -0.003020548354802412839 atol=1e-14
    @test_logs (:warn,) apco13(2.5245935e6, 0.0, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test_throws ERFAException apco13(-1e9, 0.0, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
end

@testset "apcs" begin
    date1 = 2456384.5
    date2 = 0.970031644
    pv = [[-1836024.09,1056607.72,-5998795.26],
          [-77.0361767,-133.310856,0.0971855934]]
    ebpv = [[-0.974170438,-0.211520082,-0.0917583024],
            [0.00364365824,-0.0154287319,-0.00668922024]]
    ehp = [-0.973458265,-0.209215307,-0.0906996477]
    astrom = apcs(date1, date2, pv, ebpv, ehp)
    @test astrom.pmt ≈ 13.25248468622587269 atol=1e-11
    @test astrom.eb[1] ≈ -0.9741827110630456169 atol=1e-12
    @test astrom.eb[2] ≈ -0.2115130190136085494 atol=1e-12
    @test astrom.eb[3] ≈ -0.09179840186973175487 atol=1e-12
    @test astrom.eh[1] ≈ -0.9736425571689386099 atol=1e-12
    @test astrom.eh[2] ≈ -0.2092452125849967195 atol=1e-12
    @test astrom.eh[3] ≈ -0.09075578152266466572 atol=1e-12
    @test astrom.em ≈ 0.9998233241710457140 atol=1e-12
    @test astrom.v[1] ≈ 0.2078704993282685510e-4 atol=1e-16
    @test astrom.v[2] ≈ -0.8955360106989405683e-4 atol=1e-16
    @test astrom.v[3] ≈ -0.3863338994289409097e-4 atol=1e-16
    @test astrom.bm1 ≈ 0.9999999950277561237 atol=1e-12
    @test astrom.bpn[1,1] ≈ 1 atol=1e-10
    @test astrom.bpn[2,1] ≈ 0 atol=1e-10
    @test astrom.bpn[3,1] ≈ 0 atol=1e-10
    @test astrom.bpn[1,2] ≈ 0 atol=1e-10
    @test astrom.bpn[2,2] ≈ 1 atol=1e-10
    @test astrom.bpn[3,2] ≈ 0 atol=1e-10
    @test astrom.bpn[1,3] ≈ 0 atol=1e-10
    @test astrom.bpn[2,3] ≈ 0 atol=1e-10
    @test astrom.bpn[3,3] ≈ 1 atol=1e-10
    ebpve = [[-0.211520082,-0.0917583024],
            [0.00364365824,-0.0154287319,-0.00668922024]]
    @test_throws ArgumentError apcs(date1, date2, ebpve, ebpv, ehp)
    @test_throws ArgumentError apcs(date1, date2, pv, ebpve, ehp)
    @test_throws ArgumentError apcs(date1, date2, pv, ebpv, ehp[1:2])
end

@testset "apcs13" begin
    date1 = 2456165.5
    date2 = 0.401182685
    pv = [[-6241497.16,401346.896,-1251136.04],
          [-29.264597,-455.021831,0.0266151194]]
    astrom = apcs13(date1, date2, pv)
    @test astrom.pmt ≈ 12.65133794027378508 atol=1e-11
    @test astrom.eb[1] ≈ 0.9012691529023298391 atol=1e-12
    @test astrom.eb[2] ≈ -0.4173999812023068781 atol=1e-12
    @test astrom.eb[3] ≈ -0.1809906511146821008 atol=1e-12
    @test astrom.eh[1] ≈ 0.8939939101759726824 atol=1e-12
    @test astrom.eh[2] ≈ -0.4111053891734599955 atol=1e-12
    @test astrom.eh[3] ≈ -0.1782336880637689334 atol=1e-12
    @test astrom.em ≈ 1.010428384373318379 atol=1e-12
    @test astrom.v[1] ≈ 0.4279877294121697570e-4 atol=1e-16
    @test astrom.v[2] ≈ 0.7963255087052120678e-4 atol=1e-16
    @test astrom.v[3] ≈ 0.3517564013384691531e-4 atol=1e-16
    @test astrom.bm1 ≈ 0.9999999952947980978 atol=1e-12
    @test astrom.bpn[1,1] ≈ 1 atol=1e-10
    @test astrom.bpn[2,1] ≈ 0 atol=1e-10
    @test astrom.bpn[3,1] ≈ 0 atol=1e-10
    @test astrom.bpn[1,2] ≈ 0 atol=1e-10
    @test astrom.bpn[2,2] ≈ 1 atol=1e-10
    @test astrom.bpn[3,2] ≈ 0 atol=1e-10
    @test astrom.bpn[1,3] ≈ 0 atol=1e-10
    @test astrom.bpn[2,3] ≈ 0 atol=1e-10
    @test astrom.bpn[3,3] ≈ 1 atol=1e-10
    pve = [[401346.896,-1251136.04],
          [-29.264597,-455.021831,0.0266151194]]
    @test_throws ArgumentError apcs13(date1, date2, pve)
end

@testset "aper" begin
    theta = 5.678
    pmt = 0.
    eb = zeros(Cdouble, 3)
    eh = zeros(Cdouble, 3)
    em = 0.
    v = zeros(Cdouble, 3)
    bm1 = 0.
    bpn = zeros(Cdouble, 9)
    along = 1.234
    phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb = 0., 0., 0., 0., 0., 0., 0., 0., 0.
    astrom = ERFA.ASTROM(pmt, eb, eh, em, v, bm1, bpn, along,
                         phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb)
    astrom = aper(theta, astrom)
    @test astrom.eral ≈ 6.912000000000000000 atol=1e-12
end

@testset "aper13" begin
    ut11 = 2456165.5
    ut12 = 0.401182685
    pmt = 0.
    eb = zeros(Cdouble, 3)
    eh = zeros(Cdouble, 3)
    em = 0.
    v = zeros(Cdouble, 3)
    bm1 = 0.
    bpn = zeros(Cdouble, 9)
    along = 1.234
    phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb = 0., 0., 0., 0., 0., 0., 0., 0., 0.
    astrom = ERFA.ASTROM(pmt, eb, eh, em, v, bm1, bpn, along,
                         phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb)
    astrom = aper13(ut11, ut12, astrom)
    @test astrom.eral ≈ 3.316236661789694933 atol=1e-12
end

@testset "apio" begin
    sp = -3.01974337e-11
    theta = 3.14540971
    elong = -0.527800806
    phi = -1.2345856
    hm = 2738.0
    xp = 2.47230737e-7
    yp = 1.82640464e-6
    refa = 0.000201418779
    refb = -2.36140831e-7
    astrom = apio(sp, theta, elong, phi, hm, xp, yp, refa, refb)
    @test astrom.along ≈ -0.5278008060295995734 atol=1e-12
    @test astrom.xpl ≈ 0.1133427418130752958e-5 atol=1e-17
    @test astrom.ypl ≈ 0.1453347595780646207e-5 atol=1e-17
    @test astrom.sphi ≈ -0.9440115679003211329 atol=1e-12
    @test astrom.cphi ≈ 0.3299123514971474711 atol=1e-12
    @test astrom.diurab ≈ 0.5135843661699913529e-6 atol=1e-12
    @test astrom.eral ≈ 2.617608903970400427 atol=1e-12
    @test astrom.refa ≈ 0.2014187790000000000e-3 atol=1e-15
    @test astrom.refb ≈ -0.2361408310000000000e-6 atol=1e-18
end

@testset "apio13" begin
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
    astrom = apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test astrom.along ≈ -0.5278008060295995733 atol=1e-12
    @test astrom.xpl ≈ 0.1133427418130752958e-5 atol=1e-17
    @test astrom.ypl ≈ 0.1453347595780646207e-5 atol=1e-17
    @test astrom.sphi ≈ -0.9440115679003211329 atol=1e-12
    @test astrom.cphi ≈ 0.3299123514971474711 atol=1e-12
    @test astrom.diurab ≈ 0.5135843661699913529e-6 atol=1e-12
    @test astrom.eral ≈ 2.617608909189664000 atol=1e-12
    @test astrom.refa ≈ 0.2014187785940396921e-3 atol=1e-15
    @test astrom.refb ≈ -0.2361408314943696227e-6 atol=1e-18
    @test_logs (:warn,) apio13(2.5245935e6, 0.0, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test_throws ERFAException apio13(-1e9, 0.0, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
end

@testset "atcc13" begin
   rc = 2.71
   dc = 0.174
   pr = 1e-5
   pd = 5e-6
   px = 0.1
   rv = 55.0
   date1 = 2456165.5
   date2 = 0.401182685

   ra, da = atcc13(rc, dc, pr, pd, px, rv, date1, date2)

   @test ra ≈ 2.710126504531372384 atol=1e-12
   @test da ≈ 0.1740632537628350152 atol=1e-12
end

@testset "atccq" begin
   date1 = 2456165.5
   date2 = 0.401182685
   astrom, eo = apci13(date1, date2)
   rc = 2.71
   dc = 0.174
   pr = 1e-5
   pd = 5e-6
   px = 0.1
   rv = 55.0

   ra, da = atccq(rc, dc, pr, pd, px, rv, astrom)

   @test ra ≈ 2.710126504531372384 atol=1e-12
   @test da ≈ 0.1740632537628350152 atol=1e-12
end

@testset "atci13" begin
    rc = 2.71
    dc = 0.174
    pr = 1e-5
    pd = 5e-6
    px = 0.1
    rv = 55.0
    date1 = 2456165.5
    date2 = 0.401182685
    ri, di, eo = atci13(rc, dc, pr, pd, px, rv, date1, date2)
    @test ri ≈ 2.710121572969038991 atol=1e-12
    @test di ≈ 0.1729371367218230438 atol=1e-12
    @test eo ≈ -0.002900618712657375647 atol=1e-14
end

@testset "atciq" begin
    date1 = 2456165.5
    date2 = 0.401182685
    astrom, eo = apci13(date1, date2)
    rc = 2.71
    dc = 0.174
    pr = 1e-5
    pd = 5e-6
    px = 0.1
    rv = 55.0
    ri, di = atciq(rc, dc, pr, pd, px, rv, astrom)
    @test ri ≈ 2.710121572969038991 atol=1e-12
    @test di ≈ 0.1729371367218230438 atol=1e-12
end

@testset "atciqn" begin
    date1 = 2456165.5
    date2 = 0.401182685
    astrom, eo = apci13(date1, date2)
    rc = 2.71
    dc = 0.174
    pr = 1e-5
    pd = 5e-6
    px = 0.1
    rv = 55.0
    b1 = ERFA.LDBODY(0.00028574, 3e-10,
                     [[-7.81014427,-5.60956681,-1.98079819];
                      [0.0030723249,-0.00406995477,-0.00181335842]])
    b2 = ERFA.LDBODY(0.00095435, 3e-9,
                     [[0.738098796, 4.63658692,1.9693136];
                      [-0.00755816922, 0.00126913722, 0.000727999001]])
    b3 = ERFA.LDBODY(1.0, 6e-6,
                     [[-0.000712174377, -0.00230478303, -0.00105865966];
                      [6.29235213e-6, -3.30888387e-7, -2.96486623e-7]])
    b = [b1; b2; b3]
    ri, di = atciqn(rc, dc, pr, pd, px, rv, astrom, b)
    @test ri ≈ 2.710122008105325582 atol=1e-12
    @test di ≈ 0.1729371916491459122 atol=1e-12
end

@testset "atciqz" begin
    date1 = 2456165.5
    date2 = 0.401182685
    astrom, eo = apci13(date1, date2)
    rc = 2.71
    dc = 0.174
    ri, di = atciqz(rc, dc, astrom)
    @test ri ≈ 2.709994899247599271 atol=1e-12
    @test di ≈ 0.1728740720983623469 atol=1e-12
end

@testset "atco13" begin
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
    aob, zob, hob, dob, rob, eo = atco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong,
                                         phi, hm, xp, yp, phpa, tc, rh, wl)
    @test aob ≈ 0.9251774485485515207e-1 atol=1e-12
    @test zob ≈ 1.407661405256499357 atol=1e-12
    @test hob ≈ -0.9265154431529724692e-1 atol=1e-12
    @test dob ≈ 0.1716626560072526200 atol=1e-12
    @test rob ≈ 2.710260453504961012 atol=1e-12
    @test eo ≈ -0.003020548354802412839 atol=1e-14
    @test_logs (:warn,) atco13(rc, dc, pr, pd, px, rv, 2.5245935e6, 0.0, dut1, elong, phi,
                               hm, xp, yp, phpa, tc, rh, wl)
    @test_throws ERFAException atco13(rc, dc, pr, pd, px, rv, -1e9, 0.0, dut1, elong, phi,
                                      hm, xp, yp, phpa, tc, rh, wl)
end

@testset "atic13" begin
    ri = 2.710121572969038991
    di = 0.1729371367218230438
    date1 = 2456165.5
    date2 = 0.401182685
    rc, dc, eo = atic13(ri, di, date1, date2)
    @test rc ≈ 2.710126504531374930 atol=1e-12
    @test dc ≈ 0.1740632537628342320 atol=1e-12
    @test eo ≈ -0.002900618712657375647 atol=1e-14
end

@testset "aticq" begin
    ri = 2.710121572969038991
    di = 0.1729371367218230438
    date1 = 2456165.5
    date2 = 0.401182685
    astrom, eo = apci13(date1, date2)
    rc, dc = aticq(ri, di, astrom)
    @test rc ≈ 2.710126504531374930 atol=1e-12
    @test dc ≈ 0.1740632537628342320 atol=1e-12
end

@testset "aticqn" begin
    date1 = 2456165.5
    date2 = 0.401182685
    astrom, eo = apci13(date1, date2)
    ri = 2.709994899247599271
    di = 0.1728740720983623469
    b1 = ERFA.LDBODY(0.00028574, 3e-10,
                     [[-7.81014427,-5.60956681,-1.98079819];
                      [0.0030723249,-0.00406995477,-0.00181335842]])
    b2 = ERFA.LDBODY(0.00095435, 3e-9,
                     [[0.738098796, 4.63658692,1.9693136];
                      [-0.00755816922, 0.00126913722, 0.000727999001]])
    b3 = ERFA.LDBODY(1.0, 6e-6,
                     [[-0.000712174377, -0.00230478303, -0.00105865966];
                      [6.29235213e-6, -3.30888387e-7, -2.96486623e-7]])
    b = [b1; b2; b3]
    rc, dc = aticqn(ri, di, astrom, b)
    @test rc ≈ 2.709999575032685412 atol=1e-12
    @test dc ≈ 0.1739999656317778034 atol=1e-12
end

@testset "atio13" begin
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
    aob, zob, hob, dob, rob = atio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test aob ≈ 0.9233952224895122499e-1 atol=1e-12
    @test zob ≈ 1.407758704513549991 atol=1e-12
    @test hob ≈ -0.9247619879881698140e-1 atol=1e-12
    @test dob ≈ 0.1717653435756234676 atol=1e-12
    @test rob ≈ 2.710085107988480746 atol=1e-12
    @test_logs (:warn,) atio13(ri, di, 2.5245935e6, 0.0, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test_throws ERFAException atio13(ri, di, -1e9, 0.0, dut1, elong, phi, hm, xp, yp, phpa, tc, rh, wl)
end

@testset "atioq" begin
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
    astrom = apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp,
                         phpa, tc, rh, wl)
    ri = 2.710121572969038991
    di = 0.1729371367218230438
    aob, zob, hob, dob, rob = atioq(ri, di, astrom)
    @test aob ≈ 0.9233952224895122499e-1 atol=1e-12
    @test zob ≈ 1.407758704513549991 atol=1e-12
    @test hob ≈ -0.9247619879881698140e-1 atol=1e-12
    @test dob ≈ 0.1717653435756234676 atol=1e-12
    @test rob ≈ 2.710085107988480746 atol=1e-12
end

@testset "atoc13" begin
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
    rc, dc = atoc13("r", ob1, ob2, utc1, utc2, dut1,
                         elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test rc ≈ 2.709956744659136129 atol=1e-12
    @test dc ≈ 0.1741696500898471362 atol=1e-12
    ob1 = -0.09247619879782006106
    ob2 = 0.1717653435758265198
    rc, dc = atoc13("h", ob1, ob2, utc1, utc2, dut1,
                         elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test rc ≈ 2.709956744659734086 atol=1e-12
    @test dc ≈ 0.1741696500898471362 atol=1e-12
    ob1 = 0.09233952224794989993
    ob2 = 1.407758704513722461
    rc, dc = atoc13("a", ob1, ob2, utc1, utc2, dut1,
                         elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test rc ≈ 2.709956744659734086 atol=1e-12
    @test dc ≈ 0.1741696500898471366 atol=1e-12
    rc1, dc1 = @test_logs (:warn,) atoc13("foo", ob1, ob2, utc1, utc2, dut1,
                                          elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test rc == rc1
    @test dc == dc1
    @test_logs (:warn,) atoc13("a", ob1, ob2, 2.5245935e6, 0.0, dut1, elong, phi, hm,
                               xp, yp, phpa, tc, rh, wl)
    @test_throws ERFAException atoc13("a", ob1, ob2, -1e9, 0.0, dut1, elong, phi, hm,
                                      xp, yp, phpa, tc, rh, wl)
end

@testset "atoi13" begin
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
    ri, di = atoi13("r", ob1, ob2, utc1, utc2, dut1,
                         elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test ri ≈ 2.710121574447540810 atol=1e-12
    @test di ≈ 0.1729371839116608778 atol=1e-12
    ob1 = -0.09247619879782006106
    ob2 = 0.1717653435758265198
    ri, di = atoi13("h", ob1, ob2, utc1, utc2, dut1,
                         elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test ri ≈ 2.710121574448138676 atol=1e-12
    @test di ≈ 0.1729371839116608778 atol=1e-12
    ob1 = 0.09233952224794989993
    ob2 = 1.407758704513722461
    ri, di = atoi13("a", ob1, ob2, utc1, utc2, dut1,
                         elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test ri ≈ 2.710121574448138676 atol=1e-12
    @test di ≈ 0.1729371839116608781 atol=1e-12
    @test_logs (:warn,) atoi13("a", ob1, ob2, 2.5245935e6, 0.0, dut1,
                               elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    @test_throws ERFAException atoi13("a", ob1, ob2, -1e9, 0.0, dut1,
                                      elong, phi, hm, xp, yp, phpa, tc, rh, wl)
end

@testset "atoiq" begin
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
    astrom = apio13(utc1, utc2, dut1,
                         elong, phi, hm, xp, yp, phpa, tc, rh, wl)
    ob1 = 2.710085107986886201
    ob2 = 0.1717653435758265198
    ri, di = atoiq("r", ob1, ob2, astrom)
    @test ri ≈ 2.710121574447540810 atol=1e-12
    @test di ≈ 0.17293718391166087785 atol=1e-12
    ob1 = -0.09247619879782006106
    ob2 = 0.1717653435758265198
    ri, di = atoiq("h", ob1, ob2, astrom)
    @test ri ≈ 2.710121574448138676 atol=1e-12
    @test di ≈ 0.1729371839116608778 atol=1e-12
    ob1 = 0.09233952224794989993
    ob2 = 1.407758704513722461
    ri, di = atoiq("a", ob1, ob2, astrom)
    @test ri ≈ 2.710121574448138676 atol=1e-12
    @test di ≈ 0.1729371839116608781 atol=1e-12
end

