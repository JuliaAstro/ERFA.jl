# ERFA.taitt
@testset "taitt" begin
    t1, t2 = ERFA.taitt(2453750.5, 0.892482639)
    @test isapprox(t1, 2453750.5, atol = 1e-6)
    @test isapprox(t2, 0.892855139, atol = 1e-12)
end

# ERFA.taiut1
@testset "taiut1" begin
    u1, u2 = ERFA.taiut1(2453750.5, 0.892482639, -32.6659)
    @test isapprox(u1, 2453750.5, atol = 1e-6)
    @test isapprox(u2, 0.8921045614537037037, atol = 1e-12)
end

# ERFA.taiutc
@testset "taiutc" begin
    u1, u2 = ERFA.taiutc(2453750.5, 0.892482639)
    @test isapprox(u1, 2453750.5, atol = 1e-6)
    @test isapprox(u2, 0.8921006945555555556, atol = 1e-12)
end

# ERFA.tcbtdb
@testset "tcbtdb" begin
    b1, b2 = ERFA.tcbtdb(2453750.5, 0.893019599)
    @test isapprox(b1, 2453750.5, atol = 1e-6)
    @test isapprox(b2, 0.8928551362746343397, atol = 1e-12)
end

# ERFA.tcgtt
@testset "tcgtt" begin
    t1, t2 = ERFA.tcgtt(2453750.5,  0.892862531)
    @test isapprox(t1, 2453750.5, atol = 1e-6)
    @test isapprox(t2, 0.8928551387488816828, atol = 1e-12)
end

# ERFA.tdbtcb
@testset "tdbtcb" begin
    b1, b2 = ERFA.tdbtcb(2453750.5, 0.892855137)
    @test isapprox(b1, 2453750.5, atol = 1e-6)
    @test isapprox(b2, 0.8930195997253656716, atol = 1e-12)
end

# ERFA.tdbtt
@testset "tdbtt" begin
    t1, t2 = ERFA.tdbtt(2453750.5,  0.892855137, -0.000201)
    @test isapprox(t1, 2453750.5, atol = 1e-6)
    @test isapprox(t2, 0.8928551393263888889, atol = 1e-12)
end

# ERFA.tf2a
@testset "tf2a" begin
    a = ERFA.tf2a('+', 4, 58, 20.2)
    @test isapprox(a, 1.301739278189537429, atol = 1e-12)
end

# ERFA.tf2d
@testset "tf2d" begin
    d = ERFA.tf2d('+', 23, 55, 10.9)
    @test isapprox(d, 0.9966539351851851852, atol = 1e-12)
end

# ERFA.tpors
@testset "tpors" begin
    xi = -0.03
    eta = 0.07
    ra = 1.3
    dec = 1.5

    status, az1, bz1, az2, bz2 = ERFA.tpors(xi, eta, ra, dec)

    @test az1 ≈ 1.736621577783208748 atol=1e-13
    @test bz1 ≈ 1.436736561844090323 atol=1e-13

    @test az2 ≈ 4.004971075806584490 atol=1e-13
    @test bz2 ≈ 1.565084088476417917 atol=1e-13

    @test status == 2
end

# ERFA.tporv
@testset "tporv" begin
    xi = -0.03
    eta = 0.07
    ra = 1.3
    dec = 1.5

    v = ERFA.s2c(ra, dec)

    status, vz1, vz2 = ERFA.tporv(xi, eta, v)

    @test vz1[1] ≈ -0.02206252822366888610 atol=1e-15
    @test vz1[2] ≈ 0.1318251060359645016 atol=1e-14
    @test vz1[3] ≈ 0.9910274397144543895 atol=1e-14

    @test vz2[1] ≈ -0.003712211763801968173 atol=1e-16
    @test vz2[2] ≈ -0.004341519956299836813 atol=1e-16
    @test vz2[3] ≈ 0.9999836852110587012 atol=1e-14

    @test status == 2
end

# ERFA.tpsts
@testset "tpsts" begin
    xi = -0.03
    eta = 0.07
    raz = 2.3
    decz = 1.5

    ra, dec = ERFA.tpsts(xi, eta, raz, decz)

    @test ra ≈ 0.7596127167359629775 atol=1e-14
    @test dec ≈ 1.540864645109263028 atol=1e-13
end

# ERFA.tpstv
@testset "tpstv" begin
    xi = -0.03
    eta = 0.07
    raz = 2.3
    decz = 1.5

    vz = ERFA.s2c(raz, decz)

    v = ERFA.tpstv(xi, eta, vz)

    @test v[1] ≈ 0.02170030454907376677 atol=1e-15
    @test v[2] ≈ 0.02060909590535367447 atol=1e-15
    @test v[3] ≈ 0.9995520806583523804 atol=1e-14
end

# ERFA.tpxes
@testset "tpxes" begin
    ra = 1.3
    dec = 1.55
    raz = 2.3
    decz = 1.5

    status, xi, eta = ERFA.tpxes(ra, dec, raz, decz)

    @test xi ≈ -0.01753200983236980595 atol=1e-15
    @test eta ≈ 0.05962940005778712891 atol=1e-15

    @test status == 0
end

# ERFA.tpxev
@testset "tpxev" begin
    ra = 1.3
    dec = 1.55
    raz = 2.3
    decz = 1.5
    v = ERFA.s2c(ra, dec)
    vz = ERFA.s2c(raz, decz)

    status, xi, eta = ERFA.tpxev(v, vz)

    @test xi ≈ -0.01753200983236980595 atol=1e-15
    @test eta ≈ 0.05962940005778712891 atol=1e-15

    @test status == 0
end

# ERFA.tr
@testset "tr" begin
    r = [2.0 3.0 2.0;
         3.0 2.0 3.0;
         3.0 4.0 5.0]
    rt = ERFA._tr(r)
    @test isapprox(rt[1,1], 2.0, atol = 0.0)
    @test isapprox(rt[1,2], 3.0, atol = 0.0)
    @test isapprox(rt[1,3], 3.0, atol = 0.0)
    @test isapprox(rt[2,1], 3.0, atol = 0.0)
    @test isapprox(rt[2,2], 2.0, atol = 0.0)
    @test isapprox(rt[2,3], 4.0, atol = 0.0)
    @test isapprox(rt[3,1], 2.0, atol = 0.0)
    @test isapprox(rt[3,2], 3.0, atol = 0.0)
    @test isapprox(rt[3,3], 5.0, atol = 0.0)
    @test_throws ArgumentError ERFA._tr(r[1:2,:])
    @test rt == r'
end

# ERFA.trxp
@testset "trxp" begin
    r = [2.0 3.0 2.0;
         3.0 2.0 3.0;
         3.0 4.0 5.0]
    p = [0.2,1.5,0.1]
    trp = ERFA._trxp(r, p)
    @test isapprox(trp[1], 5.2, atol = 1e-12)
    @test isapprox(trp[2], 4.0, atol = 1e-12)
    @test isapprox(trp[3], 5.4, atol = 1e-12)
    @test_throws ArgumentError ERFA._trxp(r[1:2,:], p)
    @test_throws ArgumentError ERFA._trxp(r, p[1:2])
    @test trp == r' * p
end

# ERFA.trxpv
@testset "trxpv" begin
    r = [2.0 3.0 2.0;
         3.0 2.0 3.0;
         3.0 4.0 5.0]
    pv = [[0.2,1.5,0.1],
          [1.5,0.2,0.1]]
    trpv = ERFA._trxpv(r, pv)
    @test isapprox(trpv[1][1], 5.2, atol = 1e-12)
    @test isapprox(trpv[1][2], 4.0, atol = 1e-12)
    @test isapprox(trpv[1][3], 5.4, atol = 1e-12)
    @test isapprox(trpv[2][1], 3.9, atol = 1e-12)
    @test isapprox(trpv[2][2], 5.3, atol = 1e-12)
    @test isapprox(trpv[2][3], 4.1, atol = 1e-12)
    pve = [[1.5,0.1],
          [1.5,0.2,0.1]]
    @test_throws ArgumentError ERFA._trxpv(r[1:2,:], pv)
    @test_throws ArgumentError ERFA._trxpv(r, pve)
    trpvj = [r' * pv[1], r' * pv[2]]
    for (v, vj) in zip(trpv, trpvj)
        for i in eachindex(v, vj)
            @test v[i] ≈ vj[i]
        end
    end
end

# ERFA.tttai
@testset "tttai" begin
    t1, t2 = ERFA.tttai(2453750.5, 0.892482639)
    @test isapprox(t1, 2453750.5, atol = 1e-6)
    @test isapprox(t2, 0.892110139, atol = 1e-12)
end

# ERFA.tttcg
@testset "tttcg" begin
    t1, t2 = ERFA.tttcg(2453750.5, 0.892482639)
    @test isapprox(t1, 2453750.5, atol = 1e-6)
    @test isapprox(t2, 0.8924900312508587113, atol = 1e-12)
end

# ERFA.tttdb
@testset "tttdb" begin
    t1, t2 = ERFA.tttdb(2453750.5, 0.892855139, -0.000201)
    @test isapprox(t1, 2453750.5, atol = 1e-6)
    @test isapprox(t2, 0.8928551366736111111, atol = 1e-12)
end

# ERFA.ttut1
@testset "ttut1" begin
    t1, t2 = ERFA.ttut1(2453750.5, 0.892855139, 64.8499)
    @test isapprox(t1, 2453750.5, atol = 1e-6)
    @test isapprox(t2, 0.8921045614537037037, atol = 1e-12)
end

