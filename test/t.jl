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

# ERFA.tr
@testset "tr" begin
    r = [[2.0,3.0,2.0];
         [3.0,2.0,3.0];
         [3.0,4.0,5.0]]
    rt = ERFA.tr(r)
    @test isapprox(rt[1], 2.0, atol = 0.0)
    @test isapprox(rt[2], 3.0, atol = 0.0)
    @test isapprox(rt[3], 3.0, atol = 0.0)
    @test isapprox(rt[4], 3.0, atol = 0.0)
    @test isapprox(rt[5], 2.0, atol = 0.0)
    @test isapprox(rt[6], 4.0, atol = 0.0)
    @test isapprox(rt[7], 2.0, atol = 0.0)
    @test isapprox(rt[8], 3.0, atol = 0.0)
    @test isapprox(rt[9], 5.0, atol = 0.0)
end

# ERFA.trxp
@testset "trxp" begin
    r = [[2.0,3.0,2.0];
         [3.0,2.0,3.0];
         [3.0,4.0,5.0]]
    p = [0.2,1.5,0.1]
    trp = ERFA.trxp(r, p)
    @test isapprox(trp[1], 5.2, atol = 1e-12)
    @test isapprox(trp[2], 4.0, atol = 1e-12)
    @test isapprox(trp[3], 5.4, atol = 1e-12)
end

# ERFA.trxpv
@testset "trxpv" begin
    r = [[2.0,3.0,2.0];
         [3.0,2.0,3.0];
         [3.0,4.0,5.0]]
    pv = [[0.2,1.5,0.1];
          [1.5,0.2,0.1]]
    trpv = ERFA.trxpv(r, pv)
    @test isapprox(trpv[1], 5.2, atol = 1e-12)
    @test isapprox(trpv[2], 4.0, atol = 1e-12)
    @test isapprox(trpv[3], 5.4, atol = 1e-12)
    @test isapprox(trpv[4], 3.9, atol = 1e-12)
    @test isapprox(trpv[5], 5.3, atol = 1e-12)
    @test isapprox(trpv[6], 4.1, atol = 1e-12)
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

