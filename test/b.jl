# ERFA.bi00
@testset "bi00" begin
    dp, de, dr = ERFA.bi00()
    @test isapprox(dp, -0.2025309152835086613e-6, atol = 1e-15)
    @test isapprox(de, -0.3306041454222147847e-7, atol = 1e-15)
    @test isapprox(dr, -0.7078279744199225506e-7, atol = 1e-15)
end

# ERFA.bp00
@testset "bp00" begin
    rb, rp, rbp = ERFA.bp00(2400000.5, 50123.9999)
    @test isapprox(rb[1], 0.9999999999999942498, atol = 1e-12)
    @test isapprox(rb[2], -0.7078279744199196626e-7, atol = 1e-16)
    @test isapprox(rb[3], 0.8056217146976134152e-7, atol = 1e-16)
    @test isapprox(rb[4], 0.7078279477857337206e-7, atol = 1e-16)
    @test isapprox(rb[5], 0.9999999999999969484, atol = 1e-12)
    @test isapprox(rb[6], 0.3306041454222136517e-7, atol = 1e-16)
    @test isapprox(rb[7], -0.8056217380986972157e-7, atol = 1e-16)
    @test isapprox(rb[8], -0.3306040883980552500e-7, atol = 1e-16)
    @test isapprox(rb[9], 0.9999999999999962084, atol = 1e-12)
    @test isapprox(rp[1], 0.9999995504864048241, atol = 1e-12)
    @test isapprox(rp[2], 0.8696113836207084411e-3, atol = 1e-14)
    @test isapprox(rp[3], 0.3778928813389333402e-3, atol = 1e-14)
    @test isapprox(rp[4], -0.8696113818227265968e-3, atol = 1e-14)
    @test isapprox(rp[5], 0.9999996218879365258, atol = 1e-12)
    @test isapprox(rp[6], -0.1690679263009242066e-6, atol = 1e-14)
    @test isapprox(rp[7], -0.3778928854764695214e-3, atol = 1e-14)
    @test isapprox(rp[8], -0.1595521004195286491e-6, atol = 1e-14)
    @test isapprox(rp[9], 0.9999999285984682756, atol = 1e-12)
    @test isapprox(rbp[1], 0.9999995505175087260, atol = 1e-12)
    @test isapprox(rbp[2], 0.8695405883617884705e-3, atol = 1e-14)
    @test isapprox(rbp[3], 0.3779734722239007105e-3, atol = 1e-14)
    @test isapprox(rbp[4], -0.8695405990410863719e-3, atol = 1e-14)
    @test isapprox(rbp[5], 0.9999996219494925900, atol = 1e-12)
    @test isapprox(rbp[6], -0.1360775820404982209e-6, atol = 1e-14)
    @test isapprox(rbp[7], -0.3779734476558184991e-3, atol = 1e-14)
    @test isapprox(rbp[8], -0.1925857585832024058e-6, atol = 1e-14)
    @test isapprox(rbp[9], 0.9999999285680153377, atol = 1e-12)
end

# ERFA.bp06
@testset "bp06" begin
    rb, rp, rbp = ERFA.bp06(2400000.5, 50123.9999)
    @test isapprox(rb[1], 0.9999999999999942497, atol = 1e-12)
    @test isapprox(rb[2], -0.7078368960971557145e-7, atol = 1e-14)
    @test isapprox(rb[3], 0.8056213977613185606e-7, atol = 1e-14)
    @test isapprox(rb[4], 0.7078368694637674333e-7, atol = 1e-14)
    @test isapprox(rb[5], 0.9999999999999969484, atol = 1e-12)
    @test isapprox(rb[6], 0.3305943742989134124e-7, atol = 1e-14)
    @test isapprox(rb[7], -0.8056214211620056792e-7, atol = 1e-14)
    @test isapprox(rb[8], -0.3305943172740586950e-7, atol = 1e-14)
    @test isapprox(rb[9], 0.9999999999999962084, atol = 1e-12)
    @test isapprox(rp[1], 0.9999995504864960278, atol = 1e-12)
    @test isapprox(rp[2], 0.8696112578855404832e-3, atol = 1e-14)
    @test isapprox(rp[3], 0.3778929293341390127e-3, atol = 1e-14)
    @test isapprox(rp[4], -0.8696112560510186244e-3, atol = 1e-14)
    @test isapprox(rp[5], 0.9999996218880458820, atol = 1e-12)
    @test isapprox(rp[6], -0.1691646168941896285e-6, atol = 1e-14)
    @test isapprox(rp[7], -0.3778929335557603418e-3, atol = 1e-14)
    @test isapprox(rp[8], -0.1594554040786495076e-6, atol = 1e-14)
    @test isapprox(rp[9], 0.9999999285984501222, atol = 1e-12)
    @test isapprox(rbp[1], 0.9999995505176007047, atol = 1e-12)
    @test isapprox(rbp[2], 0.8695404617348208406e-3, atol = 1e-14)
    @test isapprox(rbp[3], 0.3779735201865589104e-3, atol = 1e-14)
    @test isapprox(rbp[4], -0.8695404723772031414e-3, atol = 1e-14)
    @test isapprox(rbp[5], 0.9999996219496027161, atol = 1e-12)
    @test isapprox(rbp[6], -0.1361752497080270143e-6, atol = 1e-14)
    @test isapprox(rbp[7], -0.3779734957034089490e-3, atol = 1e-14)
    @test isapprox(rbp[8], -0.1924880847894457113e-6, atol = 1e-14)
    @test isapprox(rbp[9], 0.9999999285679971958, atol = 1e-12)
end

# ERFA.bpn2xy
@testset "bpn2xy" begin
    rbpn = [9.999962358680738e-1 -2.516417057665452e-3 -1.093569785342370e-3;
            2.516462370370876e-3  9.999968329010883e-1 4.006159587358310e-5;
            1.093465510215479e-3 -4.281337229063151e-5 9.999994012499173e-1]
    x, y = ERFA.bpn2xy(rbpn)
    @test isapprox(x, 1.093465510215479e-3, atol = 1e-12)
    @test isapprox(y, -4.281337229063151e-5, atol = 1e-12)
    @test_throws ArgumentError ERFA.bpn2xy(rbpn[1:2,:])
end

