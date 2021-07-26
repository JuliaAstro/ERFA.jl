@testset "fad03" begin
    d = fad03(0.80)
    @test isapprox(d, 1.946709205396925672, atol = 1e-12)
end

@testset "fae03" begin
    e = fae03(0.80)
    @test isapprox(e, 1.744713738913081846, atol = 1e-12)
end

@testset "faf03" begin
    f = faf03(0.80)
    @test isapprox(f, 0.2597711366745499518, atol = 1e-12)
end

@testset "faju03" begin
    l = faju03(0.80)
    @test isapprox(l, 5.275711665202481138, atol = 1e-12)
end

@testset "fal03" begin
    l = fal03(0.80)
    @test isapprox(l, 5.132369751108684150, atol = 1e-12)
end

@testset "falp03" begin
    lp = falp03(0.80)
    @test isapprox(lp, 6.226797973505507345, atol = 1e-12)
end

@testset "fama03" begin
    l = fama03(0.80)
    @test isapprox(l, 3.275506840277781492, atol = 1e-12)
end

@testset "fame03" begin
    l = fame03(0.80)
    @test isapprox(l, 5.417338184297289661, atol = 1e-12)
end

@testset "fane03" begin
    l = fane03(0.80)
    @test isapprox(l, 2.079343830860413523, atol = 1e-12)
end

@testset "faom03" begin
    l = faom03(0.80)
    @test isapprox(l, -5.973618440951302183, atol = 1e-12)
end

@testset "fapa03" begin
    l = fapa03(0.80)
    @test isapprox(l, 0.1950884762240000000e-1, atol = 1e-12)
end

@testset "fasa03" begin
    l = fasa03(0.80)
    @test isapprox(l, 5.371574539440827046, atol = 1e-12)
end

@testset "faur03" begin
    l = faur03(0.80)
    @test isapprox(l, 5.180636450180413523, atol = 1e-12)
end

@testset "fave03" begin
    l = fave03(0.80)
    @test isapprox(l, 3.424900460533758000, atol = 1e-12)
end

@testset "fk425" begin
    r1950 = 0.07626899753879587532
    d1950 = -1.137405378399605780
    dr1950 = 0.1973749217849087460e-4
    dd1950 = 0.5659714913272723189e-5
    p1950 = 0.134
    v1950 = 8.7

    r2000, d2000, dr2000, dd2000, p2000, v2000 =
    ERFA.fk425(r1950, d1950, dr1950, dd1950, p1950, v1950)

    @test r2000 ≈ 0.08757989933556446040 atol=1e-14
    @test d2000 ≈ -1.132279113042091895 atol=1e-12
    @test dr2000 ≈ 0.1953670614474396139e-4 atol=1e-17
    @test dd2000 ≈ 0.5637686678659640164e-5 atol=1e-18
    @test p2000 ≈ 0.1339919950582767871 atol=1e-13
    @test v2000 ≈ 8.736999669183529069 atol=1e-12
end

@testset "fk45z" begin
    r1950 = 0.01602284975382960982
    d1950 = -0.1164347929099906024
    bepoch = 1954.677617625256806

    r2000, d2000 = ERFA.fk45z(r1950, d1950, bepoch)

    @test r2000 ≈ 0.02719295911606862303 atol=1e-15
    @test d2000 ≈ -0.1115766001565926892 atol=1e-13
end

@testset "fk524" begin
    r2000 = 0.8723503576487275595
    d2000 = -0.7517076365138887672
    dr2000 = 0.2019447755430472323e-4
    dd2000 = 0.3541563940505160433e-5
    p2000 = 0.1559
    v2000 = 86.87

    r1950, d1950, dr1950,dd1950, p1950, v1950 =
    ERFA.fk524(r2000, d2000, dr2000, dd2000, p2000, v2000)

    @test r1950 ≈ 0.8636359659799603487 atol=1e-13
    @test d1950 ≈ -0.7550281733160843059 atol=1e-13
    @test dr1950 ≈ 0.2023628192747172486e-4 atol=1e-17
    @test dd1950 ≈ 0.3624459754935334718e-5 atol=1e-18
    @test p1950 ≈ 0.1560079963299390241 atol=1e-13
    @test v1950 ≈ 86.79606353469163751 atol=1e-11
end

@testset "fk52h" begin
    r5  =  1.76779433
    d5  = -0.2917517103
    dr5 = -1.91851572e-7
    dd5 = -5.8468475e-6
    px5 =  0.379210
    rv5 = -7.6
    rh, dh, drh, ddh, pxh, rvh = fk52h(r5, d5, dr5, dd5, px5, rv5)
    @test isapprox(rh, 1.767794226299947632, atol = 1e-14)
    @test isapprox(dh, -0.2917516070530391757, atol = 1e-14)
    @test isapprox(drh, -0.1961874125605721270e-6, atol = 1e-19)
    @test isapprox(ddh, -0.58459905176693911e-5, atol = 1e-19)
    @test isapprox(pxh, 0.37921, atol = 1e-14)
    @test isapprox(rvh, -7.6000000940000254, atol = 1e-11)
end

@testset "fk54z" begin
    r2000 = 0.02719026625066316119
    d2000 = -0.1115815170738754813
    bepoch = 1954.677308160316374

    r1950, d1950, dr1950, dd1950 = ERFA.fk54z(r2000, d2000, bepoch)

    @test r1950 ≈ 0.01602015588390065476 atol=1e-14
    @test d1950 ≈ -0.1164397101110765346 atol=1e-13
    @test dr1950 ≈ -0.1175712648471090704e-7 atol=1e-20
    @test dd1950 ≈ 0.2108109051316431056e-7 atol=1e-20
end

@testset "fk5hz" begin
    r5 =  1.76779433
    d5 = -0.2917517103
    rh, dh = fk5hz(r5, d5, 2400000.5, 54479.0)
    @test isapprox(rh, 1.767794191464423978, atol = 1e-12)
    @test isapprox(dh, -0.2917516001679884419, atol = 1e-12)
end

@testset "fw2m" begin
    gamb = -0.2243387670997992368e-5
    phib =  0.4091014602391312982
    psi  = -0.9501954178013015092e-3
    eps  =  0.4091014316587367472
    r = fw2m(gamb, phib, psi, eps)
    @test isapprox(r[1,1], 0.9999995505176007047, atol = 1e-12)
    @test isapprox(r[1,2], 0.8695404617348192957e-3, atol = 1e-12)
    @test isapprox(r[1,3], 0.3779735201865582571e-3, atol = 1e-12)
    @test isapprox(r[2,1], -0.8695404723772016038e-3, atol = 1e-12)
    @test isapprox(r[2,2], 0.9999996219496027161, atol = 1e-12)
    @test isapprox(r[2,3], -0.1361752496887100026e-6, atol = 1e-12)
    @test isapprox(r[3,1], -0.3779734957034082790e-3, atol = 1e-12)
    @test isapprox(r[3,2], -0.1924880848087615651e-6, atol = 1e-12)
    @test isapprox(r[3,3], 0.9999999285679971958, atol = 1e-12)
end

@testset "fw2xy" begin
    gamb = -0.2243387670997992368e-5
    phib =  0.4091014602391312982
    psi  = -0.9501954178013015092e-3
    eps  =  0.4091014316587367472
    x, y = fw2xy(gamb, phib, psi, eps)
    @test isapprox(x, -0.3779734957034082790e-3, atol = 1e-14)
    @test isapprox(y, -0.1924880848087615651e-6, atol = 1e-14)
end

@testset "fk5hip" begin
    r5h, s5h = fk5hip()
    @test isapprox(r5h[1], 0.9999999999999928638, atol = 1e-14)
    @test isapprox(r5h[2], 0.1110223351022919694e-6, atol = 1e-17)
    @test isapprox(r5h[3], 0.4411803962536558154e-7, atol = 1e-17)
    @test isapprox(r5h[4], -0.1110223308458746430e-6, atol = 1e-17)
    @test isapprox(r5h[5], 0.9999999999999891830, atol = 1e-14)
    @test isapprox(r5h[6], -0.9647792498984142358e-7, atol = 1e-17)
    @test isapprox(r5h[7], -0.4411805033656962252e-7, atol = 1e-17)
    @test isapprox(r5h[8], 0.9647792009175314354e-7, atol = 1e-17)
    @test isapprox(r5h[9], 0.9999999999999943728, atol = 1e-14)
    @test isapprox(s5h[1], -0.1454441043328607981e-8, atol = 1e-17)
    @test isapprox(s5h[2], 0.2908882086657215962e-8, atol = 1e-17)
    @test isapprox(s5h[3], 0.3393695767766751955e-8, atol = 1e-17)
end

