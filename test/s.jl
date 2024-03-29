@testset "s00" begin
    x = 0.5791308486706011000e-3
    y = 0.4020579816732961219e-4
    s = s00(2400000.5, 53736.0, x, y)
    @test s ≈ -0.1220036263270905693e-7 atol=1e-18
end

@testset "s00a" begin
    s = s00a(2400000.5, 52541.0)
    @test s ≈ -0.1340684448919163584e-7 atol=1e-18
end

@testset "s00b" begin
    s = s00b(2400000.5, 52541.0)
    @test s ≈ -0.1340695782951026584e-7 atol=1e-18
end

@testset "s06" begin
    x = 0.5791308486706011000e-3
    y = 0.4020579816732961219e-4
    s = s06(2400000.5, 53736.0, x, y)
    @test s ≈ -0.1220032213076463117e-7 atol=1e-18
end

@testset "s06a" begin
    s = s06a(2400000.5, 52541.0)
    @test s ≈ -0.1340680437291812383e-7 atol=1e-18
end

@testset "s2c" begin
    c = s2c(3.0123, -0.999)
    @test c[1] ≈ -0.5366267667260523906 atol=1e-12
    @test c[2] ≈ 0.0697711109765145365 atol=1e-12
    @test c[3] ≈ -0.8409302618566214041 atol=1e-12
end

@testset "s2p" begin
    p = s2p(-3.21, 0.123, 0.456)
    @test p[1] ≈ -0.4514964673880165228 atol=1e-12
    @test p[2] ≈ 0.0309339427734258688 atol=1e-12
    @test p[3] ≈ 0.0559466810510877933 atol=1e-12
end

@testset "s2pv" begin
    pv = s2pv(-3.21, 0.123, 0.456, -7.8e-6, 9.01e-6, -1.23e-5)
    @test pv[1][1] ≈ -0.4514964673880165228 atol=1e-12
    @test pv[1][2] ≈ 0.0309339427734258688 atol=1e-12
    @test pv[1][3] ≈ 0.0559466810510877933 atol=1e-12
    @test pv[2][1] ≈ 0.1292270850663260170e-4 atol=1e-16
    @test pv[2][2] ≈ 0.2652814182060691422e-5 atol=1e-16
    @test pv[2][3] ≈ 0.2568431853930292259e-5 atol=1e-16
end

@testset "s2xpv" begin
    s1 = 2.0
    s2 = 3.0
    pv = [[0.3,1.2,-2.5],
          [0.5,2.3,-0.4]]
    spv = ERFA._s2xpv(s1, s2, pv)
    @test spv[1][1] ≈ 0.6 atol=1e-12
    @test spv[1][2] ≈ 2.4 atol=1e-12
    @test spv[1][3] ≈ -5.0 atol=1e-12
    @test spv[2][1] ≈ 1.5 atol=1e-12
    @test spv[2][2] ≈ 6.9 atol=1e-12
    @test spv[2][3] ≈ -1.2 atol=1e-12
    pve = [[1.2,-2.5],
          [0.5,2.3,-0.4]]
    @test_throws ArgumentError ERFA._s2xpv(s1, s2, pve)
    @test spv == [s1 .* pv[1], s2 .* pv[2]]
end

@testset "sepp" begin
    a = [1.,0.1,0.2]
    i = [10,1,2]
    b = [-3.,1e-3,0.2]
    ab = [1.,-3.,0.1,1e-3,0.2,0.2]
    expected_s = 2.860391919024660768
    s = sepp(a, b)
    @test s ≈ expected_s atol=1e-12
    s = sepp(i, b)
    @test s ≈ expected_s atol=1e-12
    s = sepp(b, i)
    @test s ≈ expected_s atol=1e-12
    s = sepp(ab[1:2:5], ab[2:2:6])
    @test s ≈ expected_s atol=1e-12
    s = @views sepp(ab[1:2:5], ab[2:2:6])
    @test s ≈ expected_s atol=1e-12
    @test_throws ArgumentError sepp(a[1:2], b)
    @test_throws ArgumentError sepp(a, b[1:2])
end

@testset "seps" begin
    s = seps(1., .1, .2, -3.)
    @test s ≈ 2.346722016996998842 atol=1e-14
end

@testset "sp00" begin
    s = sp00(2400000.5, 52541.0)
    @test s ≈ -0.6216698469981019309e-11 atol=1e-12
end

@testset "starpm" begin
    ra1 =   0.01686756
    dec1 = -1.093989828
    pmr1 = -1.78323516e-5
    pmd1 =  2.336024047e-6
    px1 =   0.74723
    rv1 = -21.6
    ra2, dec2, pmr2, pmd2, px2, rv2 = starpm(ra1, dec1, pmr1, pmd1, px1, rv1,
                                                  2400000.5, 50083.0, 2400000.5, 53736.0)
    @test ra2 ≈ 0.01668919069414256149 atol=1e-13
    @test dec2 ≈ -1.093966454217127897 atol=1e-13
    @test pmr2 ≈ -0.1783662682153176524e-4 atol=1e-17
    @test pmd2 ≈ 0.2338092915983989595e-5 atol=1e-17
    @test px2 ≈ 0.7473533835317719243 atol=1e-13
    @test rv2 ≈ -21.59905170476417175 atol=1e-11
end

@testset "starpv" begin
    ra =   0.01686756
    dec = -1.093989828
    pmr = -1.78323516e-5
    pmd =  2.336024047e-6
    px =   0.74723
    rv = -21.6
    pv = starpv(ra, dec, pmr, pmd, px, rv)
    @test pv[1][1] ≈ 126668.5912743160601 atol=1e-10
    @test pv[1][2] ≈ 2136.792716839935195 atol=1e-12
    @test pv[1][3] ≈ -245251.2339876830091 atol=1e-10
    @test pv[2][1] ≈ -0.4051854008955659551e-2 atol=1e-13
    @test pv[2][2] ≈ -0.6253919754414777970e-2 atol=1e-15
    @test pv[2][3] ≈ 0.1189353714588109341e-1 atol=1e-13
end

@testset "sxp" begin
    s = 2.0
    p = [0.3,1.2,-2.5]
    sp = ERFA._sxp(s, p)
    @test sp[1] ≈ 0.6 atol=0.0
    @test sp[2] ≈ 2.4 atol=0.0
    @test sp[3] ≈ -5.0 atol=0.0
    @test_throws ArgumentError ERFA._sxp(s, p[1:2])
    @test sp == s .* p
end

@testset "sxpv" begin
    s = 2.0
    pv = [[0.3,1.2,-2.5],[0.5,3.2,-0.7]]
    spv = ERFA._sxpv(s, pv)
    @test spv[1][1] ≈ 0.6 atol=0.0
    @test spv[1][2] ≈ 2.4 atol=0.0
    @test spv[1][3] ≈ -5.0 atol=0.0
    @test spv[2][1] ≈ 1.0 atol=0.0
    @test spv[2][2] ≈ 6.4 atol=0.0
    @test spv[2][3] ≈ -1.4 atol=0.0
    pve = [[1.2,-2.5],[0.5,3.2,-0.7]]
    @test_throws ArgumentError ERFA._sxpv(s, pve)
    @test spv == s .* pv
end

