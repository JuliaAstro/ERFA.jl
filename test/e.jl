@testset "ee00" begin
    epsa =  0.4090789763356509900
    dpsi = -0.9630909107115582393e-5
    ee = ee00(2400000.5, 53736.0, epsa, dpsi)
    @test isapprox(ee, -0.8834193235367965479e-5, atol = 1e-18)
end

@testset "ee00a" begin
    ee = ee00a(2400000.5, 53736.0)
    @test isapprox(ee, -0.8834192459222588227e-5, atol = 1e-18)
end

@testset "ee00b" begin
    ee = ee00b(2400000.5, 53736.0)
    @test isapprox(ee, -0.8835700060003032831e-5, atol = 1e-18)
end

@testset "ee06a" begin
    ee = ee06a(2400000.5, 53736.0)
    @test isapprox(ee, -0.8834195072043790156e-5, atol = 1e-15)
end

@testset "eect00" begin
    ct = eect00(2400000.5, 53736.0)
    @test isapprox(ct, 0.2046085004885125264e-8, atol = 1e-20)
end

@testset "eo06a" begin
    eo = eo06a(2400000.5, 53736.0)
    @test isapprox(eo, -0.1332882371941833644e-2, atol = 1e-15)
end

@testset "eors" begin
    r = [0.9999989440476103608 -0.1332881761240011518e-2 -0.5790767434730085097e-3;
         0.1332858254308954453e-2 0.9999991109044505944 -0.4097782710401555759e-4;
         0.5791308472168153320e-3 0.4020595661593994396e-4 0.9999998314954572365]
    s = -0.1220040848472271978e-7
    eo = eors(r, s)
    @test isapprox(eo, -0.1332882715130744606e-2, atol = 1e-15)
    @test_throws ArgumentError eors(r[1:2,:], s)
end

@testset "eform" begin
    a, f = eform(ERFA.WGS84)
    @test isapprox(a, 6378137.0, atol = 1e-10)
    @test isapprox(f, 0.0033528106647474807, atol = 1e-18)
    a, f = eform(ERFA.GRS80)
    @test isapprox(a, 6378137.0, atol = 1e-10)
    @test isapprox(f, 0.0033528106811823189, atol = 1e-18)
    a, f = eform(ERFA.WGS72)
    @test isapprox(a, 6378135.0, atol = 1e-10)
    @test isapprox(f, 0.0033527794541675049, atol = 1e-18)
end

@testset "epb" begin
    b = epb(2415019.8135, 30103.18648)
    @test isapprox(b, 1982.418424159278580, atol = 1e-12)
end

@testset "epb2jd" begin
    dj0, dj1 = epb2jd(1957.3)
    @test isapprox(dj0, 2400000.5, atol = 1e-9)
    @test isapprox(dj1, 35948.1915101513, atol = 1e-9)
end

@testset "epj" begin
    j = epj(2451545, -7392.5)
    @test isapprox(j, 1979.760438056125941, atol = 1e-12)
end

@testset "epj2jd" begin
    dj0, dj1 = epj2jd(1996.8)
    @test isapprox(dj0, 2400000.5, atol = 1e-9)
    @test isapprox(dj1, 50375.7, atol = 1e-9)
end

@testset "epv00" begin
    pvh, pvb = epv00(2400000.5, 53411.52501161)
    @test isapprox(pvh[1][1], -0.7757238809297706813, atol = 1e-14)
    @test isapprox(pvh[1][2], 0.5598052241363340596, atol = 1e-14)
    @test isapprox(pvh[1][3], 0.2426998466481686993, atol = 1e-14)
    @test isapprox(pvh[2][1], -0.1091891824147313846e-1, atol = 1e-15)
    @test isapprox(pvh[2][2], -0.1247187268440845008e-1, atol = 1e-15)
    @test isapprox(pvh[2][3], -0.5407569418065039061e-2, atol = 1e-15)
    @test isapprox(pvb[1][1], -0.7714104440491111971, atol = 1e-14)
    @test isapprox(pvb[1][2], 0.5598412061824171323, atol = 1e-14)
    @test isapprox(pvb[1][3], 0.2425996277722452400, atol = 1e-14)
    @test isapprox(pvb[2][1], -0.1091874268116823295e-1, atol = 1e-15)
    @test isapprox(pvb[2][2], -0.1246525461732861538e-1, atol = 1e-15)
    @test isapprox(pvb[2][3], -0.5404773180966231279e-2, atol = 1e-15)
end

@testset "eqeq94" begin
    ee = eqeq94(2400000.5, 41234.0)
    @test isapprox(ee, 0.5357758254609256894e-4, atol = 1e-17)
end

@testset "era00" begin
    era = era00(2400000.5, 54388.0)
    @test isapprox(era, 0.4022837240028158102, atol = 1e-12)
end

@testset "eceq06" begin
    dr, dd = eceq06(2456165.5, 0.401182685, 5.1, -0.9)
    @test isapprox(dr, 5.533459733613627767, atol = 1e-14)
    @test isapprox(dd, -1.246542932554480576, atol = 1e-14)
end

@testset "eqec06" begin
    dl, db = eqec06(1234.5, 2440000.5, 1.234, 0.987)
    @test isapprox(dl, 1.342509918994654619, atol = 1e-14)
    @test isapprox(db, 0.5926215259704608132, atol = 1e-14)
end

@testset "ecm06" begin
    rm = ecm06(2456165.5, 0.401182685)
    @test isapprox(rm[1,1], 0.9999952427708701137, atol = 1e-14)
    @test isapprox(rm[1,2], -0.2829062057663042347e-2, atol = 1e-14)
    @test isapprox(rm[1,3], -0.1229163741100017629e-2, atol = 1e-14)
    @test isapprox(rm[2,1], 0.3084546876908653562e-2, atol = 1e-14)
    @test isapprox(rm[2,2], 0.9174891871550392514, atol = 1e-14)
    @test isapprox(rm[2,3], 0.3977487611849338124, atol = 1e-14)
    @test isapprox(rm[3,1], 0.2488512951527405928e-5, atol = 1e-14)
    @test isapprox(rm[3,2], -0.3977506604161195467, atol = 1e-14)
    @test isapprox(rm[3,3], 0.9174935488232863071, atol = 1e-14)
end

