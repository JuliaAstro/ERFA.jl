@testset "gc2gd" begin
    xyz = [2.e6, 3.e6, 5.244e6]
    e, p, h = ERFA.gc2gd(1, xyz)
    @test isapprox(e, 0.98279372324732907, atol = 1e-14)
    @test isapprox(p, 0.97160184819075459, atol = 1e-14)
    @test isapprox(h, 331.41724614260599, atol = 1e-8)
    e, p, h = ERFA.gc2gd(2, xyz)
    @test isapprox(e, 0.98279372324732907, atol = 1e-14)
    @test isapprox(p, 0.97160184820607853, atol = 1e-14)
    @test isapprox(h, 331.41731754844348, atol = 1e-8)
    e, p, h = ERFA.gc2gd(3, xyz)
    @test isapprox(e, 0.98279372324732907, atol = 1e-14)
    @test isapprox(p, 0.97160181811015119, atol = 1e-14)
    @test isapprox(h, 333.27707261303181, atol = 1e-8)
    @test_throws ArgumentError ERFA.gc2gd(1, xyz[1:2])
end

@testset "gc2gde" begin
    a = 6378136.0
    f = 0.0033528
    xyz = [2e6, 3e6, 5.244e6]
    e, p, h = ERFA.gc2gde(a, f, xyz)
    @test isapprox(e, 0.98279372324732907, atol = 1e-14)
    @test isapprox(p, 0.97160183775704115, atol = 1e-14)
    @test isapprox(h, 332.36862495764397, atol = 1e-8)
    @test_throws ArgumentError ERFA.gc2gde(a, f, xyz[1:2])
end

@testset "gd2gc" begin
    e = 3.1
    p = -0.5
    h = 2500.0
    xyz = ERFA.gd2gc(1, e, p, h)
    @test isapprox(xyz[1], -5599000.5577049947, atol = 1e-7)
    @test isapprox(xyz[2], 233011.67223479203, atol = 1e-7)
    @test isapprox(xyz[3], -3040909.4706983363, atol = 1e-7)
    xyz = ERFA.gd2gc(2, e, p, h)
    @test isapprox(xyz[1], -5599000.5577260984, atol = 1e-7)
    @test isapprox(xyz[2], 233011.6722356703, atol = 1e-7)
    @test isapprox(xyz[3], -3040909.4706095476, atol = 1e-7)
    xyz = ERFA.gd2gc(3, e, p, h)
    @test isapprox(xyz[1], -5598998.7626301490, atol = 1e-7)
    @test isapprox(xyz[2], 233011.5975297822, atol = 1e-7)
    @test isapprox(xyz[3], -3040908.6861467111, atol = 1e-7)
end

@testset "gd2gce" begin
    a = 6378136.0
    f = 0.0033528
    e = 3.1
    p = -0.5
    h = 2500.0
    xyz = ERFA.gd2gce(a, f, e, p, h)
    @test isapprox(xyz[1], -5598999.6665116328, atol = 1e-7)
    @test isapprox(xyz[2], 233011.63514630572, atol = 1e-7)
    @test isapprox(xyz[3], -3040909.0517314132, atol = 1e-7)
end

@testset "gmst00" begin
    g = ERFA.gmst00(2400000.5, 53736.0, 2400000.5, 53736.0)
    @test isapprox(g, 1.754174972210740592, atol = 1e-14)
end

@testset "gmst06" begin
    g = ERFA.gmst06(2400000.5, 53736.0, 2400000.5, 53736.0)
    @test isapprox(g, 1.754174971870091203, atol = 1e-14)
end

@testset "gmst82" begin
    g = ERFA.gmst82(2400000.5, 53736.0)
    @test isapprox(g, 1.754174981860675096, atol = 1e-14)
end

@testset "gst00a" begin
    g = ERFA.gst00a(2400000.5, 53736.0, 2400000.5, 53736.0)
    @test isapprox(g, 1.754166138018281369, atol = 1e-14)
end

@testset "gst00b" begin
    g = ERFA.gst00b(2400000.5, 53736.0)
    @test isapprox(g, 1.754166136510680589, atol = 1e-14)
end

@testset "gst06" begin
    rnpb = [0.9999989440476103608 -0.1332881761240011518e-2 -0.5790767434730085097e-3;
            0.1332858254308954453e-2 0.9999991109044505944 -0.4097782710401555759e-4;
            0.5791308472168153320e-3 0.4020595661593994396e-4 0.9999998314954572365]
    g = ERFA.gst06(2400000.5, 53736.0, 2400000.5, 53736.0, rnpb)
    @test isapprox(g, 1.754166138018167568, atol = 1e-14)
    @test_throws ArgumentError ERFA.gst06(2400000.5, 53736.0, 2400000.5, 53736.0, rnpb[1:2,:])
end

@testset "gst06a" begin
    g = ERFA.gst06a(2400000.5, 53736.0, 2400000.5, 53736.0)
    @test isapprox(g, 1.754166137675019159, atol = 1e-14)
end

@testset "gst94" begin
    g = ERFA.gst94(2400000.5, 53736.0)
    @test isapprox(g, 1.754166136020645203, atol = 1e-14)
end

@testset "g2icrs" begin
    dr, dd = ERFA.g2icrs(5.5850536063818546461558105, -0.7853981633974483096156608)
    @test isapprox(dr, 5.9338074302227188048671, atol = 1e-14)
    @test isapprox(dd, -1.1784870613579944551541, atol = 1e14)
end

