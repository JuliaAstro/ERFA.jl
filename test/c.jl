@testset "c2i00a" begin
    rc2i = ERFA.c2i00a(2400000.5, 53736.0)
    @test isapprox(rc2i[1,1], 0.9999998323037165557, atol = 1e-12)
    @test isapprox(rc2i[1,2], 0.5581526348992140183e-9, atol = 1e-12)
    @test isapprox(rc2i[1,3], -0.5791308477073443415e-3, atol = 1e-12)
    @test isapprox(rc2i[2,1], -0.2384266227870752452e-7, atol = 1e-12)
    @test isapprox(rc2i[2,2], 0.9999999991917405258, atol = 1e-12)
    @test isapprox(rc2i[2,3], -0.4020594955028209745e-4, atol = 1e-12)
    @test isapprox(rc2i[3,1], 0.5791308472168152904e-3, atol = 1e-12)
    @test isapprox(rc2i[3,2], 0.4020595661591500259e-4, atol = 1e-12)
    @test isapprox(rc2i[3,3], 0.9999998314954572304, atol = 1e-12)
end

@testset "c2i00b" begin
    rc2i = ERFA.c2i00b(2400000.5, 53736.0)
    @test isapprox(rc2i[1,1], 0.9999998323040954356, atol = 1e-12)
    @test isapprox(rc2i[1,2], 0.5581526349131823372e-9, atol = 1e-12)
    @test isapprox(rc2i[1,3], -0.5791301934855394005e-3, atol = 1e-12)
    @test isapprox(rc2i[2,1], -0.2384239285499175543e-7, atol = 1e-12)
    @test isapprox(rc2i[2,2], 0.9999999991917574043, atol = 1e-12)
    @test isapprox(rc2i[2,3], -0.4020552974819030066e-4, atol = 1e-12)
    @test isapprox(rc2i[3,1], 0.5791301929950208873e-3, atol = 1e-12)
    @test isapprox(rc2i[3,2], 0.4020553681373720832e-4, atol = 1e-12)
    @test isapprox(rc2i[3,3], 0.9999998314958529887, atol = 1e-12)
end

@testset "c2i06a" begin
    rc2i = ERFA.c2i06a(2400000.5, 53736.0)
    @test isapprox(rc2i[1,1], 0.9999998323037159379, atol = 1e-12)
    @test isapprox(rc2i[1,2], 0.5581121329587613787e-9, atol = 1e-12)
    @test isapprox(rc2i[1,3], -0.5791308487740529749e-3, atol = 1e-12)
    @test isapprox(rc2i[2,1], -0.2384253169452306581e-7, atol = 1e-12)
    @test isapprox(rc2i[2,2], 0.9999999991917467827, atol = 1e-12)
    @test isapprox(rc2i[2,3], -0.4020579392895682558e-4, atol = 1e-12)
    @test isapprox(rc2i[3,1], 0.5791308482835292617e-3, atol = 1e-12)
    @test isapprox(rc2i[3,2], 0.4020580099454020310e-4, atol = 1e-12)
    @test isapprox(rc2i[3,3], 0.9999998314954628695, atol = 1e-12)
end

@testset "c2ibpn" begin
    rbpn = [9.999962358680738e-1 -2.516417057665452e-3 -1.093569785342370e-3;
            2.516462370370876e-3 9.999968329010883e-1 4.006159587358310e-5;
            1.093465510215479e-3 -4.281337229063151e-5 9.999994012499173e-1]
    rc2i = ERFA.c2ibpn(2400000.5, 50123.9999, rbpn)
    @test isapprox(rc2i[1,1], 0.9999994021664089977, atol = 1e-12)
    @test isapprox(rc2i[1,2], -0.3869195948017503664e-8, atol = 1e-12)
    @test isapprox(rc2i[1,3], -0.1093465511383285076e-2, atol = 1e-12)
    @test isapprox(rc2i[2,1], 0.5068413965715446111e-7, atol = 1e-12)
    @test isapprox(rc2i[2,2], 0.9999999990835075686, atol = 1e-12)
    @test isapprox(rc2i[2,3], 0.4281334246452708915e-4, atol = 1e-12)
    @test isapprox(rc2i[3,1], 0.1093465510215479000e-2, atol = 1e-12)
    @test isapprox(rc2i[3,2], -0.4281337229063151000e-4, atol = 1e-12)
    @test isapprox(rc2i[3,3], 0.9999994012499173103, atol = 1e-12)
    @test_throws ArgumentError ERFA.c2ibpn(2400000.5, 50123.9999, rbpn[1:2,:])
end

@testset "c2s" begin
    t, p = ERFA.c2s([100.,-50.,25.])
    @test isapprox(t, -0.4636476090008061162, atol = 1e-15)
    @test isapprox(p, 0.2199879773954594463, atol = 1e-15)
    @test_throws ArgumentError ERFA.c2s([100.,-50.])
end

@testset "c2ixy" begin
    x = 0.5791308486706011000e-3
    y = 0.4020579816732961219e-4
    rc2i = ERFA.c2ixy(2400000.5, 53736., x, y)
    @test isapprox(rc2i[1,1], 0.9999998323037157138, atol = 1e-12)
    @test isapprox(rc2i[1,2], 0.5581526349032241205e-9, atol = 1e-12)
    @test isapprox(rc2i[1,3], -0.5791308491611263745e-3, atol = 1e-12)
    @test isapprox(rc2i[2,1], -0.2384257057469842953e-7, atol = 1e-12)
    @test isapprox(rc2i[2,2], 0.9999999991917468964, atol = 1e-12)
    @test isapprox(rc2i[2,3], -0.4020579110172324363e-4, atol = 1e-12)
    @test isapprox(rc2i[3,1], 0.5791308486706011000e-3, atol = 1e-12)
    @test isapprox(rc2i[3,2], 0.4020579816732961219e-4, atol = 1e-12)
    @test isapprox(rc2i[3,3], 0.9999998314954627590, atol = 1e-12)
end

@testset "c2ixys" begin
    x =  0.5791308486706011000e-3
    y =  0.4020579816732961219e-4
    s = -0.1220040848472271978e-7
    rc2i = ERFA.c2ixys(x, y, s)
    @test isapprox(rc2i[1,1], 0.9999998323037157138, atol = 1e-12)
    @test isapprox(rc2i[1,2], 0.5581984869168499149e-9, atol = 1e-12)
    @test isapprox(rc2i[1,3], -0.5791308491611282180e-3, atol = 1e-12)
    @test isapprox(rc2i[2,1], -0.2384261642670440317e-7, atol = 1e-12)
    @test isapprox(rc2i[2,2], 0.9999999991917468964, atol = 1e-12)
    @test isapprox(rc2i[2,3], -0.4020579110169668931e-4, atol = 1e-12)
    @test isapprox(rc2i[3,1], 0.5791308486706011000e-3, atol = 1e-12)
    @test isapprox(rc2i[3,2], 0.4020579816732961219e-4, atol = 1e-12)
    @test isapprox(rc2i[3,3], 0.9999998314954627590, atol = 1e-12)
end

@testset "c2t00a" begin
    tta = 2400000.5
    uta = 2400000.5
    ttb = 53736.0
    utb = 53736.0
    xp = 2.55060238e-7
    yp = 1.860359247e-6
    rc2t = ERFA.c2t00a(tta, ttb, uta, utb, xp, yp)
    @test isapprox(rc2t[1,1], -0.1810332128307182668, atol = 1e-12)
    @test isapprox(rc2t[1,2], 0.9834769806938457836, atol = 1e-12)
    @test isapprox(rc2t[1,3], 0.6555535638688341725e-4, atol = 1e-12)
    @test isapprox(rc2t[2,1], -0.9834768134135984552, atol = 1e-12)
    @test isapprox(rc2t[2,2], -0.1810332203649520727, atol = 1e-12)
    @test isapprox(rc2t[2,3], 0.5749801116141056317e-3, atol = 1e-12)
    @test isapprox(rc2t[3,1], 0.5773474014081406921e-3, atol = 1e-12)
    @test isapprox(rc2t[3,2], 0.3961832391770163647e-4, atol = 1e-12)
    @test isapprox(rc2t[3,3], 0.9999998325501692289, atol = 1e-12)
end

@testset "c2t00b" begin
    tta = 2400000.5
    uta = 2400000.5
    ttb = 53736.0
    utb = 53736.0
    xp = 2.55060238e-7
    yp = 1.860359247e-6
    rc2t = ERFA.c2t00b(tta, ttb, uta, utb, xp, yp)
    @test isapprox(rc2t[1,1], -0.1810332128439678965, atol = 1e-12)
    @test isapprox(rc2t[1,2], 0.9834769806913872359, atol = 1e-12)
    @test isapprox(rc2t[1,3], 0.6555565082458415611e-4, atol = 1e-12)
    @test isapprox(rc2t[2,1], -0.9834768134115435923, atol = 1e-12)
    @test isapprox(rc2t[2,2], -0.1810332203784001946, atol = 1e-12)
    @test isapprox(rc2t[2,3], 0.5749793922030017230e-3, atol = 1e-12)
    @test isapprox(rc2t[3,1], 0.5773467471863534901e-3, atol = 1e-12)
    @test isapprox(rc2t[3,2], 0.3961790411549945020e-4, atol = 1e-12)
    @test isapprox(rc2t[3,3], 0.9999998325505635738, atol = 1e-12)
end

@testset "c2t06a" begin
    tta = 2400000.5
    uta = 2400000.5
    ttb = 53736.0
    utb = 53736.0
    xp = 2.55060238e-7
    yp = 1.860359247e-6
    rc2t = ERFA.c2t06a(tta, ttb, uta, utb, xp, yp)
    @test isapprox(rc2t[1,1], -0.1810332128305897282, atol = 1e-12)
    @test isapprox(rc2t[1,2], 0.9834769806938592296, atol = 1e-12)
    @test isapprox(rc2t[1,3], 0.6555550962998436505e-4, atol = 1e-12)
    @test isapprox(rc2t[2,1], -0.9834768134136214897, atol = 1e-12)
    @test isapprox(rc2t[2,2], -0.1810332203649130832, atol = 1e-12)
    @test isapprox(rc2t[2,3], 0.5749800844905594110e-3, atol = 1e-12)
    @test isapprox(rc2t[3,1], 0.5773474024748545878e-3, atol = 1e-12)
    @test isapprox(rc2t[3,2], 0.3961816829632690581e-4, atol = 1e-12)
    @test isapprox(rc2t[3,3], 0.9999998325501747785, atol = 1e-12)
end

@testset "c2tcio" begin
    c = [0.9999998323037164738 0.5581526271714303683e-9 -0.5791308477073443903e-3;
         -0.2384266227524722273e-7 0.9999999991917404296 -0.4020594955030704125e-4;
         0.5791308472168153320e-3 .4020595661593994396e-4 0.9999998314954572365]
    era = 1.75283325530307
    p = [0.9999999999999674705 -0.1367174580728847031e-10 0.2550602379999972723e-6;
         0.1414624947957029721e-10 0.9999999999982694954 -0.1860359246998866338e-5;
         -0.2550602379741215275e-6 0.1860359247002413923e-5 0.9999999999982369658]
    rc2t = ERFA.c2tcio(c, era, p)
    @test isapprox(rc2t[1,1], -0.1810332128307110439, atol = 1e-12)
    @test isapprox(rc2t[1,2], 0.9834769806938470149, atol = 1e-12)
    @test isapprox(rc2t[1,3], 0.6555535638685466874e-4, atol = 1e-12)
    @test isapprox(rc2t[2,1], -0.9834768134135996657, atol = 1e-12)
    @test isapprox(rc2t[2,2], -0.1810332203649448367, atol = 1e-12)
    @test isapprox(rc2t[2,3], 0.5749801116141106528e-3, atol = 1e-12)
    @test isapprox(rc2t[3,1], 0.5773474014081407076e-3, atol = 1e-12)
    @test isapprox(rc2t[3,2], 0.3961832391772658944e-4, atol = 1e-12)
    @test isapprox(rc2t[3,3], 0.9999998325501691969, atol = 1e-12)
    @test_throws ArgumentError ERFA.c2tcio(c[1:2,:], era, p)
    @test_throws ArgumentError ERFA.c2tcio(c, era, p[1:2,:])
end

@testset "c2teqx" begin
    c = [0.9999989440476103608 -0.1332881761240011518e-2 -0.5790767434730085097e-3;
         0.1332858254308954453e-2 0.9999991109044505944 -0.4097782710401555759e-4;
         0.5791308472168153320e-3 0.4020595661593994396e-4 0.9999998314954572365]
    gst = 1.754166138040730516
    p = [0.9999999999999674705 -0.1367174580728847031e-10 0.2550602379999972723e-6;
         0.1414624947957029721e-10 0.9999999999982694954 -0.1860359246998866338e-5;
         -0.2550602379741215275e-6 0.1860359247002413923e-5 0.9999999999982369658]
    rc2t = ERFA.c2teqx(c, gst, p)
    @test isapprox(rc2t[1,1], -0.1810332128528685730, atol = 1e-12)
    @test isapprox(rc2t[1,2], 0.9834769806897685071, atol = 1e-12)
    @test isapprox(rc2t[1,3], 0.6555535639982634449e-4, atol = 1e-12)
    @test isapprox(rc2t[2,1], -0.9834768134095211257, atol = 1e-12)
    @test isapprox(rc2t[2,2], -0.1810332203871023800, atol = 1e-12)
    @test isapprox(rc2t[2,3], 0.5749801116126438962e-3, atol = 1e-12)
    @test isapprox(rc2t[3,1], 0.5773474014081539467e-3, atol = 1e-12)
    @test isapprox(rc2t[3,2], 0.3961832391768640871e-4, atol = 1e-12)
    @test isapprox(rc2t[3,3], 0.9999998325501691969, atol = 1e-12)
    @test_throws ArgumentError ERFA.c2teqx(c[1:2,:], gst, p)
    @test_throws ArgumentError ERFA.c2teqx(c, gst, p[1:2,:])
end

@testset "c2tpe" begin
    tta = 2400000.5
    uta = 2400000.5
    ttb = 53736.0
    utb = 53736.0
    deps =  0.4090789763356509900
    dpsi = -0.9630909107115582393e-5
    xp = 2.55060238e-7
    yp = 1.860359247e-6
    rc2t = ERFA.c2tpe(tta, ttb, uta, utb, dpsi, deps, xp, yp)
    @test isapprox(rc2t[1,1], -0.1813677995763029394, atol = 1e-12)
    @test isapprox(rc2t[1,2], 0.9023482206891683275, atol = 1e-12)
    @test isapprox(rc2t[1,3], -0.3909902938641085751, atol = 1e-12)
    @test isapprox(rc2t[2,1], -0.9834147641476804807, atol = 1e-12)
    @test isapprox(rc2t[2,2], -0.1659883635434995121, atol = 1e-12)
    @test isapprox(rc2t[2,3], 0.7309763898042819705e-1, atol = 1e-12)
    @test isapprox(rc2t[3,1], 0.1059685430673215247e-2, atol = 1e-12)
    @test isapprox(rc2t[3,2], 0.3977631855605078674, atol = 1e-12)
    @test isapprox(rc2t[3,3], 0.9174875068792735362, atol = 1e-12)
end

@testset "c2txy" begin
    tta = 2400000.5
    uta = 2400000.5
    ttb = 53736.0
    utb = 53736.0
    x = 0.5791308486706011000e-3
    y = 0.4020579816732961219e-4
    xp = 2.55060238e-7
    yp = 1.860359247e-6
    rc2t = ERFA.c2txy(tta, ttb, uta, utb, x, y, xp, yp)
    @test isapprox(rc2t[1,1], -0.1810332128306279253, atol = 1e-12)
    @test isapprox(rc2t[1,2], 0.9834769806938520084, atol = 1e-12)
    @test isapprox(rc2t[1,3], 0.6555551248057665829e-4, atol = 1e-12)
    @test isapprox(rc2t[2,1], -0.9834768134136142314, atol = 1e-12)
    @test isapprox(rc2t[2,2], -0.1810332203649529312, atol = 1e-12)
    @test isapprox(rc2t[2,3], 0.5749800843594139912e-3, atol = 1e-12)
    @test isapprox(rc2t[3,1], 0.5773474028619264494e-3, atol = 1e-12)
    @test isapprox(rc2t[3,2], 0.3961816546911624260e-4, atol = 1e-12)
    @test isapprox(rc2t[3,3], 0.9999998325501746670, atol = 1e-12)
end

@testset "cal2jd" begin
    dmj0, dmj = ERFA.cal2jd(2003, 6, 1)
    @test isapprox(dmj0, 2400000.5, atol = 1e-9)
    @test isapprox(dmj, 52791.0, atol = 1e-9)
end

@testset "cp" begin
    a = randn(3)
    b1 = ERFA._cp(a)
    b2 = copy(a)
    b1[1] = 1.0
    b2[1] = 1.0
    @test a != b1
    @test a != b2
    @test b1 == b2
end

@testset "cpv" begin
    a = [randn(3), randn(3)]
    b1 = ERFA._cpv(a)
    b2 = deepcopy(a)
    b1[1][1] = 1.0
    b2[1][1] = 1.0
    @test a != b1
    @test a != b2
    @test b1 == b2
end

@testset "cr" begin
    a = randn(3, 3)
    b1 = ERFA._cr(a)
    b2 = copy(a)
    b1[1,1] = 1.0
    b2[1,1] = 1.0
    @test a != b1
    @test a != b2
    @test b1 == b2
end

