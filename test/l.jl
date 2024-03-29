@testset "ld" begin
    bm = 0.00028574
    p = [-0.763276255, -0.608633767, -0.216735543]
    q = [-0.763276255, -0.608633767, -0.216735543]
    e = [0.76700421, 0.605629598, 0.211937094]
    em = 8.91276983
    dlim = 3e-10
    p1 = ld(bm, p, q, e, em, dlim)
    @test p1[1] ≈ -0.7632762548968159627 atol=1e-12
    @test p1[2] ≈ -0.6086337670823762701 atol=1e-12
    @test p1[3] ≈ -0.2167355431320546947 atol=1e-12
    @test_throws ArgumentError ld(bm, p[1:2], q, e, em, dlim)
    @test_throws ArgumentError ld(bm, p, q[1:2], e, em, dlim)
    @test_throws ArgumentError ld(bm, p, q, e[1:2], em, dlim)
end

@testset "ldn" begin
    sc = [-0.763276255, -0.608633767, -0.216735543]
    ob = [-0.974170437, -0.2115201, -0.0917583114]
    pv1 = [-7.81014427,-5.60956681,-1.98079819,
           0.0030723249,-0.00406995477,-0.00181335842]
    pv2 = [0.738098796, 4.63658692,1.9693136,
           -0.00755816922, 0.00126913722, 0.000727999001]
    pv3 = [-0.000712174377, -0.00230478303, -0.00105865966,
           6.29235213e-6, -3.30888387e-7, -2.96486623e-7]
    b1 = ERFA.LDBODY(0.00028574, 3e-10, pv1)
    b2 = ERFA.LDBODY(0.00095435, 3e-9, pv2)
    b3 = ERFA.LDBODY(1.0, 6e-6, pv3)
    l = [b1, b2, b3]
    sn = ldn(l, ob, sc)
    @test sn[1] ≈ -0.7632762579693333866 atol=1e-12
    @test sn[2] ≈ -0.6086337636093002660 atol=1e-12
    @test sn[3] ≈ -0.2167355420646328159 atol=1e-12
    @test_throws ArgumentError ldn(l, ob[1:2], sc)
    @test_throws ArgumentError ldn(l, ob, sc[1:2])
end

@testset "ldsun" begin
    p = [-0.763276255, -0.608633767, -0.216735543]
    e = [-0.973644023, -0.20925523, -0.0907169552]
    em = 0.999809214
    p1 = ldsun(p, e, em)
    @test p1[1] ≈ -0.7632762580731413169 atol=1e-12
    @test p1[2] ≈ -0.6086337635262647900 atol=1e-12
    @test p1[3] ≈ -0.2167355419322321302 atol=1e-12
    @test_throws ArgumentError ldsun(p[1:2], e, em)
    @test_throws ArgumentError ldsun(p, e[1:2], em)
end

@testset "ltp" begin
    rp = ltp(1666.666)
    @test rp[1,1] ≈ 0.9967044141159213819 atol=1e-14
    @test rp[1,2] ≈ 0.7437801893193210840e-1 atol=1e-14
    @test rp[1,3] ≈ 0.3237624409345603401e-1 atol=1e-14
    @test rp[2,1] ≈ -0.7437802731819618167e-1 atol=1e-14
    @test rp[2,2] ≈ 0.9972293894454533070 atol=1e-14
    @test rp[2,3] ≈ -0.1205768842723593346e-2 atol=1e-14
    @test rp[3,1] ≈ -0.3237622482766575399e-1 atol=1e-14
    @test rp[3,2] ≈ -0.1206286039697609008e-2 atol=1e-14
    @test rp[3,3] ≈ 0.9994750246704010914 atol=1e-14
end

@testset "ltpb" begin
    rp = ltpb(1666.666)
    @test rp[1,1] ≈ 0.9967044167723271851 atol=1e-14
    @test rp[1,2] ≈ 0.7437794731203340345e-1 atol=1e-14
    @test rp[1,3] ≈ 0.3237632684841625547e-1 atol=1e-14
    @test rp[2,1] ≈ -0.7437795663437177152e-1 atol=1e-14
    @test rp[2,2] ≈ 0.9972293947500013666 atol=1e-14
    @test rp[2,3] ≈ -0.1205741865911243235e-2 atol=1e-14
    @test rp[3,1] ≈ -0.3237630543224664992e-1 atol=1e-14
    @test rp[3,2] ≈ -0.1206316791076485295e-2 atol=1e-14
    @test rp[3,3] ≈ 0.9994750220222438819 atol=1e-14
end

@testset "ltpecl" begin
    vec = ltpecl(-1500.0)
    @test vec[1] ≈ 0.4768625676477096525e-3 atol=1e-14
    @test vec[2] ≈ -0.4052259533091875112 atol=1e-14
    @test vec[3] ≈ 0.9142164401096448012 atol=1e-14
end

@testset "ltpequ" begin
    vec = ltpequ(-2500.0)
    @test vec[1] ≈ -0.3586652560237326659 atol=1e-14
    @test vec[2] ≈ -0.1996978910771128475 atol=1e-14
    @test vec[3] ≈ 0.9118552442250819624 atol=1e-14
end

@testset "ltecm" begin
    rm = ltecm(-3000.0)
    @test rm[1,1] ≈ 0.3564105644859788825 atol=1e-14
    @test rm[1,2] ≈ 0.8530575738617682284 atol=1e-14
    @test rm[1,3] ≈ 0.3811355207795060435 atol=1e-14
    @test rm[2,1] ≈ -0.9343283469640709942 atol=1e-14
    @test rm[2,2] ≈ 0.3247830597681745976 atol=1e-14
    @test rm[2,3] ≈ 0.1467872751535940865 atol=1e-14
    @test rm[3,1] ≈ 0.1431636191201167793e-2 atol=1e-14
    @test rm[3,2] ≈ -0.4084222566960599342 atol=1e-14
    @test rm[3,3] ≈ 0.9127919865189030899 atol=1e-14
end

@testset "lteceq" begin
    dr, dd = lteceq(2500.0, 1.5, 0.6)
    @test dr ≈ 1.275156021861921167 atol=1e-14
    @test dd ≈ 0.9966573543519204791 atol=1e-14
end

@testset "lteqec" begin
    dl, db = lteqec(-1500.0, 1.234, 0.987)
    @test dl ≈ 0.5039483649047114859 atol=1e-14
    @test db ≈ 0.5848534459726224882 atol=1e-14
end

