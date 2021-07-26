@testset "h2fk5" begin
    rh  =  1.767794352
    dh  = -0.2917512594
    drh = -2.76413026e-6
    ddh = -5.92994449e-6
    pxh =  0.379210
    rvh = -7.6
    r5, d5, dr5, dd5, px5, rv5 = h2fk5(rh, dh, drh, ddh, pxh, rvh)
    @test r5 ≈ 1.767794455700065506 atol=1e-13
    @test d5 ≈ -0.2917513626469638890 atol=1e-13
    @test dr5 ≈ -0.27597945024511204e-5 atol=1e-18
    @test dd5 ≈ -0.59308014093262838e-5 atol=1e-18
    @test px5 ≈ 0.37921 atol=1e-13
    @test rv5 ≈ -7.6000001309071126 atol=1e-10
end

@testset "hd2ae" begin
    h = 1.1
    d = 1.2
    p = 0.3

    a, e = ERFA.hd2ae(h, d, p)

    @test a ≈ 5.916889243730066194 atol=1e-13
    @test e ≈ 0.4472186304990486228 atol=1e-14
end

@testset "hd2pa" begin
    h = 1.1
    d = 1.2
    p = 0.3

    q = ERFA.hd2pa(h, d, p)

    @test q ≈ 1.906227428001995580 atol=1e-13
end

@testset "hfk5z" begin
    rh =  1.767794352
    dh = -0.2917512594
    r5, d5, dr5, dd5 = hfk5z(rh, dh, 2400000.5, 54479.0)
    @test r5 ≈ 1.767794490535581026 atol=1e-13
    @test d5 ≈ -0.2917513695320114258 atol=1e-14
    @test dr5 ≈ 0.4335890983539243029e-8 atol=1e-22
    @test dd5 ≈ -0.8569648841237745902e-9 atol=1e-23
end

