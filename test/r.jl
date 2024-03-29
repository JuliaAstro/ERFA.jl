@testset "refco" begin
    phpa = 800.0
    tc = 10.0
    rh = 0.9
    wl = 0.4
    refa, refb = refco(phpa, tc, rh, wl)
    @test refa ≈ 0.2264949956241415009e-3 atol=1e-15
    @test refb ≈ -0.2598658261729343970e-6 atol=1e-18
end

@testset "rm2v" begin
    w = rm2v([0.0 -0.8 -0.6;
                   0.8 -0.36 0.48;
                   0.6 0.48 -0.64])
    @test w[1] ≈ 0.0 atol=1e-12
    @test w[2] ≈ 1.413716694115406957 atol=1e-12
    @test w[3] ≈ -1.884955592153875943 atol=1e-12
    @test_throws ArgumentError rm2v([0.0 -0.8 -0.6])
end

@testset "rv2m" begin
    r = rv2m([0.0, 1.41371669, -1.88495559])
    @test r[1,1] ≈ -0.7071067782221119905 atol=1e-14
    @test r[1,2] ≈ -0.5656854276809129651 atol=1e-14
    @test r[1,3] ≈ -0.4242640700104211225 atol=1e-14
    @test r[2,1] ≈ 0.5656854276809129651 atol=1e-14
    @test r[2,2] ≈ -0.0925483394532274246 atol=1e-14
    @test r[2,3] ≈ -0.8194112531408833269 atol=1e-14
    @test r[3,1] ≈ 0.4242640700104211225 atol=1e-14
    @test r[3,2] ≈ -0.8194112531408833269 atol=1e-14
    @test r[3,3] ≈ 0.3854415612311154341 atol=1e-14
    @test_throws ArgumentError rv2m([1.41371669, -1.88495559])
end

@testset "rx" begin
    phi = 0.3456789
    r = [2.0 3.0 2.0;
         3.0 2.0 3.0;
         3.0 4.0 5.0]
    r = rx(phi, r)
    @test r[1,1] ≈ 2.0 atol=0.0
    @test r[1,2] ≈ 3.0 atol=0.0
    @test r[1,3] ≈ 2.0 atol=0.0
    @test r[2,1] ≈ 3.839043388235612460 atol=1e-12
    @test r[2,2] ≈ 3.237033249594111899 atol=1e-12
    @test r[2,3] ≈ 4.516714379005982719 atol=1e-12
    @test r[3,1] ≈ 1.806030415924501684 atol=1e-12
    @test r[3,2] ≈ 3.085711545336372503 atol=1e-12
    @test r[3,3] ≈ 3.687721683977873065 atol=1e-12
    @test_throws ArgumentError rx(phi, r[1:2,:])
end

@testset "rxp" begin
    r = [2.0 3.0 2.0;
         3.0 2.0 3.0;
         3.0 4.0 5.0]
    p = [0.2,1.5,0.1]
    rp = ERFA._rxp(r, p)
    @test rp[1] ≈ 5.1 atol=1e-12
    @test rp[2] ≈ 3.9 atol=1e-12
    @test rp[3] ≈ 7.1 atol=1e-12
    rp_exp = r * p
    @testset for i in eachindex(rp, rp_exp)
        @test rp[i] ≈ rp_exp[i]
    end
    for i in eachindex(rp, rp_exp)
        @test rp[i] ≈ rp_exp[i]
    end
    @test_throws ArgumentError ERFA._rxp(r[1:2,:], p)
    @test_throws ArgumentError ERFA._rxp(r, p[1:2])
end

@testset "rxpv" begin
    r = [2.0 3.0 2.0;
         3.0 2.0 3.0;
         3.0 4.0 5.0]
    pv = [[0.2,1.5,0.1],
          [1.5,0.2,0.1]]
    rpv = ERFA._rxpv(r, pv)
    @test rpv[1][1] ≈ 5.1 atol=1e-12
    @test rpv[1][2] ≈ 3.9 atol=1e-12
    @test rpv[1][3] ≈ 7.1 atol=1e-12
    @test rpv[2][1] ≈ 3.8 atol=1e-12
    @test rpv[2][2] ≈ 5.2 atol=1e-12
    @test rpv[2][3] ≈ 5.8 atol=1e-12
    @test_throws ArgumentError ERFA._rxpv(r[1:2,:], pv)
    pve = [[1.5,0.1],
          [1.5,0.2,0.1]]
    @test_throws ArgumentError ERFA._rxpv(r, pve)
    rpj = [r * pv[1], r * pv[2]]
    for (v, vj) in zip(rpv, rpj)
        for i in eachindex(v, vj)
            @test v[i] ≈ vj[i]
        end
    end
end

@testset "rxr" begin
    a = [2.0 3.0 2.0;
         3.0 2.0 3.0;
         3.0 4.0 5.0]
    b = [1.0 2.0 2.0;
         4.0 1.0 1.0;
         3.0 0.0 1.0]
    atb = ERFA._rxr(a, b)
    @test atb[1,1] ≈ 20.0 atol=1e-12
    @test atb[1,2] ≈ 7.0 atol=1e-12
    @test atb[1,3] ≈ 9.0 atol=1e-12
    @test atb[2,1] ≈ 20.0 atol=1e-12
    @test atb[2,2] ≈ 8.0 atol=1e-12
    @test atb[2,3] ≈ 11.0 atol=1e-12
    @test atb[3,1] ≈ 34.0 atol=1e-12
    @test atb[3,2] ≈ 10.0 atol=1e-12
    @test atb[3,3] ≈ 15.0 atol=1e-12
    @test_throws ArgumentError ERFA._rxr(a[1:2,:], b)
    @test_throws ArgumentError ERFA._rxr(a, b[1:2,:])
    @test atb == a * b
end

@testset "ry" begin
    theta = 0.3456789
    r = [2.0 3.0 2.0;
         3.0 2.0 3.0;
         3.0 4.0 5.0]
    r = ry(theta, r)
    @test r[1,1] ≈ 0.8651847818978159930 atol=1e-12
    @test r[1,2] ≈ 1.467194920539316554 atol=1e-12
    @test r[1,3] ≈ 0.1875137911274457342 atol=1e-12
    @test r[2,1] ≈ 3 atol=1e-12
    @test r[2,2] ≈ 2 atol=1e-12
    @test r[2,3] ≈ 3 atol=1e-12
    @test r[3,1] ≈ 3.500207892850427330 atol=1e-12
    @test r[3,2] ≈ 4.779889022262298150 atol=1e-12
    @test r[3,3] ≈ 5.381899160903798712 atol=1e-12
    @test_throws ArgumentError ry(theta, r[1:2,:])
end

@testset "rz" begin
    psi = 0.3456789
    r = [2.0 3.0 2.0;
         3.0 2.0 3.0;
         3.0 4.0 5.0]
    r = rz(psi, r)
    @test r[1,1] ≈ 2.898197754208926769 atol=1e-12
    @test r[1,2] ≈ 3.500207892850427330 atol=1e-12
    @test r[1,3] ≈ 2.898197754208926769 atol=1e-12
    @test r[2,1] ≈ 2.144865911309686813 atol=1e-12
    @test r[2,2] ≈ 0.865184781897815993 atol=1e-12
    @test r[2,3] ≈ 2.144865911309686813 atol=1e-12
    @test r[3,1] ≈ 3.0 atol=1e-12
    @test r[3,2] ≈ 4.0 atol=1e-12
    @test r[3,3] ≈ 5.0 atol=1e-12
    @test_throws ArgumentError rz(psi, r[1:2,:])
end

