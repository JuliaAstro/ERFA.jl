@testset "dat" begin
    d = dat(2003, 6, 1, 0.0)
    @test isapprox(d, 32.0, atol = 1e-9)
    d = dat(2008, 1, 17, 0.0)
    @test isapprox(d, 33.0, atol = 1e-9)
end

@testset "d2dtf" begin
    y, m, d, H, M, S, F = d2dtf("UTC", 5, 2400000.5, 49533.99999)
    @test (y, m, d, H, M, S, F) == (1994, 6, 30, 23, 59, 60, 13599)
end

@testset "d2tf" begin
    @test d2tf(4, -0.987654321) == ('-', 23, 42, 13, 3333)
end

@testset "dtdb" begin
    d = dtdb(2448939.5, 0.123, 0.76543, 5.0123, 5525.242, 3190.0)
    @test isapprox(d, -0.1280368005936998991e-2, atol = 1e-17)
end

@testset "dtf2d" begin
    jd1, jd2 = dtf2d("UTC", 1994, 6, 30, 23, 59, 60.13599)
    @test isapprox(jd1 + jd2, 2449534.49999, atol = 1e-6)
end

