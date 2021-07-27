@testset "dat" begin
    d = dat(2003, 6, 1, 0.0)
    @test d ≈ 32.0 atol=1e-9
    d = dat(2008, 1, 17, 0.0)
    @test d ≈ 33.0 atol=1e-9
    @test_logs (:warn,) dat(2100, 1, 17, 0.0)
    @test_throws ERFAException dat(1e9, 1, 1, 0.0)
    @test_throws ERFAException dat(-1e9, 1, 1, 0.0)
    @test_throws ERFAException dat(2000, 13, 1, 0.0)
    @test_throws ERFAException dat(2000, 1, 32, 0.0)
    @test_throws ERFAException dat(2000, 1, 1, 3.0)
end

@testset "d2dtf" begin
    y, m, d, H, M, S, F = d2dtf("UTC", 5, 2400000.5, 49533.99999)
    @test (y, m, d, H, M, S, F) == (1994, 6, 30, 23, 59, 60, 13599)
    @test_logs (:warn,) d2dtf("UTC", 5, 0, 0)
    @test_throws ERFAException d2dtf("TAI", 5, 2e9, 0.0)
end

@testset "d2tf" begin
    @test d2tf(4, -0.987654321) == ('-', 23, 42, 13, 3333)
end

@testset "dtdb" begin
    d = dtdb(2448939.5, 0.123, 0.76543, 5.0123, 5525.242, 3190.0)
    @test d ≈ -0.1280368005936998991e-2 atol=1e-17
end

@testset "dtf2d" begin
    jd1, jd2 = dtf2d("UTC", 1994, 6, 30, 23, 59, 60.13599)
    @test jd1 + jd2 ≈ 2449534.49999 atol=1e-6
    @test_logs (:warn,) dtf2d("UTC", 1950, 6, 30, 23, 59, 59.13599)
    @test_logs (:warn,) dtf2d("UTC", 2021, 6, 30, 23, 59, 60.13599)
    @test_logs (:warn,) dtf2d("UTC", 1950, 6, 30, 23, 59, 60.13599)
    @test_throws ERFAException dtf2d("TAI", -2e9, 6, 30, 23, 59, 59.0)
    @test_throws ERFAException dtf2d("TAI", 2021, 13, 30, 23, 59, 59.0)
    @test_throws ERFAException dtf2d("TAI", 2021, 12, 32, 23, 59, 59.0)
    @test_throws ERFAException dtf2d("TAI", 2021, 12, 31, 24, 59, 59.0)
    @test_throws ERFAException dtf2d("TAI", 2021, 12, 31, 23, 60, 59.0)
    @test_throws ERFAException dtf2d("TAI", 2021, 12, 31, 23, 59, -59.0)
end

