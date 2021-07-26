@testset "jd2cal" begin
    y, m, d, fd = ERFA.jd2cal(2400000.5, 50123.9999)
    @test (y, m, d) == (1996, 2, 10)
    @test isapprox(fd, 0.9999, atol = 1e-7)
end

@testset "jdcalf" begin
    y, m, d, fd = ERFA.jdcalf(4, 2400000.5, 50123.9999)
    @test (y, m, d, fd) == (1996, 2, 10, 9999)
end

