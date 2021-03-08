@testset "zp" begin
    @test iszero(ERFA._zp(randn(3)))
    @test iszero(fill!(randn(3), 0.0))
end

@testset "zpv" begin
    pv1 = [randn(3), randn(3)]
    pv2 = [randn(3), randn(3)]
    @test iszero(ERFA._zpv(pv1))
    @test iszero(fill!.(pv2, 0.0))
end

@testset "zr" begin
    r1 = randn(3, 3)
    r2 = randn(3, 3)
    @test iszero(ERFA._zr(r1))
    @test iszero(fill!(r2, 0.0))
end

