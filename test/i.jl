using LinearAlgebra: I

@testset "icrs2g" begin
    dl, db = icrs2g(5.9338074302227188048671087, -1.1784870613579944551540570)
    @test dl ≈ 5.5850536063818546461558 atol=1e-14
    @test db ≈ -0.7853981633974483096157 atol=1e14
end

@testset "ir" begin
    r1 = ERFA._ir()
    r2 = Array{Float64}(I, 3, 3)
    @test r1 == r2
end
