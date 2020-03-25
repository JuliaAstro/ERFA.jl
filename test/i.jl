# ERFA.icrs2g
@testset "icrs2g" begin
    dl, db = ERFA.icrs2g(5.9338074302227188048671087, -1.1784870613579944551540570)
    @test isapprox(dl, 5.5850536063818546461558, atol = 1e-14)
    @test isapprox(db, -0.7853981633974483096157, atol = 1e14)
end

