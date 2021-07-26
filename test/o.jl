@testset "obl06" begin
    obl = ERFA.obl06(2400000.5, 54388.0)
    @test isapprox(obl, 0.4090749229387258204, atol = 1e-16)
end

@testset "obl80" begin
    obl = ERFA.obl80(2400000.5, 54388.0)
    @test isapprox(obl, 0.409075134764381621, atol = 1e-16)
end

