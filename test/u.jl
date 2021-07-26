@testset "ut1tai" begin
    a1, a2 = ERFA.ut1tai(2453750.5, 0.892104561, -32.6659)
    @test isapprox(a1, 2453750.5, atol = 1e-6)
    @test isapprox(a2, 0.8924826385462962963, atol = 1e-12)
end

@testset "ut1tt" begin
    a1, a2 = ERFA.ut1tt(2453750.5, 0.892104561, 64.8499)
    @test isapprox(a1, 2453750.5, atol = 1e-6)
    @test isapprox(a2, 0.8928551385462962963, atol = 1e-15)
end

@testset "ut1utc" begin
    a1, a2 = ERFA.ut1utc(2453750.5, 0.892104561, 0.3341)
    @test isapprox(a1, 2453750.5, atol = 1e-6)
    @test isapprox(a2, 0.8921006941018518519, atol = 1e-13)
end

@testset "utctai" begin
    u1, u2 = ERFA.utctai(2453750.5, 0.892100694)
    @test isapprox(u1, 2453750.5, atol = 1e-6)
    @test isapprox(u2, 0.8924826384444444444, atol = 1e-13)
end

@testset "utcut1" begin
    u1, u2 = ERFA.utcut1(2453750.5, 0.892100694, 0.3341)
    @test isapprox(u1, 2453750.5, atol = 1e-6)
    @test isapprox(u2, 0.8921045608981481481, atol = 1e-13)
end

