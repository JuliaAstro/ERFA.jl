@testset "xy06" begin
    x, y = xy06(2400000.5, 53736.0)
    @test x ≈ 0.5791308486706010975e-3 atol=1e-16
    @test y ≈ 0.4020579816732958141e-4 atol=1e-17
end

@testset "xys00a" begin
    x, y, s = xys00a(2400000.5, 53736.0)
    @test x ≈ 0.5791308472168152904e-3 atol=1e-16
    @test y ≈ 0.4020595661591500259e-4 atol=3e-17  # originally eps=1e-17; relaxed for 32-bit windows
    @test s ≈ -0.1220040848471549623e-7 atol=1e-20
end

@testset "xys00b" begin
    x, y, s = xys00b(2400000.5, 53736.0)
    @test x ≈ 0.5791301929950208873e-3 atol=1e-16
    @test y ≈ 0.4020553681373720832e-4 atol=1e-16
    @test s ≈ -0.1220027377285083189e-7 atol=1e-19
end

@testset "xys06a" begin
    x, y, s = xys06a(2400000.5, 53736.0)
    @test x ≈ 0.5791308482835292617e-3 atol=1e-16
    @test y ≈ 0.4020580099454020310e-4 atol=1e-15
    @test s ≈ -0.1220032294164579896e-7 atol=1e-19
end

