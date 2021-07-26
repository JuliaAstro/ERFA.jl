@testset "num00a" begin
    rmatn = num00a(2400000.5, 53736.0)
    @test rmatn[1,1] ≈ 0.9999999999536227949 atol=1e-12
    @test rmatn[1,2] ≈ 0.8836238544090873336e-5 atol=1e-12
    @test rmatn[1,3] ≈ 0.3830835237722400669e-5 atol=1e-12
    @test rmatn[2,1] ≈ -0.8836082880798569274e-5 atol=1e-12
    @test rmatn[2,2] ≈ 0.9999999991354655028 atol=1e-12
    @test rmatn[2,3] ≈ -0.4063240865362499850e-4 atol=1e-12
    @test rmatn[3,1] ≈ -0.3831194272065995866e-5 atol=1e-12
    @test rmatn[3,2] ≈ 0.4063237480216291775e-4 atol=1e-12
    @test rmatn[3,3] ≈ 0.9999999991671660338 atol=1e-12
end

@testset "num00b" begin
    rmatn = num00b(2400000.5, 53736.0)
    @test rmatn[1,1] ≈ 0.9999999999536069682 atol=1e-12
    @test rmatn[1,2] ≈ 0.8837746144871248011e-5 atol=1e-12
    @test rmatn[1,3] ≈ 0.3831488838252202945e-5 atol=1e-12
    @test rmatn[2,1] ≈ -0.8837590456632304720e-5 atol=1e-12
    @test rmatn[2,2] ≈ 0.9999999991354692733 atol=1e-12
    @test rmatn[2,3] ≈ -0.4063198798559591654e-4 atol=1e-12
    @test rmatn[3,1] ≈ -0.3831847930134941271e-5 atol=1e-12
    @test rmatn[3,2] ≈ 0.4063195412258168380e-4 atol=1e-12
    @test rmatn[3,3] ≈ 0.9999999991671806225 atol=1e-12
end

@testset "num06a" begin
    rmatn = num06a(2400000.5, 53736.)
    @test rmatn[1,1] ≈ 0.9999999999536227668 atol=1e-12
    @test rmatn[1,2] ≈ 0.8836241998111535233e-5 atol=1e-12
    @test rmatn[1,3] ≈ 0.3830834608415287707e-5 atol=1e-12
    @test rmatn[2,1] ≈ -0.8836086334870740138e-5 atol=1e-12
    @test rmatn[2,2] ≈ 0.9999999991354657474 atol=1e-12
    @test rmatn[2,3] ≈ -0.4063240188248455065e-4 atol=1e-12
    @test rmatn[3,1] ≈ -0.3831193642839398128e-5 atol=1e-12
    @test rmatn[3,2] ≈ 0.4063236803101479770e-4 atol=1e-12
    @test rmatn[3,3] ≈ 0.9999999991671663114 atol=1e-12
end

@testset "numat" begin
    epsa =  0.4090789763356509900
    dpsi = -0.9630909107115582393e-5
    deps =  0.4063239174001678826e-4
    rmatn = numat(epsa, dpsi, deps)
    @test rmatn[1,1] ≈ 0.9999999999536227949 atol=1e-12
    @test rmatn[1,2] ≈ 0.8836239320236250577e-5 atol=1e-12
    @test rmatn[1,3] ≈ 0.3830833447458251908e-5 atol=1e-12
    @test rmatn[2,1] ≈ -0.8836083657016688588e-5 atol=1e-12
    @test rmatn[2,2] ≈ 0.9999999991354654959 atol=1e-12
    @test rmatn[2,3] ≈ -0.4063240865361857698e-4 atol=1e-12
    @test rmatn[3,1] ≈ -0.3831192481833385226e-5 atol=1e-12
    @test rmatn[3,2] ≈ 0.4063237480216934159e-4 atol=1e-12
    @test rmatn[3,3] ≈ 0.9999999991671660407 atol=1e-12
end

@testset "nut00a" begin
    dpsi, deps = nut00a(2400000.5, 53736.0)
    @test dpsi ≈ -0.9630909107115518431e-5 atol=1e-13
    @test deps ≈ 0.4063239174001678710e-4 atol=1e-13
end

@testset "nut00b" begin
    dpsi, deps = nut00b(2400000.5, 53736.0)
    @test dpsi ≈ -0.9632552291148362783e-5 atol=1e-13
    @test deps ≈ 0.4063197106621159367e-4 atol=1e-13
end

@testset "nut06a" begin
    dpsi, deps = nut06a(2400000.5, 53736.0)
    @test dpsi ≈ -0.9630912025820308797e-5 atol=1e-13
    @test deps ≈ 0.4063238496887249798e-4 atol=1e-13
end

@testset "nut80" begin
    dpsi, deps = nut80(2400000.5, 53736.0)
    @test dpsi ≈ -0.9643658353226563966e-5 atol=1e-13
    @test deps ≈ 0.4060051006879713322e-4 atol=1e-13
end

@testset "nutm80" begin
    rmatn = nutm80(2400000.5, 53736.)
    @test rmatn[1,1] ≈ 0.9999999999534999268 atol=1e-12
    @test rmatn[1,2] ≈ 0.8847935789636432161e-5 atol=1e-12
    @test rmatn[1,3] ≈ 0.3835906502164019142e-5 atol=1e-12
    @test rmatn[2,1] ≈ -0.8847780042583435924e-5 atol=1e-12
    @test rmatn[2,2] ≈ 0.9999999991366569963 atol=1e-12
    @test rmatn[2,3] ≈ -0.4060052702727130809e-4 atol=1e-12
    @test rmatn[3,1] ≈ -0.3836265729708478796e-5 atol=1e-12
    @test rmatn[3,2] ≈ 0.4060049308612638555e-4 atol=1e-12
    @test rmatn[3,3] ≈ 0.9999999991684415129 atol=1e-12
end

