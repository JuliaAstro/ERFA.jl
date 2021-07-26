using ERFA

using StaticArrays: @SVector
using Test

@testset "ERFA" begin
    @testset "Calendar Tools" begin
        u1, u2 = dtf2d("UTC", 2010, 7, 24, 11, 18, 7.318)
        a1, a2 = utctai(u1, u2)
        t1, t2 = taitt(a1, a2)
        @test d2dtf("tt", 3, t1, t2) == (2010, 7, 24, 11, 19, 13, 502)

        iy = 2008; imo = 2; id = 29
        ihour = 23; imin = 59; sec = 59.9
        d1, d2 = cal2jd(iy, imo, id)
        d = tf2d('+', ihour, imin, sec)
        d2 += d
        @test d1 == 2400000.5
        @test isapprox(d2, 54525.999999, atol = 5e-7)
        iy, imo, id, fd = jd2cal(d1, d2)
        @test (iy, imo, id) == (2008, 2, 29)
        @test isapprox(fd, 0.999999, atol = 5e-7)
        @test jdcalf(3, d1, d2) == (2008, 3, 1, 0)

        d = 2457073.05631
        e = epb(0., d)
        @test isapprox(e, 2015.1365941021, atol = 5e-11)
        d1, d2 = epb2jd(e)
        d = d1 + d2
        @test isapprox(d, 2457073.056310000, atol = 5e-10)
        e = epj(0., d)
        @test isapprox(e, 2015.1349933196, atol = 5e-11)
        d1, d2 = epj2jd(e)
        d = d1 + d2
        @test isapprox(d, 2457073.056310000, atol = 5e-10)
    end
    @testset "StaticArrays" begin
        base = gc2gde(6378.137 * 1000, 1 / 298.257223563, [6378.137 * 1000 + 100, 0, 0])
        static = gc2gde(6378.137 * 1000, 1 / 298.257223563, @SVector[6378.137 * 1000 + 100, 0, 0])
        @test base == static
    end

    include("a.jl")
    include("b.jl")
    include("c.jl")
    include("d.jl")
    include("e.jl")
    include("f.jl")
    include("g.jl")
    include("h.jl")
    include("i.jl")
    include("j.jl")
    include("l.jl")
    include("n.jl")
    include("o.jl")
    include("p.jl")
    include("r.jl")
    include("s.jl")
    include("t.jl")
    include("u.jl")
    include("x.jl")
    include("z.jl")
end
