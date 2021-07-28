module ERFA

using LinearAlgebra: I, cross, dot, norm, normalize

using ERFA_jll

export ERFAException
export DPI, D2PI, DR2D, DD2R, DR2AS, DAS2R, DS2R, TURNAS, DMAS2R, DTY, DAYSEC, DJY, DJC
export DJ00, DJM0, DJM00, DJM77, TTMTAI, DAU, CMPS, AULT, DC, ELG, ELB, TDB0, SRS
export WGS84, GRS80, WGS72

include("exports.jl")
include("erfa_common.jl")

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
include("m.jl")
include("n.jl")
include("o.jl")
include("p.jl")
include("r.jl")
include("s.jl")
include("t.jl")
include("u.jl")
include("x.jl")
include("z.jl")

end # module
