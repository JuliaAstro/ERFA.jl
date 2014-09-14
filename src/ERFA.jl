module ERFA

export
    PosVel

export
    eraA2af,
    eraA2tf,
    eraBp00,
    eraBp06,
    eraCal2jd,
    eraC2i00a,
    eraC2i00b,
    eraC2i06a,
    eraDat,
    eraD2dtf,
    eraD2tf,
    eraDtf2d,
    eraEe00a,
    eraEe00b,
    eraEe06a,
    eraEect00,
    eraEo06a,
    eraEpb,
    eraEpb2jd,
    eraEpj,
    eraEpj2jd,
    eraEpv00,
    eraEqeq94,
    eraEra00,
    eraFad03,
    eraFae03,
    eraFaf03,
    eraFaju03,
    eraFal03,
    eraFalp03,
    eraFama03,
    eraFame03,
    eraFane03,
    eraFaom03,
    eraFapa03,
    eraFasa03,
    eraFaur03,
    eraFave03,
    eraGmst82,
    eraGst00b,
    eraGst94,
    eraJd2cal,
    eraJdcalf,
    eraLDBODY,
    eraLdn,
    eraNum00a,
    eraNum00b,
    eraNum06a,
    eraNumat,
    eraNut00a,
    eraNut00b,
    eraNut06a,
    eraNut80,
    eraNutm80,
    eraObl06,
    eraObl80,
    eraPlan94,
    eraPmat00,
    eraPmat06,
    eraPmat76,
    eraPnm00a,
    eraPnm00b,
    eraPnm06a,
    eraPnm80,
    eraTaitt,
    eraTaiut1,
    eraTaiutc,
    eraTcbtdb,
    eraTcgtt,
    eraTdbtcb,
    eraTdbtt,
    eraTf2d,
    eraTttai,
    eraTttcg,
    eraTttdb,
    eraTtut1,
    eraUt1tai,
    eraUt1tt,
    eraUt1utc,
    eraUtctai,
    eraUtcut1

include("../deps/deps.jl")

immutable PosVel
    p1::Float64
    p2::Float64
    p3::Float64
    v1::Float64
    v2::Float64
    v3::Float64
end

immutable eraLDBODY
    bm::Float64
    dl::Float64
    pv::PosVel
end

function eraCal2jd(iy::Integer, imo::Integer, id::Integer)
    r1 = [0.]
    r2 = [0.]
    i = ccall((:eraCal2jd,liberfa), Cint,
              (Cint,Cint,Cint,Ptr{Float64},Ptr{Float64}),
              iy, imo, id, r1, r2)
    @assert i == 0
    r1[1], r2[1]
end

function eraJd2cal(d1::Real, d2::Real)
    iy = Int32[0]
    imo = Int32[0]
    id = Int32[0]
    fd = [0.]
    i = ccall((:eraJd2cal,liberfa), Cint,
              (Float64,Float64,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Float64}),
              d1, d2, iy, imo, id, fd)
    @assert i == 0
    iy[1], imo[1], id[1], fd[1]
end

function eraJdcalf(ndp::Integer, d1::Real, d2::Real)
    iymdf = Int32[0, 0, 0, 0]
    i = ccall((:eraJdcalf,liberfa), Cint,
              (Cint,Float64,Float64,Ptr{Cint}),
              ndp, d1, d2, iymdf)
    @assert i == 0
    iymdf[1], iymdf[2], iymdf[3], iymdf[4]
end

function eraDat(iy::Integer, im::Integer, id::Integer, fd::Real)
    d = [0.]
    i = ccall((:eraDat, liberfa), Cint,
              (Cint, Cint, Cint, Float64, Ptr{Float64}),
              iy, im, id, fd, d)
    @assert i == 0
    d[1]
end

function eraD2dtf(scale::ByteString, ndp::Integer, d1::Real, d2::Real)
    iy = Int32[0]
    imo = Int32[0]
    id = Int32[0]
    ihmsf = Int32[0, 0, 0, 0]
    i = ccall((:eraD2dtf,liberfa), Cint,
              (Ptr{Cchar},Cint,Float64,Float64,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),
              scale, ndp, d1, d2, iy, imo, id, ihmsf)
    @assert i == 0
    iy[1], imo[1], id[1], ihmsf[1], ihmsf[2], ihmsf[3], ihmsf[4]
end

function eraDtf2d(scale::ByteString, iy::Integer, imo::Integer, id::Integer, ih::Integer, imi::Integer, sec::Real)
    r1 = [0.]
    r2 = [0.]
    i = ccall((:eraDtf2d,liberfa), Cint,
              (Ptr{Cchar},Cint,Cint,Cint,Cint,Cint,Float64,Ptr{Float64},Ptr{Float64}),
              scale, iy, imo, id, ih, imi, sec, r1, r2)
    @assert i == 0
    r1[1], r2[1]
end

function eraEpv00(date1::Float64, date2::Float64)
    pvh = zeros(6)
    pvb = zeros(6)
    i = ccall((:eraEpv00, liberfa),
              Cint,
              (Float64, Float64, Ptr{Float64}, Ptr{Float64}),
              date1, date2, pvh, pvb)
    if i == 1
        warn("date outside the range 1900-2100 AD")
    end
    pvh, pvb
end

function eraLdn(l::Array{eraLDBODY}, ob::Array{Float64}, sc::Array{Float64})
    sn = zeros(3)
    n = length(l)
    ccall((:eraLdn, liberfa),
          Void,
          (Cint, Ptr{eraLDBODY}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          n, l, ob, sc, sn)
    sn
end

function eraNumat(epsa::Real, dpsi::Real, deps::Real)
    rmatn = zeros(9)
    ccall((:eraNumat,liberfa),
          Void,
          (Float64, Float64, Float64, Ptr{Float64}),
          epsa, dpsi, deps, rmatn)
    rmatn
end

function eraPlan94(date1::Float64, date2::Float64, np::Int64)
    pv = zeros(6)
    i = ccall((:eraPlan94, liberfa),
              Cint,
              (Float64, Float64, Int64, Ptr{Float64}),
              date1, date2, np, pv)
    if i == -1
        error("illegal np,  not in range(1,8) for planet")
    elseif i == 1
        warn("year outside range(1000:3000)")
        return pv
    elseif i == 2
        error("computation failed to converge")
    elseif i == 0
        return pv
    end
end

for f in (:eraA2af,
          :eraA2tf,
          :eraD2tf)
    @eval begin
        function ($f)(ndp::Int64, a::Float64)
            s = "+"
            i = Int32[0, 0, 0, 0]
            ccall(($(Expr(:quote,f)),liberfa),
                  Void,
                  (Int64, Float64, Ptr{ASCIIString}, Ptr{Cint}),
                  ndp, a, &s, i)
            s[1], i[1], i[2], i[3], i[4]
        end
    end
end
    
for f in (:eraBp00,
          :eraBp06)
    @eval begin
        function ($f)(a::Float64, b::Float64)
            rb = zeros(9)
            rp = zeros(9)
            rbp = zeros(9)
            ccall(($(Expr(:quote,f)),liberfa),
                  Void,
                  (Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
                  a, b, rb, rp, rbp)
            rb, rp, rbp
        end
    end
end

for f in (:eraC2i00a,
          :eraC2i00b,
          :eraC2i06a,
          :eraNum00a,
          :eraNum00b,
          :eraNum06a,
          :eraNutm80,
          :eraPmat00,
          :eraPmat06,
          :eraPmat76,
          :eraPnm00a,
          :eraPnm00b,
          :eraPnm06a,
          :eraPnm80)
    @eval begin
        function ($f)(a::Float64, b::Float64)
            r = zeros(9)
            ccall(($(Expr(:quote,f)),liberfa),
                  Void,
                  (Float64, Float64, Ptr{Float64}),
                  a, b, r)
            r
        end
    end
end

for f in (:eraNut00a,
          :eraNut00b,
          :eraNut06a,
          :eraNut80)
    @eval begin
        function ($f)(a::Float64, b::Float64)
            r1 = [0.]
            r2 = [0.]
            ccall(($(Expr(:quote,f)),liberfa),
                  Void,
                  (Float64, Float64, Ptr{Float64}, Ptr{Float64}),
                  a, b, r1, r2)
            r1[1], r2[1]
        end
    end
end

function eraTf2d(s::Char, ihour::Integer, imin::Integer, sec::Real)
    days = [0.]
    i = ccall((:eraTf2d,liberfa), Cint,
              (Cchar,Cint,Cint,Float64,Ptr{Float64}),
              s, ihour, imin, sec, days)
    @assert i == 0
    days[1]
end

for f in (:eraFad03,
          :eraFae03,
          :eraFaf03,
          :eraFaju03,
          :eraFal03,
          :eraFalp03,
          :eraFama03,
          :eraFame03,
          :eraFane03,
          :eraFaom03,
          :eraFapa03,
          :eraFasa03,
          :eraFaur03,
          :eraFave03)
    @eval ($f)(d::Real) = ccall(($(Expr(:quote,f)),liberfa), Float64, (Float64,), d)
end

for f in (:eraEe00a,
          :eraEe00b,
          :eraEe06a,
          :eraEect00,
          :eraEo06a,
          :eraEpb,
          :eraEpj,
          :eraEqeq94,
          :eraEra00,
          :eraGmst82,
          :eraGst00b,
          :eraGst94,
          :eraObl06,
          :eraObl80)
    @eval ($f)(d1::Real, d2::Real) = ccall(($(Expr(:quote,f)),liberfa), Float64, (Float64,Float64), d1, d2)
end

for f in (:eraTaitt,
          :eraTaiutc,
          :eraTcbtdb,
          :eraTcgtt,
          :eraTdbtcb,
          :eraTttai,
          :eraTttcg,
          :eraUtctai)
    @eval begin
        function ($f)(a::Float64, b::Float64)
            r1 = [0.]
            r2 = [0.]
            i = ccall(($(Expr(:quote,f)),liberfa), Cint,
                      (Float64, Float64, Ptr{Float64}, Ptr{Float64}),
                      a, b, r1, r2)
            @assert i == 0
            r1[1], r2[1]
        end
    end
end

for f in (:eraEpb2jd,
          :eraEpj2jd)
    @eval begin
        function ($f)(d::Real)
            r1 = [0.]
            r2 = [0.]
            ccall(($(Expr(:quote,f)),liberfa), Void,
                  (Float64, Ptr{Float64}, Ptr{Float64}),
                  d, r1, r2)
            r1[1], r2[1]
        end
    end
end

for f in (:eraTaiut1,
          :eraTdbtt,
          :eraTttdb,
          :eraTtut1,
          :eraUt1tai,
          :eraUt1tt,
          :eraUt1utc,
          :eraUtcut1)
    @eval begin
        function ($f)(a::Float64, b::Float64, c::Float64)
            r1 = [0.]
            r2 = [0.]
            i = ccall(($(Expr(:quote,f)),liberfa), Cint,
                      (Float64, Float64, Float64, Ptr{Float64}, Ptr{Float64}),
                      a, b, c, r1, r2)
            @assert i == 0
            r1[1], r2[1]
        end
    end
end

end # module
