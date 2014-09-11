module ERFA

export
    eraCal2jd,
    eraDat,
    eraD2dtf,
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
    eraPmat00,
    eraPmat06,
    eraPmat76,
    eraNum00a,
    eraNum00b,
    eraNum06a,
    eraNut00a,
    eraNut00b,
    eraNut06a,
    eraNut80,
    eraNutm80,
    eraObl80,
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

for f in (:eraNutm80,
          :eraPmat00,
          :eraPmat06,
          :eraPmat76)
    @eval begin
        function ($f)(a::Float64, b::Float64)
            r = [zeros(3), zeros(3), zeros(3)]
            ccall(($(Expr(:quote,f)),liberfa),
                  Void,
                  (Float64, Float64, Ptr{Float64}),
                  a, b, r)
            r
        end
    end
end

for f in (:eraNum00a,
          :eraNum00b,
          :eraNum06a,
          :eraNut00a,
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
          :eraGst94)
    @eval ($f)(d1::Real, d2::Real) = ccall(($(Expr(:quote,f)),liberfa), Float64, (Float64,Float64), d1, d2)
end

for f in (:eraTaitt,
          :eraTaiutc,
          :eraTcbtdb,
          :eraTcgtt,
          :eraTdbtcb,
          :eraTttai,
          :eraTttcg,
          :eraUt1tai,
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
