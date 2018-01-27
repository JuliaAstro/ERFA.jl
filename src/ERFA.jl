__precompile__()

module ERFA

import Base.getindex

const depsfile = joinpath(dirname(dirname(@__FILE__)),"deps","deps.jl")
if isfile(depsfile)
    include(depsfile)
else
    error("ERFA is not properly installed. Please run Pkg.build(\"ERFA\")")
end

include("erfa_common.jl")

function ASTROM(pmt, eb::AbstractArray, eh::AbstractArray, em, v::AbstractArray, bm1, bpn::AbstractArray, along, phi, xpl, ypl, sphi, cphi, diurab, eral, refa, refb)
    ASTROM(pmt,
              (eb[1], eb[2], eb[3]),
              (eh[1], eh[2], eh[3]),
              em,
              (v[1], v[2], v[3]),
              bm1,
              (bpn[1],bpn[2],bpn[3],bpn[4],bpn[5],bpn[6],bpn[7],bpn[8],bpn[9]),
              along,
              phi,
              xpl,
              ypl,
              sphi,
              cphi,
              diurab,
              eral,
              refa,
              refb)
end

function ab(pnat,v,s,bm1)
    ppr = zeros(3)
    ccall((:eraAb,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Cdouble,Cdouble,Ptr{Cdouble}),
          pnat,v,s,bm1,ppr)
    ppr
end

function apcg(date1,date2,ebpv,ehp)
    astrom = ASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApcg,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ref{ASTROM}),
          date1,date2,ebpv,ehp,astrom)
    astrom
end

function apcg13(date1,date2)
    astrom = ASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApcg13,liberfa),Void,
          (Cdouble,Cdouble,Ref{ASTROM}),
          date1,date2,astrom)
    astrom
end

function apci(date1,date2,ebpv,ehp,x,y,s)
    astrom = ASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApci,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Cdouble,Cdouble,Cdouble,Ref{ASTROM}),
          date1,date2,ebpv,ehp,x,y,s,astrom)
    astrom
end

function apci13(date1,date2)
    astrom = ASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    eo = [0.]
    ccall((:eraApci13,liberfa),Void,
          (Cdouble,Cdouble,Ref{ASTROM},Ptr{Cdouble}),
          date1,date2,astrom,eo)
    astrom, eo[1]
end

function apco(date1,date2,ebpv,ehp,x,y,s,theta,elong,phi,hm,xp,yp,sp,refa,refb)
    astrom = ASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApco,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ref{ASTROM}),
          date1,date2,ebpv,ehp,x,y,s,theta,elong,phi,hm,xp,yp,sp,refa,refb,astrom)
    astrom
end

function apco13(utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl)
    astrom = ASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    eo = [0.]
    i = ccall((:eraApco13,liberfa),Cint,
              (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ref{ASTROM},Ptr{Cdouble}),
              utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl,astrom,eo)
    if i == -1
        error("unacceptable date")
    elseif i == +1
        warn("dubious year")
    end
    astrom, eo[1]
end

function apcs(date1,date2,pv,ebpv,ehp)
    astrom = ASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApcs,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ref{ASTROM}),
          date1,date2,pv,ebpv,ehp,astrom)
    astrom
end

function apcs13(date1,date2,pv)
    astrom = ASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApcs13,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ref{ASTROM}),
          date1,date2,pv,astrom)
    astrom
end

function aper(theta,astrom)
    ccall((:eraAper,liberfa),Void,
          (Cdouble,Ref{ASTROM}),
          theta,astrom)
    astrom
end

function aper13(ut11,ut12,astrom)
    ccall((:eraAper13,liberfa),Void,
          (Cdouble,Cdouble,Ref{ASTROM}),
          ut11,ut12,astrom)
    astrom
end

function apio(sp,theta,elong,phi,hm,xp,yp,refa,refb)
    astrom = ASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApio,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ref{ASTROM}),
          sp,theta,elong,phi,hm,xp,yp,refa,refb,astrom)
    astrom
end

function apio13(utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl)
    astrom = ASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    i = ccall((:eraApio13,liberfa),Cint,
              (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ref{ASTROM}),
              utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl,astrom)
    if i == -1
        error("unacceptable date")
    elseif i == +1
        warn("dubious year")
    end
    astrom
end

function atci13(rc,dc,pr,pd,px,rv,date1,date2)
    ri = [0.]
    di = [0.]
    eo = [0.]
    ccall((:eraAtci13,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          rc,dc,pr,pd,px,rv,date1,date2,ri,di,eo)
    ri[1], di[1], eo[1]
end

function atciq(rc,dc,pr,pd,px,rv,astrom)
    ri = [0.]
    di = [0.]
    ccall((:eraAtciq,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ref{ASTROM},Ptr{Cdouble},Ptr{Cdouble}),
          rc,dc,pr,pd,px,rv,astrom,ri,di)
    ri[1], di[1]
end

function atciqn(rc,dc,pr,pd,px,rv,astrom,b::Array{LDBODY})
    ri = [0.]
    di = [0.]
    n = length(b)
    ccall((:eraAtciqn,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ref{ASTROM},Cint,Ptr{LDBODY},Ptr{Cdouble},Ptr{Cdouble}),
          rc,dc,pr,pd,px,rv,astrom,n,b,ri,di)
    ri[1], di[1]
end

function atciqz(rc,dc,astrom)
    ri = [0.]
    di = [0.]
    ccall((:eraAtciqz,liberfa),Void,
          (Cdouble,Cdouble,Ref{ASTROM},Ptr{Cdouble},Ptr{Cdouble}),
          rc,dc,astrom,ri,di)
    ri[1], di[1]
end

function atco13(rc,dc,pr,pd,px,rv,utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl)
    aob = [0.]
    zob = [0.]
    hob = [0.]
    dob = [0.]
    rob = [0.]
    eo = [0.]
    i = ccall((:eraAtco13,liberfa),Cint,
              (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
              rc,dc,pr,pd,px,rv,utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl,aob,zob,hob,dob,rob,eo)
    if i == -1
        error("unacceptable date")
    elseif i == +1
        warn("dubious year")
    end
    aob[1],zob[1],hob[1],dob[1],rob[1],eo[1]
end

function atic13(ri,di,date1,date2)
    rc = [0.]
    dc = [0.]
    eo = [0.]
    ccall((:eraAtic13,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          ri,di,date1,date2,rc,dc,eo)
    rc[1],dc[1],eo[1]
end

function aticq(ri,di,astrom)
    rc = [0.]
    dc = [0.]
    ccall((:eraAticq,liberfa),Void,
          (Cdouble,Cdouble,Ref{ASTROM},Ptr{Cdouble},Ptr{Cdouble}),
          ri,di,astrom,rc,dc)
    rc[1],dc[1]
end

function aticqn(ri,di,astrom,b::Array{LDBODY})
    rc = [0.]
    dc = [0.]
    n = length(b)
    ccall((:eraAticqn,liberfa),Void,
          (Cdouble,Cdouble,Ref{ASTROM},Cint,Ptr{LDBODY},Ptr{Cdouble},Ptr{Cdouble}),
          ri,di,astrom,n,b,rc,dc)
    rc[1],dc[1]
end

function atio13(ri,di,utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl)
    aob = [0.]
    zob = [0.]
    hob = [0.]
    dob = [0.]
    rob = [0.]
    i = ccall((:eraAtio13,liberfa),Cint,
              (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
              ri,di,utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl,aob,zob,hob,dob,rob)
    if i == -1
        error("unacceptable date")
    elseif i == +1
        warn("dubious year")
    end
    aob[1],zob[1],hob[1],dob[1],rob[1]
end

function atioq(ri,di,astrom)
    aob = [0.]
    zob = [0.]
    hob = [0.]
    dob = [0.]
    rob = [0.]
    ccall((:eraAtioq,liberfa),Void,
          (Cdouble,Cdouble,Ref{ASTROM},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          ri,di,astrom,aob,zob,hob,dob,rob)
    aob[1],zob[1],hob[1],dob[1],rob[1]
end

function atoc13(typeofcoordinates,ob1,ob2,utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl)
    rc = [0.]
    dc = [0.]
    if !(typeofcoordinates in ('R', 'r', 'H', 'h', 'A', 'a'))
        typeofcoordinates = 'A'
    end
    i = ccall((:eraAtoc13,liberfa),Cint,
              (Ref{Char},Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
              typeofcoordinates,ob1,ob2,utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl,rc,dc)
    if i == -1
        error("unacceptable date")
    elseif i == +1
        warn("dubious year")
    end
    rc[1],dc[1]
end

function atoi13(typeofcoordinates,ob1,ob2,utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl)
    ri = [0.]
    di = [0.]
    i = ccall((:eraAtoi13,liberfa),Cint,
              (Ref{Char},Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
              typeofcoordinates,ob1,ob2,utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl,ri,di)
    if i == -1
        error("unacceptable date")
    elseif i == +1
        warn("dubious year")
    end
    ri[1],di[1]
end

function atoiq(typeofcoordinates,ob1,ob2,astrom)
    ri = [0.]
    di = [0.]
    ccall((:eraAtoiq,liberfa),
          Void,(Ref{Char},Cdouble,Cdouble,Ref{ASTROM},Ptr{Cdouble},Ptr{Cdouble}),
          typeofcoordinates,ob1,ob2,astrom,ri,di)
    ri[1],di[1]
end

function bi00()
    dpsibi = [0.]
    depsbi = [0.]
    dra = [0.]
    ccall((:eraBi00,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          dpsibi,depsbi,dra)
    dpsibi[1], depsbi[1], dra[1]
end

function bpn2xy(rbpn)
    x = [0.]
    y = [0.]
    ccall((:eraBpn2xy,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          rbpn,x,y)
    x[1], y[1]
end

function c2ibpn(date1,date2,rbpn)
    rc2i = zeros((3,3))
    ccall((:eraC2ibpn,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          date1,date2,rbpn,rc2i)
    rc2i
end

function c2s(p)
    theta = [0.]
    phi = [0.]
    ccall((:eraC2s,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          p,theta,phi)
    theta[1], phi[1]
end

function cal2jd(iy, imo, id)
    r1 = [0.]
    r2 = [0.]
    i = ccall((:eraCal2jd,liberfa), Cint,
              (Cint,Cint,Cint,Ptr{Cdouble},Ptr{Cdouble}),
              iy, imo, id, r1, r2)
    @assert i == 0
    r1[1], r2[1]
end

function cp(p)
    c = zeros(3)
    ccall((:eraCp,liberfa),Void,(Ptr{Cdouble},Ptr{Cdouble}),
          p,c)
    c
end

function cpv(pv)
    c = zeros((2,3))
    ccall((:eraCpv,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          pv,c)
    c
end

function cr(p)
    r = zeros((3,3))
    ccall((:eraCr,liberfa),Void,(Ptr{Cdouble},Ptr{Cdouble}),
          p,r)
    r
end

function dtdb(date1,date2,ut,elong,u,v)
    ccall((:eraDtdb,liberfa),Cdouble,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble),
          date1,date2,ut,elong,u,v)
end

function eform(n)
    a = [0.]
    f = [0.]
    i = ccall((:eraEform,liberfa),Cint,
              (Cint,Ptr{Cdouble},Ptr{Cdouble}),
              n,a,f)
    if i == -1
        error("illegal identifier")
    end
    a[1], f[1]
end

function eors(rnpb,s)
    ccall((:eraEors,liberfa),Cdouble,
          (Ptr{Cdouble},Cdouble),
          rnpb,s)
end

function jd2cal(d1, d2)
    iy = Int32[0]
    imo = Int32[0]
    id = Int32[0]
    fd = [0.]
    i = ccall((:eraJd2cal,liberfa), Cint,
              (Cdouble,Cdouble,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cdouble}),
              d1, d2, iy, imo, id, fd)
    @assert i == 0
    iy[1], imo[1], id[1], fd[1]
end

function jdcalf(ndp, d1, d2)
    iymdf = Int32[0, 0, 0, 0]
    i = ccall((:eraJdcalf,liberfa), Cint,
              (Cint,Cdouble,Cdouble,Ptr{Cint}),
              ndp, d1, d2, iymdf)
    @assert i == 0
    iymdf[1], iymdf[2], iymdf[3], iymdf[4]
end

function dat(iy, im, id, fd)
    d = [0.]
    i = ccall((:eraDat, liberfa), Cint,
              (Cint, Cint, Cint, Cdouble, Ptr{Cdouble}),
              iy, im, id, fd, d)
    @assert i == 0
    d[1]
end

function d2dtf(scale::AbstractString, ndp, d1, d2)
    iy = Int32[0]
    imo = Int32[0]
    id = Int32[0]
    ihmsf = Int32[0, 0, 0, 0]
    i = ccall((:eraD2dtf,liberfa), Cint,
              (Ptr{Cchar},Cint,Cdouble,Cdouble,Ptr{Cint},Ptr{Cint},Ptr{Cint},Ptr{Cint}),
              scale, ndp, d1, d2, iy, imo, id, ihmsf)
    @assert i == 0
    iy[1], imo[1], id[1], ihmsf[1], ihmsf[2], ihmsf[3], ihmsf[4]
end

function dtf2d(scale::AbstractString, iy, imo, id, ih, imi, sec)
    r1 = [0.]
    r2 = [0.]
    i = ccall((:eraDtf2d,liberfa), Cint,
              (Ptr{Cchar},Cint,Cint,Cint,Cint,Cint,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
              scale, iy, imo, id, ih, imi, sec, r1, r2)
    @assert i == 0
    r1[1], r2[1]
end

function epv00(date1, date2)
    pvh = zeros((2,3))
    pvb = zeros((2,3))
    i = ccall((:eraEpv00, liberfa),
              Cint,
              (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),
              date1, date2, pvh, pvb)
    if i == 1
        warn("date outside the range 1900-2100 AD")
    end
    pvh, pvb
end

function fk5hip()
    r5h = zeros((3,3))
    s5h = zeros(3)
    ccall((:eraFk5hip,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          r5h,s5h)
    r5h,s5h
end

function fk5hz(r5,d5,date1,date2)
    rh = [0.]
    dh = [0.]
    ccall((:eraFk5hz,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          r5,d5,date1,date2,rh,dh)
    rh[1],dh[1]
end

function fw2xy(gamb,phib,psi,eps)
    x =[0.]
    y = [0.]
    ccall((:eraFw2xy,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          gamb,phib,psi,eps,x,y)
    x[1], y[1]
end

function gc2gd(n,xyz)
    elong = [0.]
    phi = [0.]
    height = [0.]
    i = ccall((:eraGc2gd,liberfa),Cint,
              (Cint,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
              n,xyz,elong,phi,height)
    if i == -1
        error("illegal identifier")
    elseif i == -2
        error("internal error")
    end
    elong[1],phi[1],height[1]
end

function gc2gde(a,f,xyz)
    elong = [0.]
    phi = [0.]
    height = [0.]
    i = ccall((:eraGc2gde,liberfa),Cint,
              (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
              a,f,xyz,elong,phi,height)
    if i == -1
        error("illegal f")
    elseif i == -2
        error("internal a")
    end
    elong[1],phi[1],height[1]
end

function gd2gc(n,elong,phi,height)
    xyz = zeros(3)
    i = ccall((:eraGd2gc,liberfa),Cint,
              (Cint,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
              n,elong,phi,height,xyz)
    if i == -1
        error("illegal identifier")
    elseif i == -2
        error("illegal case")
    end
    xyz
end

function gd2gce(a,f,elong,phi,height)
    xyz = zeros(3)
    i = ccall((:eraGd2gce,liberfa),Cint,
              (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
              a,f,elong,phi,height,xyz)
    if i == -1
        error("illegal case")
    end
    xyz
end

function gst06(uta,utb,tta,ttb,rnpb)
    ccall((:eraGst06,liberfa),Cdouble,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
          uta,utb,tta,ttb,rnpb)
end

function hfk5z(rh,dh,date1,date2)
    r5 = [0.]
    d5 = [0.]
    dr5 = [0.]
    dd5 = [0.]
    ccall((:eraHfk5z,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          rh,dh,date1,date2,r5,d5,dr5,dd5)
    r5[1],d5[1],dr5[1],dd5[1]
end

function LDBODY(bm, dl, pv::AbstractArray)
    LDBODY(bm, dl, (pv...))
end

function ld(bm,p,q,e,em,dlim)
    p1 = zeros(3)
    ccall((:eraLd,liberfa),Void,
          (Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Cdouble,Cdouble,Ptr{Cdouble}),
          bm,p,q,e,em,dlim,p1)
    p1
end

function ldn(l::Array{LDBODY}, ob::AbstractArray, sc::AbstractArray)
    sn = zeros(3)
    n = length(l)
    ccall((:eraLdn, liberfa),Void,
          (Cint, Ptr{LDBODY}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
          n, l, ob, sc, sn)
    sn
end

function ldsun(p,e,em)
    p1 = zeros(3)
    ccall((:eraLdsun,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Cdouble,Ptr{Cdouble}),
          p,e,em,p1)
    p1
end

function pmpx(rc,dc,pr,pd,px,rv,pmt,vob)
    pco = zeros(3)
    ccall((:eraPmpx,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          rc,dc,pr,pd,px,rv,pmt,vob,pco)
    pco
end

function numat(epsa, dpsi, deps)
    rmatn = zeros((3,3))
    ccall((:eraNumat,liberfa),Void,
          (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
          epsa, dpsi, deps, rmatn)
    rmatn
end

function p06e(date1,date2)
    eps0 = [0.]
    psia = [0.]
    oma = [0.]
    bpa = [0.]
    bqa = [0.]
    pia = [0.]
    bpia = [0.]
    epsa = [0.]
    chia = [0.]
    za = [0.]
    zetaa = [0.]
    thetaa = [0.]
    pa = [0.]
    gam = [0.]
    phi = [0.]
    psi = [0.]
    ccall((:eraP06e,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          date1,date2,eps0,psia,oma,bpa,bqa,pia,bpia,epsa,chia,za,zetaa,thetaa,pa,gam,phi,psi)
    eps0[1],psia[1],oma[1],bpa[1],bqa[1],pia[1],bpia[1],epsa[1],chia[1],za[1],zetaa[1],thetaa[1],pa[1],gam[1],phi[1],psi[1]
end

function p2s(p)
    theta = [0.]
    phi = [0.]
    r = [0.]
    ccall((:eraP2s,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          p,theta,phi,r)
    theta[1], phi[1], r[1]
end

function p2pv(p)
    pv = zeros((2,3))
    ccall((:eraP2pv,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          p,pv)
    pv
end

function pb06(date1,date2)
    bzeta = [0.]
    bz = [0.]
    btheta = [0.]
    ccall((:eraPb06,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          date1,date2,bzeta,bz,btheta)
    bzeta[1],bz[1],btheta[1]
end


function pfw06(date1,date2)
    gamb = [0.]
    phib = [0.]
    psib = [0.]
    epsa = [0.]
    ccall((:eraPfw06,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          date1,date2,gamb,phib,psib,epsa)
    gamb[1],phib[1],psib[1],epsa[1]
end

function plan94(date1, date2, np)
    pv = zeros((2,3))
    i = ccall((:eraPlan94, liberfa),Cint,
              (Cdouble, Cdouble, Cint, Ptr{Cdouble}),
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

function pm(p)
    ccall((:eraPm,liberfa),Cdouble,(Ptr{Cdouble},),p)
end

function pmsafe(ra1,dec1,pmr1,pmd1,px1,rv1,ep1a,ep1b,ep2a,ep2b)
    ra2 = [0.]
    dec2 = [0.]
    pmr2 = [0.]
    pmd2 = [0.]
    px2 = [0.]
    rv2 = [0.]
    i = ccall((:eraPmsafe,liberfa),Cint,
              (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
              ra1,dec1,pmr1,pmd1,px1,rv1,ep1a,ep1b,ep2a,ep2b,ra2,dec2,pmr2,pmd2,px2,rv2)
    if i == -1
        error("system error")
    elseif i == 1
        warn("distance overridden")
    elseif i == 2
        warn("excessive velocity")
    elseif i == 4
        error("solution didn't converge")
    end
    ra2[1],dec2[1],pmr2[1],pmd2[1],px2[1],rv2[1]
end

function pn(p)
    r = [0.]
    u = zeros(3)
    ccall((:eraPn,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          p,r,u)
    r[1],u
end

function ppsp(a,s,b)
    apsb = zeros(3)
    ccall((:eraPpsp,liberfa),Void,
          (Ptr{Cdouble},Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          a,s,b,apsb)
    apsb
end

function prec76(ep01,ep02,ep11,ep12)
    zeta = [0.]
    z = [0.]
    theta = [0.]
    ccall((:eraPrec76,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          ep01,ep02,ep11,ep12,zeta,z,theta)
    zeta[1],z[1],theta[1]
end

function pv2s(pv)
    theta = [0.]
    phi = [0.]
    r = [0.]
    td = [0.]
    pd = [0.]
    rd = [0.]
   ccall((:eraPv2s,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          pv,theta,phi,r,td,pd,rd)
    theta[1], phi[1], r[1], td[1], pd[1], rd[1]
end

function pv2p(pv)
    p = zeros(3)
    ccall((:eraPv2p,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          pv,p)
    p
end

function pvdpv(a,b)
    adb = zeros(2)
    ccall((:eraPvdpv,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          a,b,adb)
    adb
end

function pvm(pv)
    s = [0.]
    r = [0.]
    ccall((:eraPvm,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          pv,r,s)
    r[1], s[1]
end

function pvstar(pv)
    ra = [0.]
    dec = [0.]
    pmr = [0.]
    pmd = [0.]
    px = [0.]
    rv = [0.]
    i = ccall((:eraPvstar,liberfa),Cint,
              (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
              pv,ra,dec,pmr,pmd,px,rv)
    if i == -1
        error("superluminal speed")
    elseif i == -2
        warn("null position vector")
        return ra[1],dec[1],pmr[1],pmd[1],px[1],rv[1]
    end
    ra[1],dec[1],pmr[1],pmd[1],px[1],rv[1]
end

function pvtob(elong,phi,height,xp,yp,sp,theta)
    pv = zeros((2,3))
    ccall((:eraPvtob,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
          elong,phi,height,xp,yp,sp,theta,pv)
    pv
end

function refco(phpa,tk,rh,wl)
    refa = [0.]
    refb = [0.]
    ccall((:eraRefco,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          phpa,tk,rh,wl,refa,refb)
    refa[1],refb[1]
end

function pvu(dt,pv)
    upv = zeros((2,3))
    ccall((:eraPvu,liberfa),Void,
          (Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          dt,pv,upv)
    upv
end

function pvup(dt,pv)
    p = zeros(3)
    ccall((:eraPvup,liberfa),Void,
          (Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          dt,pv,p)
    p
end

function rm2v(r)
    w = zeros(3)
    ccall((:eraRm2v,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          r,w)
    w
end

function rv2m(w)
    r = zeros((3,3))
    ccall((:eraRv2m,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          w,r)
    r
end

function rxr(a,b)
    atb = zeros((3,3))
    ccall((:eraRxr,liberfa),
          Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          a,b,atb)
    atb
end

function s2c(theta,phi)
    c = zeros(3)
    ccall((:eraS2c,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble}),
          theta,phi,c)
    c
end

function s2p(theta,phi,r)
    p = zeros(3)
    ccall((:eraS2p,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
          theta,phi,r,p)
    p
end

function s2pv(theta,phi,r,td,pd,rd)
    pv = zeros((2,3))
    ccall((:eraS2pv,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
          theta,phi,r,td,pd,rd,pv)
    pv
end

function s2xpv(s1,s2,pv)
    spv = zeros((2,3))
    ccall((:eraS2xpv,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          s1,s2,pv,spv)
    spv
end

function starpm(ra1,dec1,pmr1,pmd1,px1,rv1,ep1a,ep1b,ep2a,ep2b)
    ra2 = [0.]
    dec2 = [0.]
    pmr2 = [0.]
    pmd2 = [0.]
    px2 = [0.]
    rv2 = [0.]
    i = ccall((:eraStarpm,liberfa),Cint,
              (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
              ra1,dec1,pmr1,pmd1,px1,rv1,ep1a,ep1b,ep2a,ep2b,ra2,dec2,pmr2,pmd2,px2,rv2)
    if i == -1
        error("system error")
    elseif i == 1
        warn("distance overridden")
        return ra2[1],dec2[1],pmr2[1],pmd2[1],px2[1],rv2[1]
    elseif i == 2
        warn("excessive velocity")
        return ra2[1],dec2[1],pmr2[1],pmd2[1],px2[1],rv2[1]
    elseif i == 4
        error("solution didn't converge")
    end
    ra2[1],dec2[1],pmr2[1],pmd2[1],px2[1],rv2[1]
end

function starpv(ra,dec,pmr,pmd,px,rv)
    pv = zeros((2,3))
    i = ccall((:eraStarpv,liberfa),Cint,
              (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
              ra,dec,pmr,pmd,px,rv,pv)
    if i == 1
        warn("distance overridden")
        return pv
    elseif i == 2
        warn("excessive speed ")
        return pv
    elseif i == 4
        error("solution didn't converge")
    end
    pv
end

function sxp(s,p)
    sp = zeros(3)
    ccall((:eraSxp,liberfa),Void,
          (Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          s,p,sp)
    sp
end

function sxpv(s,pv)
    spv = zeros((2,3))
    ccall((:eraSxpv,liberfa),Void,
          (Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          s,pv,spv)
    spv
end

function tr(r)
    rt = zeros((3,3))
    ccall((:eraTr,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          r,rt)
    rt
end

function xy06(date1,date2)
    x = [0.]
    y = [0.]
    ccall((:eraXy06,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          date1,date2,x,y)
    x[1], y[1]
end

for name in ("af2a",
          "tf2a",
          "tf2d")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(s,ideg,iamin,asec)
            rad = [0.]
            i = ccall(($fc,liberfa),Cint,
                       (Cchar,Cint,Cint,Cdouble,Ptr{Cdouble}),
                       s,ideg,iamin,asec,rad)
            @assert i == 0
            rad[1]
        end
    end
end

function ir()
    r = zeros((3,3))
    ccall((:eraIr,liberfa),Void,(Ptr{Cdouble},),r)
    r
end

function zp()
    p = zeros(3)
    #ccall((:eraZp,liberfa),Void,(Ptr{Cdouble},),p)
    p
end

function zpv()
    pv = zeros((2,3))
    #ccall((:eraZpv,liberfa),Void,(Ptr{Cdouble},),pv)
    pv
end

function zr()
    r = zeros((3,3))
    #ccall((:eraZr,liberfa),Void,(Ptr{Cdouble},),r)
    r
end

for name in ("a2af",
          "a2tf",
          "d2tf")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(ndp, a)
            s = Vector{UInt8}(1)
            s[1] = '+'
            i = Int32[0, 0, 0, 0]
            ccall(($fc,liberfa),Void,
                  (Cint, Cdouble, Ptr{UInt8}, Ptr{Cint}),
                  ndp, a, s, i)
            Char(s[1]), i[1], i[2], i[3], i[4]
        end
    end
end

for name in ("anp",
          "anpm")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(a)
            ccall(($fc,liberfa),Cdouble,(Cdouble,),a)
        end
    end
end

for name in ("bp00",
          "bp06")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(a, b)
            rb = zeros((3,3))
            rp = zeros((3,3))
            rbp = zeros((3,3))
            ccall(($fc,liberfa),
                  Void,
                  (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
                  a, b, rb, rp, rbp)
            rb, rp, rbp
        end
    end
end


for name in ("c2ixys",
          "pom00")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(x, y, s)
            r = zeros((3,3))
            ccall(($fc,liberfa),Void,
                  (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
                  x, y, s, r)
            r
        end
    end
end

for name in ("c2tcio",
          "c2teqx")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(rc2i, era, rpom)
            rc2t = zeros((3,3))
            ccall(($fc,liberfa),Void,
                  (Ptr{Cdouble},Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
                  rc2i,era,rpom,rc2t)
            rc2t
        end
    end
end

for name in ("c2ixy",
          "fw2m")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(x, y, s, t)
            r = zeros((3,3))
            ccall(($fc,liberfa),Void,
                  (Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
                  x, y, s, t, r)
            r
        end
    end
end

for name in ("c2t00a",
          "c2t00b",
          "c2t06a")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(tta,ttb,uta,utb,xp,yp)
            rc2t =zeros((3,3))
            ccall(($fc,liberfa),Void,
                  (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
                  tta,ttb,uta,utb,xp,yp,rc2t)
            rc2t
        end
    end
end


for name in ("c2tpe",
          "c2txy")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(tta,ttb,uta,utb,x,y,xp,yp)
            rc2t =zeros((3,3))
            ccall(($fc,liberfa),Void,
                  (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
                  tta,ttb,uta,utb,x,y,xp,yp,rc2t)
            rc2t
        end
    end
end

for name in ("c2i00a",
          "c2i00b",
          "c2i06a",
          "ecm06",
          "num00a",
          "num00b",
          "num06a",
          "nutm80",
          "pmat00",
          "pmat06",
          "pmat76",
          "pnm00a",
          "pnm00b",
          "pnm06a",
          "pnm80")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(a, b)
            r = zeros((3,3))
            ccall(($fc,liberfa),Void,
                  (Cdouble, Cdouble, Ptr{Cdouble}),
                  a, b, r)
            r
        end
    end
end

for name in ("pn00",
          "pn06")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(date1, date2, dpsi, deps)
            epsa = [0.]
            rb = zeros((3,3))
            rp = zeros((3,3))
            rbp = zeros((3,3))
            rn = zeros((3,3))
            rbpn = zeros((3,3))
            ccall(($fc, liberfa),Void,
                  (Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
                  date1, date2, dpsi, deps, epsa, rb, rp, rbp, rn, rbpn)
            epsa[1], rb, rp, rbp, rn, rbpn
        end
    end
end

for name in ("pmp",
          "ppp",
          "pxp")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(a,b)
            ab = zeros(3)
            ccall(($fc,liberfa),Void,
                  (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  a,b,ab)
            ab
        end
    end
end

for name in ("pvmpv",
          "pvppv",
          "pvxpv")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(a,b)
            ab = zeros((2,3))
            ccall(($fc,liberfa),Void,
                  (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  a,b,ab)
            ab
        end
    end
end

for name in ("pn00a",
          "pn00b",
          "pn06a")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(date1, date2)
            dpsi = [0.]
            deps = [0.]
            epsa = [0.]
            rb = zeros((3,3))
            rp = zeros((3,3))
            rbp = zeros((3,3))
            rn = zeros((3,3))
            rbpn = zeros((3,3))
            ccall(($fc, liberfa),Void,
                  (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
                  date1, date2, dpsi, deps, epsa, rb, rp, rbp, rn, rbpn)
            dpsi[1], deps[1], epsa[1], rb, rp, rbp, rn, rbpn
        end
    end
end

for name in ("nut00a",
          "nut00b",
          "nut06a",
          "nut80",
          "pr00",
          "g2icrs",
          "icrs2g")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(a, b)
            r1 = [0.]
            r2 = [0.]
            ccall(($fc,liberfa),Void,
                  (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),
                  a, b, r1, r2)
            r1[1], r2[1]
        end
    end
end

for name in ("fk52h",
          "h2fk5")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(ra,dec,dra,ddec,px,rv)
            r = [0.]
            d = [0.]
            dr = [0.]
            dd = [0.]
            p = [0.]
            v = [0.]
            ccall(($fc,liberfa),Void,
                  (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  ra,dec,dra,ddec,px,rv,r,d,dr,dd,p,v)
            r[1],d[1],dr[1],dd[1],p[1],v[1]
        end
    end
end

for name in ("fad03",
          "fae03",
          "faf03",
          "faju03",
          "fal03",
          "falp03",
          "fama03",
          "fame03",
          "fane03",
          "faom03",
          "fapa03",
          "fasa03",
          "faur03",
          "fave03")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval ($f)(d) = ccall(($fc,liberfa), Cdouble, (Cdouble,), d)
end

for name in ("ee00",
          "gmst00",
          "gmst06",
          "gst00a",
          "gst06a",
          "s00",
          "s06")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval ($f)(d1,d2,t1,t2) = ccall(($fc,liberfa), Cdouble, (Cdouble,Cdouble,Cdouble,Cdouble), d1,d2,t1,t2)
end

for name in ("ee00a",
          "ee00b",
          "ee06a",
          "eect00",
          "eo06a",
          "epb",
          "epj",
          "eqeq94",
          "era00",
          "gmst82",
          "gst00b",
          "gst94",
          "obl06",
          "obl80",
          "s00a",
          "s00b",
          "s06a",
          "sp00")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval ($f)(d1, d2) = ccall(($fc,liberfa), Cdouble, (Cdouble,Cdouble), d1, d2)
end

for name in ("taitt",
          "taiutc",
          "tcbtdb",
          "tcgtt",
          "tdbtcb",
          "tttai",
          "tttcg",
          "utctai")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(a, b)
            r1 = [0.]
            r2 = [0.]
            i = ccall(($fc,liberfa), Cint,
                      (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),
                      a, b, r1, r2)
            @assert i == 0
            r1[1], r2[1]
        end
    end
end

for name in ("epb2jd",
          "epj2jd")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(d)
            r1 = [0.]
            r2 = [0.]
            ccall(($fc,liberfa), Void,
                  (Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),
                  d, r1, r2)
            r1[1], r2[1]
        end
    end
end

for name in ("taiut1",
          "tdbtt",
          "tttdb",
          "ttut1",
          "ut1tai",
          "ut1tt",
          "ut1utc",
          "utcut1")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(a, b, c)
            r1 = [0.]
            r2 = [0.]
            i = ccall(($fc,liberfa), Cint,
                      (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),
                      a, b, c, r1, r2)
            @assert i == 0
            r1[1], r2[1]
        end
    end
end

for name in ("pap",
          "pdp",
          "sepp")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(a, b)
            ccall(($fc,liberfa),Cdouble,(Ptr{Cdouble},Ptr{Cdouble}),a,b)
        end
    end
end

for name in ("pas",
          "seps")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(al,ap,bl,bp)
            ccall(($fc,liberfa),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),al,ap,bl,bp)
        end
    end
end

for name in ("rxp",
          "trxp")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(r, p)
            rp = zeros(3)
            ccall(($fc,liberfa), Void,
                  (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  r,p,rp)
            rp
        end
    end
end

for name in ("rx",
          "ry",
          "rz")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(a, r)
            ccall(($fc,liberfa),Void,
                  (Cdouble,Ptr{Cdouble}),
                  a,r)
            r
        end
    end
end

for name in ("rxpv",
          "trxpv")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(r, p)
            rp = zeros((2,3))
            ccall(($fc,liberfa), Void,
                  (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  r,p,rp)
            rp
        end
    end
end

for name in ("xys00a",
          "xys00b",
          "xys06a")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(date1, date2)
            x = [0.]
            y = [0.]
            s = [0.]
            ccall(($fc,liberfa), Void,
                  (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  date1,date2,x,y,s)
            x[1], y[1], s[1]
        end
    end
end

for name in ("ltecm",
          "ltp",
          "ltpb")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(epj)
            rp = zeros((3,3))
            ccall(($fc,liberfa), Void,
                  (Cdouble,Ptr{Cdouble}),
                  epj,rp)
            rp
        end
    end
end

for name in ("ltpecl",
          "ltpequ")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(epj)
            vec = zeros(3)
            ccall(($fc,liberfa), Void,
                  (Cdouble,Ptr{Cdouble}),
                  epj,vec)
            vec
        end
    end
end

for name in ("eceq06",
          "eqec06")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(date1,date2,d1,d2)
            r1 = [0.0]
            r2 = [0.0]
            ccall(($fc,liberfa), Void,
                  (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
                  date1, date2, d1, d2, r1,r2)
            r1[1], r2[1]
        end
    end
end

for name in ("lteceq",
          "lteqec")
    f = Symbol(name)
    fc = "era" * titlecase(name)
    @eval begin
        function ($f)(epj,d1,d2)
            r1 = [0.0]
            r2 = [0.0]
            ccall(($fc,liberfa), Void,
                  (Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
                  epj, d1, d2, r1,r2)
            r1[1], r2[1]
        end
    end
end

end # module
