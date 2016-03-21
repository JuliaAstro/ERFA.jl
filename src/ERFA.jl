module ERFA

using Compat
import Base.getindex

export
    ERFA_DPI,
    ERFA_D2PI,
    ERFA_DR2D,
    ERFA_DD2R,
    ERFA_DR2AS,
    ERFA_DAS2R,
    ERFA_DS2R,
    ERFA_TURNAS,
    ERFA_DMAS2R,
    ERFA_DTY,
    ERFA_DAYSEC,
    ERFA_DJ,
    ERFA_DJC,
    ERFA_DJM,
    ERFA_DJ00,
    ERFA_DJM0,
    ERFA_DJM00,
    ERFA_DJM77,
    ERFA_TTMTAI,
    ERFA_DAU,
    ERFA_CMP,
    ERFA_AULT,
    ERFA_DC,
    ERFA_ELG,
    ERFA_ELB,
    ERFA_TDB0,
    ERFA_SRS,
    ERFA_WGS84,
    ERFA_GRS80,
    ERFA_WGS72

export
    eraASTROM,
    eraLDBODY,
    eraA2af,
    eraA2tf,
    eraAb,
    eraAf2a,
    eraAnp,
    eraAnpm,
    eraApcg,
    eraApcg13,
    eraApci,
    eraApci13,
    eraApco,
    eraApco13,
    eraApcs,
    eraApcs13,
    eraAper,
    eraAper13,
    eraApio,
    eraApio13,
    eraAtci13,
    eraAtciq,
    eraAtciqn,
    eraAtciqz,
    eraAtco13,
    eraAtic13,
    eraAticq,
    eraAticqn,
    eraAtio13,
    eraAtioq,
    eraAtoc13,
    eraAtoi13,
    eraAtoiq,
    eraBi00,
    eraBp00,
    eraBp06,
    eraBpn2xy,
    eraCal2jd,
    eraC2s,
    eraC2i00a,
    eraC2i00b,
    eraC2i06a,
    eraC2ibpn,
    eraC2ixy,
    eraC2ixys,
    eraC2t00a,
    eraC2t00b,
    eraC2t06a,
    eraC2tcio,
    eraC2teqx,
    eraC2tpe,
    eraC2txy,
    eraCp,
    eraCpv,
    eraCr,
    eraDat,
    eraD2dtf,
    eraD2tf,
    eraDtdb,
    eraDtf2d,
    eraEe00,
    eraEe00a,
    eraEe00b,
    eraEe06a,
    eraEect00,
    eraEform,
    eraEo06a,
    eraEors,
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
    eraFk52h,
    eraFk5hip,
    eraFk5hz,
    eraFw2m,
    eraFw2xy,
    eraG2icrs,
    eraGc2gd,
    eraGc2gde,
    eraGd2gc,
    eraGd2gce,
    eraGmst00,
    eraGmst06,
    eraGmst82,
    eraGst00a,
    eraGst00b,
    eraGst06,
    eraGst06a,
    eraGst94,
    eraH2fk5,
    eraHfk5z,
    eraIcrs2g,
    eraIr,
    eraJd2cal,
    eraJdcalf,
    eraLd,
    eraLdn,
    eraLdsun,
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
    eraP06e,
    eraP2pv,
    eraP2s,
    eraPap,
    eraPas,
    eraPb06,
    eraPdp,
    eraPfw06,
    eraPlan94,
    eraPm,
    eraPmp,
    eraPmat00,
    eraPmat06,
    eraPmat76,
    eraPmpx,
    eraPmsafe,
    eraPn,
    eraPn00,
    eraPn00a,
    eraPn00b,
    eraPn06,
    eraPn06a,
    eraPnm00a,
    eraPnm00b,
    eraPnm06a,
    eraPnm80,
    eraPom00,
    eraPpp,
    eraPpsp,
    eraPr00,
    eraPrec76,
    eraPv2p,
    eraPv2s,
    eraPvdpv,
    eraPvm,
    eraPvmpv,
    eraPvppv,
    eraPvstar,
    eraPvtob,
    eraPvu,
    eraPvup,
    eraPvxpv,
    eraPxp,
    eraRefco,
    eraRm2v,
    eraRv2m,
    eraRx,
    eraRy,
    eraRz,
    eraRxp,
    eraRxpv,
    eraRxr,
    eraS00,
    eraS00a,
    eraS00b,
    eraS06,
    eraS06a,
    eraS2c,
    eraS2p,
    eraS2pv,
    eraS2xpv,
    eraSp00,
    eraSepp,
    eraSeps,
    eraStarpm,
    eraStarpv,
    eraSxp,
    eraSxpv,
    eraTaitt,
    eraTaiut1,
    eraTaiutc,
    eraTcbtdb,
    eraTcgtt,
    eraTdbtcb,
    eraTdbtt,
    eraTf2a,
    eraTf2d,
    eraTr,
    eraTrxp,
    eraTrxpv,
    eraTttai,
    eraTttcg,
    eraTttdb,
    eraTtut1,
    eraUt1tai,
    eraUt1tt,
    eraUt1utc,
    eraUtctai,
    eraUtcut1,
    eraXy06,
    eraXys00a,
    eraXys00b,
    eraXys06a,
    eraZp,
    eraZpv,
    eraZr

include("../deps/deps.jl")
include("erfa_common.jl")

function getindex(A::Array_3_Cdouble, i::Integer)
    if i === 1
        return A.d1
    elseif i === 2
        return A.d2
    elseif i === 3
        return A.d3
    end
end

function getindex(A::Array_3_Array_3_Cdouble, i::Integer)
    if i === 1
        return getindex(A.d1,1)
    elseif i === 2
        return getindex(A.d1,2)
    elseif i === 3
        return getindex(A.d1,3)
    elseif i === 4
        return getindex(A.d2,1)
    elseif i === 5
        return getindex(A.d2,2)
    elseif i === 6
        return getindex(A.d2,3)
    elseif i === 7
        return getindex(A.d3,1)
    elseif i === 8
        return getindex(A.d3,2)
    elseif i === 9
        return getindex(A.d3,3)
    end
end

function eraASTROM(pmt::Float64, eb::Array{Float64}, eh::Array{Float64}, em::Float64, v::Array{Float64}, bm1::Float64, bpn::Array{Float64}, along::Float64, phi::Float64, xpl::Float64, ypl::Float64, sphi::Float64, cphi::Float64, diurab::Float64, eral::Float64, refa::Float64, refb::Float64) 
    eraASTROM(pmt,
              Array_3_Cdouble(eb[1], eb[2], eb[3]),
              Array_3_Cdouble(eh[1], eh[2], eh[3]),
              em,
              Array_3_Cdouble(v[1], v[2], v[3]),
              bm1,
              Array_3_Array_3_Cdouble(Array_3_Cdouble(bpn[1],bpn[2],bpn[3]),
                                      Array_3_Cdouble(bpn[4],bpn[5],bpn[6]),
                                      Array_3_Cdouble(bpn[7],bpn[8],bpn[9])),
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

function eraAb(pnat::Array{Cdouble},v::Array{Cdouble},s::Cdouble,bm1::Cdouble)
    ppr = zeros(3)
    ccall((:eraAb,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Cdouble,Cdouble,Ptr{Cdouble}),
          pnat,v,s,bm1,ppr)
    ppr
end

function eraApcg(date1::Cdouble,date2::Cdouble,ebpv::Array{Cdouble},ehp::Array{Cdouble})
    astrom = eraASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApcg,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{eraASTROM}),
          date1,date2,ebpv,ehp,&astrom)
    astrom
end

function eraApcg13(date1::Cdouble,date2::Cdouble)
    astrom = eraASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApcg13,liberfa),Void,
          (Cdouble,Cdouble,Ptr{eraASTROM}),
          date1,date2,&astrom)
    astrom
end

function eraApci(date1::Cdouble,date2::Cdouble,ebpv::Array{Cdouble},ehp::Array{Cdouble},x::Cdouble,y::Cdouble,s::Cdouble)
    astrom = eraASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApci,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Cdouble,Cdouble,Cdouble,Ptr{eraASTROM}),
          date1,date2,ebpv,ehp,x,y,s,&astrom)
    astrom
end

function eraApci13(date1::Cdouble,date2::Cdouble)
    astrom = eraASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    eo = [0.]
    ccall((:eraApci13,liberfa),Void,
          (Cdouble,Cdouble,Ptr{eraASTROM},Ptr{Cdouble}),
          date1,date2,&astrom,eo)
    astrom, eo[1]
end

function eraApco(date1::Cdouble,date2::Cdouble,ebpv::Array{Cdouble},ehp::Array{Cdouble},x::Cdouble,y::Cdouble,s::Cdouble,theta::Cdouble,elong::Cdouble,phi::Cdouble,hm::Cdouble,xp::Cdouble,yp::Cdouble,sp::Cdouble,refa::Cdouble,refb::Cdouble)
    astrom = eraASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApco,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{eraASTROM}),
          date1,date2,ebpv,ehp,x,y,s,theta,elong,phi,hm,xp,yp,sp,refa,refb,&astrom)
    astrom
end

function eraApco13(utc1::Cdouble,utc2::Cdouble,dut1::Cdouble,elong::Cdouble,phi::Cdouble,hm::Cdouble,xp::Cdouble,yp::Cdouble,phpa::Cdouble,tk::Cdouble,rh::Cdouble,wl::Cdouble)
    astrom = eraASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    eo = [0.]
    i = ccall((:eraApco13,liberfa),Cint,
              (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{eraASTROM},Ptr{Cdouble}),
              utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl,&astrom,eo)
    if i == -1
        error("unacceptable date")
    elseif i == +1
        warn("dubious year")
    end
    astrom, eo[1]
end

function eraApcs(date1::Cdouble,date2::Cdouble,pv::Array{Cdouble},ebpv::Array{Cdouble},ehp::Array{Cdouble})
    astrom = eraASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApcs,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{eraASTROM}),
          date1,date2,pv,ebpv,ehp,&astrom)
    astrom
end

function eraApcs13(date1::Cdouble,date2::Cdouble,pv::Array{Cdouble})
    astrom = eraASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApcs13,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{eraASTROM}),
          date1,date2,pv,&astrom)
    astrom
end

function eraAper(theta::Cdouble,astrom::eraASTROM)
    ccall((:eraAper,liberfa),Void,
          (Cdouble,Ptr{eraASTROM}),
          theta,&astrom)
    astrom
end

function eraAper13(ut11::Cdouble,ut12::Cdouble,astrom::eraASTROM)
    ccall((:eraAper13,liberfa),Void,
          (Cdouble,Cdouble,Ptr{eraASTROM}),
          ut11,ut12,&astrom)
    astrom
end

function eraApio(sp::Cdouble,theta::Cdouble,elong::Cdouble,phi::Cdouble,hm::Cdouble,xp::Cdouble,yp::Cdouble,refa::Cdouble,refb::Cdouble)
    astrom = eraASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    ccall((:eraApio,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{eraASTROM}),
          sp,theta,elong,phi,hm,xp,yp,refa,refb,&astrom)
    astrom
end

function eraApio13(utc1::Cdouble,utc2::Cdouble,dut1::Cdouble,elong::Cdouble,phi::Cdouble,hm::Cdouble,xp::Cdouble,yp::Cdouble,phpa::Cdouble,tk::Cdouble,rh::Cdouble,wl::Cdouble)    
    astrom = eraASTROM(0.0,zeros(3),zeros(3),0.0,zeros(3),0.0,zeros((3,3)),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    i = ccall((:eraApio13,liberfa),Cint,
              (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{eraASTROM}),
              utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl,&astrom)
    if i == -1
        error("unacceptable date")
    elseif i == +1
        warn("dubious year")
    end
    astrom
end

function eraAtci13(rc::Cdouble,dc::Cdouble,pr::Cdouble,pd::Cdouble,px::Cdouble,rv::Cdouble,date1::Cdouble,date2::Cdouble)
    ri = [0.]
    di = [0.]
    eo = [0.]
    ccall((:eraAtci13,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          rc,dc,pr,pd,px,rv,date1,date2,ri,di,eo)
    ri[1], di[1], eo[1]
end

function eraAtciq(rc::Cdouble,dc::Cdouble,pr::Cdouble,pd::Cdouble,px::Cdouble,rv::Cdouble,astrom::eraASTROM)
    ri = [0.]
    di = [0.]
    ccall((:eraAtciq,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{eraASTROM},Ptr{Cdouble},Ptr{Cdouble}),
          rc,dc,pr,pd,px,rv,&astrom,ri,di)
    ri[1], di[1]
end

function eraAtciqn(rc::Cdouble,dc::Cdouble,pr::Cdouble,pd::Cdouble,px::Cdouble,rv::Cdouble,astrom::eraASTROM,b::Array{eraLDBODY})
    ri = [0.]
    di = [0.]
    n = length(b)
    ccall((:eraAtciqn,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{eraASTROM},Cint,Ptr{eraLDBODY},Ptr{Cdouble},Ptr{Cdouble}),
          rc,dc,pr,pd,px,rv,&astrom,n,b,ri,di)
    ri[1], di[1]
end

function eraAtciqz(rc::Cdouble,dc::Cdouble,astrom::eraASTROM)
    ri = [0.]
    di = [0.]
    ccall((:eraAtciqz,liberfa),Void,
          (Cdouble,Cdouble,Ptr{eraASTROM},Ptr{Cdouble},Ptr{Cdouble}),
          rc,dc,&astrom,ri,di)
    ri[1], di[1]
end

function eraAtco13(rc::Cdouble,dc::Cdouble,pr::Cdouble,pd::Cdouble,px::Cdouble,rv::Cdouble,utc1::Cdouble,utc2::Cdouble,dut1::Cdouble,elong::Cdouble,phi::Cdouble,hm::Cdouble,xp::Cdouble,yp::Cdouble,phpa::Cdouble,tk::Cdouble,rh::Cdouble,wl::Cdouble)
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

function eraAtic13(ri::Cdouble,di::Cdouble,date1::Cdouble,date2::Cdouble)
    rc = [0.]
    dc = [0.]
    eo = [0.]
    ccall((:eraAtic13,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          ri,di,date1,date2,rc,dc,eo)
    rc[1],dc[1],eo[1]
end

function eraAticq(ri::Cdouble,di::Cdouble,astrom::eraASTROM)
    rc = [0.]
    dc = [0.]
    ccall((:eraAticq,liberfa),Void,
          (Cdouble,Cdouble,Ptr{eraASTROM},Ptr{Cdouble},Ptr{Cdouble}),
          ri,di,&astrom,rc,dc)
    rc[1],dc[1]
end

function eraAticqn(ri::Cdouble,di::Cdouble,astrom::eraASTROM,b::Array{eraLDBODY})
    rc = [0.]
    dc = [0.]
    n = length(b)
    ccall((:eraAticqn,liberfa),Void,
          (Cdouble,Cdouble,Ptr{eraASTROM},Cint,Ptr{eraLDBODY},Ptr{Cdouble},Ptr{Cdouble}),
          ri,di,&astrom,n,b,rc,dc)
    rc[1],dc[1]
end

function eraAtio13(ri::Cdouble,di::Cdouble,utc1::Cdouble,utc2::Cdouble,dut1::Cdouble,elong::Cdouble,phi::Cdouble,hm::Cdouble,xp::Cdouble,yp::Cdouble,phpa::Cdouble,tk::Cdouble,rh::Cdouble,wl::Cdouble)
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

function eraAtioq(ri::Cdouble,di::Cdouble,astrom::eraASTROM)
    aob = [0.]
    zob = [0.]
    hob = [0.]
    dob = [0.]
    rob = [0.]
    ccall((:eraAtioq,liberfa),Void,
          (Cdouble,Cdouble,Ptr{eraASTROM},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          ri,di,&astrom,aob,zob,hob,dob,rob)
    aob[1],zob[1],hob[1],dob[1],rob[1]
end

function eraAtoc13(typeofcoordinates::Char,ob1::Cdouble,ob2::Cdouble,utc1::Cdouble,utc2::Cdouble,dut1::Cdouble,elong::Cdouble,phi::Cdouble,hm::Cdouble,xp::Cdouble,yp::Cdouble,phpa::Cdouble,tk::Cdouble,rh::Cdouble,wl::Cdouble)
    rc = [0.]
    dc = [0.]
    if !(typeofcoordinates in ('R', 'r', 'H', 'h', 'A', 'a'))
        typeofcoordinates = 'A'
    end
    i = ccall((:eraAtoc13,liberfa),Cint,
              (Ptr{Char},Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
              &typeofcoordinates,ob1,ob2,utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl,rc,dc)
    if i == -1
        error("unacceptable date")
    elseif i == +1
        warn("dubious year")
    end
    rc[1],dc[1]
end

function eraAtoi13(typeofcoordinates::Char,ob1::Cdouble,ob2::Cdouble,utc1::Cdouble,utc2::Cdouble,dut1::Cdouble,elong::Cdouble,phi::Cdouble,hm::Cdouble,xp::Cdouble,yp::Cdouble,phpa::Cdouble,tk::Cdouble,rh::Cdouble,wl::Cdouble)
    ri = [0.]
    di = [0.]
    i = ccall((:eraAtoi13,liberfa),Cint,
              (Ptr{Char},Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
              &typeofcoordinates,ob1,ob2,utc1,utc2,dut1,elong,phi,hm,xp,yp,phpa,tk,rh,wl,ri,di)
    if i == -1
        error("unacceptable date")
    elseif i == +1
        warn("dubious year")
    end
    ri[1],di[1]
end

function eraAtoiq(typeofcoordinates::Char,ob1::Cdouble,ob2::Cdouble,astrom::eraASTROM)
    ri = [0.]
    di = [0.]
    ccall((:eraAtoiq,liberfa),
          Void,(Ptr{Char},Cdouble,Cdouble,Ptr{eraASTROM},Ptr{Cdouble},Ptr{Cdouble}),
          &typeofcoordinates,ob1,ob2,&astrom,ri,di)
    ri[1],di[1]
end

function eraBi00()
    dpsibi = [0.]
    depsbi = [0.]
    dra = [0.]
    ccall((:eraBi00,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          dpsibi,depsbi,dra)
    dpsibi[1], depsbi[1], dra[1]
end

function eraBpn2xy(rbpn::Array{Cdouble})
    x = [0.]
    y = [0.]
    ccall((:eraBpn2xy,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          rbpn,x,y)
    x[1], y[1]
end

function eraC2ibpn(date1::Cdouble,date2::Cdouble,rbpn::Array{Cdouble})
    rc2i = zeros((3,3))
    ccall((:eraC2ibpn,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          date1,date2,rbpn,rc2i)
    rc2i
end

function eraC2s(p::Array{Cdouble})
    theta = [0.]
    phi = [0.]
    ccall((:eraC2s,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          p,theta,phi)
    theta[1], phi[1]
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

function eraCp(p::Array{Cdouble})
    c = zeros(3)
    ccall((:eraCp,liberfa),Void,(Ptr{Cdouble},Ptr{Cdouble}),
          p,c)
    c
end

function eraCpv(pv::Array{Cdouble})
    c = zeros((2,3))
    ccall((:eraCpv,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          pv,c)
    c
end

function eraCr(p::Array{Cdouble})
    r = zeros((3,3))
    ccall((:eraCr,liberfa),Void,(Ptr{Cdouble},Ptr{Cdouble}),
          p,r)
    r
end

function eraDtdb(date1::Cdouble,date2::Cdouble,ut::Cdouble,elong::Cdouble,u::Cdouble,v::Cdouble)
    ccall((:eraDtdb,liberfa),Cdouble,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble),
          date1,date2,ut,elong,u,v)
end

function eraEform(n::Integer)
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

function eraEors(rnpb::Array{Cdouble},s::Cdouble)
    ccall((:eraEors,liberfa),Cdouble,
          (Ptr{Cdouble},Cdouble),
          rnpb,s)
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
    pvh = zeros((2,3))
    pvb = zeros((2,3))
    i = ccall((:eraEpv00, liberfa),
              Cint,
              (Float64, Float64, Ptr{Float64}, Ptr{Float64}),
              date1, date2, pvh, pvb)
    if i == 1
        warn("date outside the range 1900-2100 AD")
    end
    pvh, pvb
end

function eraFk5hip()
    r5h = zeros((3,3))
    s5h = zeros(3)
    ccall((:eraFk5hip,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          r5h,s5h)
    r5h,s5h
end

function eraFk5hz(r5::Cdouble,d5::Cdouble,date1::Cdouble,date2::Cdouble)
    rh = [0.]
    dh = [0.]
    ccall((:eraFk5hz,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          r5,d5,date1,date2,rh,dh)
    rh[1],dh[1]
end

function eraFw2xy(gamb::Cdouble,phib::Cdouble,psi::Cdouble,eps::Cdouble)
    x =[0.]
    y = [0.]
    ccall((:eraFw2xy,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          gamb,phib,psi,eps,x,y)
    x[1], y[1]
end

function eraGc2gd(n::Integer,xyz::Array{Cdouble})
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

function eraGc2gde(a::Cdouble,f::Cdouble,xyz::Array{Cdouble})
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

function eraGd2gc(n::Integer,elong::Cdouble,phi::Cdouble,height::Cdouble)
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

function eraGd2gce(a::Cdouble,f::Cdouble,elong::Cdouble,phi::Cdouble,height::Cdouble)
    xyz = zeros(3)
    i = ccall((:eraGd2gce,liberfa),Cint,
              (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
              a,f,elong,phi,height,xyz)
    if i == -1
        error("illegal case")
    end
    xyz
end

function eraGst06(uta::Cdouble,utb::Cdouble,tta::Cdouble,ttb::Cdouble,rnpb::Array{Cdouble})
    ccall((:eraGst06,liberfa),Cdouble,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
          uta,utb,tta,ttb,rnpb)
end

function eraHfk5z(rh::Cdouble,dh::Cdouble,date1::Cdouble,date2::Cdouble)
    r5 = [0.]
    d5 = [0.]
    dr5 = [0.]
    dd5 = [0.]
    ccall((:eraHfk5z,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          rh,dh,date1,date2,r5,d5,dr5,dd5)
    r5[1],d5[1],dr5[1],dd5[1]
end

function eraLDBODY(bm::Cdouble, dl::Cdouble, pv::Array{Float64})
    p = Array_3_Cdouble(pv[1], pv[2], pv[3])
    v = Array_3_Cdouble(pv[4], pv[5], pv[6])
    eraLDBODY(bm, dl, Array_2_Array_3_Cdouble(p, v))
end

function eraLd(bm::Cdouble,p::Array{Cdouble},q::Array{Cdouble},e::Array{Cdouble},em::Cdouble,dlim::Cdouble)
    p1 = zeros(3)
    ccall((:eraLd,liberfa),Void,
          (Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Cdouble,Cdouble,Ptr{Cdouble}),
          bm,p,q,e,em,dlim,p1)
    p1
end

function eraLdn(l::Array{eraLDBODY}, ob::Array{Float64}, sc::Array{Float64})
    sn = zeros(3)
    n = length(l)
    ccall((:eraLdn, liberfa),Void,
          (Cint, Ptr{eraLDBODY}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
          n, l, ob, sc, sn)
    sn
end

function eraLdsun(p::Array{Cdouble},e::Array{Cdouble},em::Cdouble)
    p1 = zeros(3)
    ccall((:eraLdsun,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Cdouble,Ptr{Cdouble}),
          p,e,em,p1)
    p1
end

function eraPmpx(rc::Cdouble,dc::Cdouble,pr::Cdouble,pd::Cdouble,px::Cdouble,rv::Cdouble,pmt::Cdouble,vob::Array{Cdouble})
    pco = zeros(3)
    ccall((:eraPmpx,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          rc,dc,pr,pd,px,rv,pmt,vob,pco)
    pco
end

function eraNumat(epsa::Real, dpsi::Real, deps::Real)
    rmatn = zeros((3,3))
    ccall((:eraNumat,liberfa),Void,
          (Float64, Float64, Float64, Ptr{Float64}),
          epsa, dpsi, deps, rmatn)
    rmatn
end

function eraP06e(date1::Cdouble,date2::Cdouble)
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

function eraP2s(p::Array{Cdouble})
    theta = [0.]
    phi = [0.]
    r = [0.]
    ccall((:eraP2s,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          p,theta,phi,r)
    theta[1], phi[1], r[1]
end

function eraP2pv(p::Array{Cdouble})
    pv = zeros((2,3))
    ccall((:eraP2pv,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          p,pv)
    pv
end

function eraPb06(date1::Cdouble,date2::Cdouble)
    bzeta = [0.]
    bz = [0.]
    btheta = [0.]
    ccall((:eraPb06,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          date1,date2,bzeta,bz,btheta)
    bzeta[1],bz[1],btheta[1]
end


function eraPfw06(date1::Cdouble,date2::Cdouble)
    gamb = [0.]
    phib = [0.]
    psib = [0.]
    epsa = [0.]
    ccall((:eraPfw06,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          date1,date2,gamb,phib,psib,epsa)
    gamb[1],phib[1],psib[1],epsa[1]
end

function eraPlan94(date1::Float64, date2::Float64, np::Integer)
    pv = zeros((2,3))
    i = ccall((:eraPlan94, liberfa),Cint,
              (Float64, Float64, Cint, Ptr{Float64}),
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

function eraPm(p::Array{Cdouble})
    ccall((:eraPm,liberfa),Cdouble,(Ptr{Cdouble},),p)
end

function eraPmsafe(ra1::Cdouble,dec1::Cdouble,pmr1::Cdouble,pmd1::Cdouble,px1::Cdouble,rv1::Cdouble,ep1a::Cdouble,ep1b::Cdouble,ep2a::Cdouble,ep2b::Cdouble)
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

function eraPn(p::Array{Cdouble})
    r = [0.]
    u = zeros(3)
    ccall((:eraPn,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          p,r,u)
    r[1],u
end

function eraPpsp(a::Array{Cdouble},s::Cdouble,b::Array{Cdouble})
    apsb = zeros(3)
    ccall((:eraPpsp,liberfa),Void,
          (Ptr{Cdouble},Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          a,s,b,apsb)
    apsb
end

function eraPrec76(ep01::Cdouble,ep02::Cdouble,ep11::Cdouble,ep12::Cdouble)
    zeta = [0.]
    z = [0.]
    theta = [0.]
    ccall((:eraPrec76,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          ep01,ep02,ep11,ep12,zeta,z,theta)
    zeta[1],z[1],theta[1]
end

function eraPv2s(pv::Array{Cdouble})
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

function eraPv2p(pv::Array{Cdouble})
    p = zeros(3)
    ccall((:eraPv2p,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          pv,p)
    p
end

function eraPvdpv(a::Array{Cdouble},b::Array{Cdouble})
    adb = zeros(2)
    ccall((:eraPvdpv,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          a,b,adb)
    adb
end

function eraPvm(pv::Array{Cdouble})
    s = [0.]
    r = [0.]
    ccall((:eraPvm,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          pv,r,s)
    r[1], s[1]
end

function eraPvstar(pv::Array{Cdouble})
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

function eraPvtob(elong::Cdouble,phi::Cdouble,height::Cdouble,xp::Cdouble,yp::Cdouble,sp::Cdouble,theta::Cdouble)
    pv = zeros((2,3))
    ccall((:eraPvtob,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
          elong,phi,height,xp,yp,sp,theta,pv)
    pv
end

function eraRefco(phpa::Cdouble,tk::Cdouble,rh::Cdouble,wl::Cdouble)
    refa = [0.]
    refb = [0.]
    ccall((:eraRefco,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          phpa,tk,rh,wl,refa,refb)
    refa[1],refb[1]
end

function eraPvu(dt::Cdouble,pv::Array{Cdouble})
    upv = zeros((2,3))
    ccall((:eraPvu,liberfa),Void,
          (Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          dt,pv,upv)
    upv
end

function eraPvup(dt::Cdouble,pv::Array{Cdouble})
    p = zeros(3)
    ccall((:eraPvup,liberfa),Void,
          (Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          dt,pv,p)
    p
end

function eraRm2v(r::Array{Cdouble})
    w = zeros(3)
    ccall((:eraRm2v,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          r,w)
    w
end

function eraRv2m(w::Array{Cdouble})
    r = zeros((3,3))
    ccall((:eraRv2m,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          w,r)
    r
end

function eraRxr(a::Array{Cdouble},b::Array{Cdouble})
    atb = zeros((3,3))
    ccall((:eraRxr,liberfa),
          Void,
          (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
          a,b,atb)
    atb
end

function eraS2c(theta::Cdouble,phi::Cdouble)
    c = zeros(3)
    ccall((:eraS2c,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble}),
          theta,phi,c)
    c
end

function eraS2p(theta::Cdouble,phi::Cdouble,r::Cdouble)
    p = zeros(3)
    ccall((:eraS2p,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
          theta,phi,r,p)
    p
end

function eraS2pv(theta::Cdouble,phi::Cdouble,r::Cdouble,td::Cdouble,pd::Cdouble,rd::Cdouble)
    pv = zeros((2,3))
    ccall((:eraS2pv,liberfa),Void,
          (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
          theta,phi,r,td,pd,rd,pv)
    pv
end

function eraS2xpv(s1::Cdouble,s2::Cdouble,pv::Array{Cdouble})
    spv = zeros((2,3))
    ccall((:eraS2xpv,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          s1,s2,pv,spv)
    spv
end

function eraStarpm(ra1::Cdouble,dec1::Cdouble,pmr1::Cdouble,pmd1::Cdouble,px1::Cdouble,rv1::Cdouble,ep1a::Cdouble,ep1b::Cdouble,ep2a::Cdouble,ep2b::Cdouble)
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

function eraStarpv(ra::Cdouble,dec::Cdouble,pmr::Cdouble,pmd::Cdouble,px::Cdouble,rv::Cdouble)
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

function eraSxp(s::Cdouble,p::Array{Cdouble})
    sp = zeros(3)
    ccall((:eraSxp,liberfa),Void,
          (Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          s,p,sp)
    sp
end

function eraSxpv(s::Cdouble,pv::Array{Cdouble})
    spv = zeros((2,3))
    ccall((:eraSxpv,liberfa),Void,
          (Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          s,pv,spv)
    spv
end

function eraTr(r::Array{Cdouble})
    rt = zeros((3,3))
    ccall((:eraTr,liberfa),Void,
          (Ptr{Cdouble},Ptr{Cdouble}),
          r,rt)
    rt
end

function eraXy06(date1::Real,date2::Real)
    x = [0.]
    y = [0.]
    ccall((:eraXy06,liberfa),Void,
          (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
          date1,date2,x,y)
    x[1], y[1]
end

for f in (:eraAf2a,
          :eraTf2a,
          :eraTf2d)
    @eval begin
        function ($f)(s::Char,ideg::Integer,iamin::Integer,asec::Real)
            rad = [0.]
            i = ccall(($(Expr(:quote,f)),liberfa),Cint,
                       (Cchar,Cint,Cint,Cdouble,Ptr{Cdouble}),
                       s,ideg,iamin,asec,rad)
            @assert i == 0
            rad[1]
        end
    end              
end

function eraIr()
    r = zeros((3,3))
    ccall((:eraIr,liberfa),Void,(Ptr{Cdouble},),r)
    r
end

function eraZp()
    p = zeros(3)
    #ccall((:eraZp,liberfa),Void,(Ptr{Cdouble},),p)
    p
end

function eraZpv()
    pv = zeros((2,3))
    #ccall((:eraZpv,liberfa),Void,(Ptr{Cdouble},),pv)
    pv
end

function eraZr()
    r = zeros((3,3))
    #ccall((:eraZr,liberfa),Void,(Ptr{Cdouble},),r)
    r
end

for f in (:eraA2af,
          :eraA2tf,
          :eraD2tf)
    @eval begin
        function ($f)(ndp::Integer, a::Float64)
            s = Array(@compat(UInt8), 1)
            s[1] = '+'
            i = Int32[0, 0, 0, 0]
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Cint, Float64, Ptr{UInt8}, Ptr{Cint}),
                  ndp, a, s, i)
            @compat(Char(s[1])), i[1], i[2], i[3], i[4]
        end
    end
end

for f in (:eraAnp,
          :eraAnpm)
    @eval begin
        function ($f)(a::Cdouble)
            ccall(($(Expr(:quote,f)),liberfa),Cdouble,(Cdouble,),a)
        end
    end
end

for f in (:eraBp00,
          :eraBp06)
    @eval begin
        function ($f)(a::Float64, b::Float64)
            rb = zeros((3,3))
            rp = zeros((3,3))
            rbp = zeros((3,3))
            ccall(($(Expr(:quote,f)),liberfa),
                  Void,
                  (Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
                  a, b, rb, rp, rbp)
            rb, rp, rbp
        end
    end
end


for f in (:eraC2ixys,
          :eraPom00)
    @eval begin
        function ($f)(x::Cdouble, y::Cdouble, s::Cdouble)
            r = zeros((3,3))
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
                  x, y, s, r)
            r
        end
    end
end

for f in (:eraC2tcio,
          :eraC2teqx)
    @eval begin
        function ($f)(rc2i::Array{Cdouble}, era::Cdouble, rpom::Array{Cdouble})
            rc2t = zeros((3,3))
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Ptr{Cdouble},Cdouble,Ptr{Cdouble},Ptr{Cdouble}),
                  rc2i,era,rpom,rc2t)
            rc2t
        end
    end
end

for f in (:eraC2ixy,
          :eraFw2m)
    @eval begin
        function ($f)(x::Cdouble, y::Cdouble, s::Cdouble, t::Cdouble)
            r = zeros((3,3))
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
                  x, y, s, t, r)
            r
        end
    end
end

for f in (:eraC2t00a,
          :eraC2t00b,
          :eraC2t06a)
    @eval begin
        function ($f)(tta::Cdouble,ttb::Cdouble,uta::Cdouble,utb::Cdouble,xp::Cdouble,yp::Cdouble)
            rc2t =zeros((3,3))
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
                  tta,ttb,uta,utb,xp,yp,rc2t)
            rc2t
        end
    end
end


for f in (:eraC2tpe,
          :eraC2txy)
    @eval begin
        function ($f)(tta::Cdouble,ttb::Cdouble,uta::Cdouble,utb::Cdouble,x::Cdouble,y::Cdouble,xp::Cdouble,yp::Cdouble)
            rc2t =zeros((3,3))
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble}),
                  tta,ttb,uta,utb,x,y,xp,yp,rc2t)
            rc2t
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
            r = zeros((3,3))
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Float64, Float64, Ptr{Float64}),
                  a, b, r)
            r
        end
    end
end

for f in (:eraPn00,
          :eraPn06)
    @eval begin
        function ($f)(date1::Float64, date2::Float64, dpsi::Float64, deps::Float64)
            epsa = [0.]
            rb = zeros((3,3))
            rp = zeros((3,3))
            rbp = zeros((3,3))
            rn = zeros((3,3))
            rbpn = zeros((3,3))
            ccall(($(Expr(:quote,f)), liberfa),Void,
                  (Float64, Float64, Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
                  date1, date2, dpsi, deps, epsa, rb, rp, rbp, rn, rbpn)
            epsa[1], rb, rp, rbp, rn, rbpn
        end
    end
end

for f in (:eraPmp,
          :eraPpp,
          :eraPxp)
    @eval begin
        function ($f)(a::Array{Cdouble},b::Array{Cdouble})
            ab = zeros(3)
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  a,b,ab)
            ab
        end
    end
end

for f in (:eraPvmpv,
          :eraPvppv,
          :eraPvxpv)
    @eval begin
        function ($f)(a::Array{Cdouble},b::Array{Cdouble})
            ab = zeros((2,3))
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  a,b,ab)
            ab
        end
    end
end

for f in (:eraPn00a,
          :eraPn00b,
          :eraPn06a)
    @eval begin
        function ($f)(date1::Float64, date2::Float64)
            dpsi = [0.]
            deps = [0.]
            epsa = [0.]
            rb = zeros((3,3))
            rp = zeros((3,3))
            rbp = zeros((3,3))
            rn = zeros((3,3))
            rbpn = zeros((3,3))
            ccall(($(Expr(:quote,f)), liberfa),Void,
                  (Float64, Float64, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}),
                  date1, date2, dpsi, deps, epsa, rb, rp, rbp, rn, rbpn)
            dpsi[1], deps[1], epsa[1], rb, rp, rbp, rn, rbpn
        end
    end
end

for f in (:eraNut00a,
          :eraNut00b,
          :eraNut06a,
          :eraNut80,
          :eraPr00,
          :eraG2icrs,
          :eraIcrs2g)
    @eval begin
        function ($f)(a::Float64, b::Float64)
            r1 = [0.]
            r2 = [0.]
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Float64, Float64, Ptr{Float64}, Ptr{Float64}),
                  a, b, r1, r2)
            r1[1], r2[1]
        end
    end
end

for f in (:eraFk52h,
          :eraH2fk5)
    @eval begin
        function ($f)(ra::Cdouble,dec::Cdouble,dra::Cdouble,ddec::Cdouble,px::Cdouble,rv::Cdouble)
            r = [0.]
            d = [0.]
            dr = [0.]
            dd = [0.]
            p = [0.]
            v = [0.]
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  ra,dec,dra,ddec,px,rv,r,d,dr,dd,p,v)
            r[1],d[1],dr[1],dd[1],p[1],v[1]
        end
    end
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

for f in (:eraEe00,
          :eraGmst00,
          :eraGmst06,
          :eraGst00a,
          :eraGst06a,
          :eraS00,
          :eraS06)
    @eval ($f)(d1::Real,d2::Real,t1::Real,t2::Real) = ccall(($(Expr(:quote,f)),liberfa), Float64, (Float64,Float64,Float64,Float64), d1,d2,t1,t2)
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
          :eraObl80,
          :eraS00a,
          :eraS00b,
          :eraS06a,
          :eraSp00)
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

for f in (:eraPap,
          :eraPdp,
          :eraSepp)
    @eval begin
        function ($f)(a::Array{Cdouble}, b::Array{Cdouble})
            ccall(($(Expr(:quote,f)),liberfa),Cdouble,(Ptr{Cdouble},Ptr{Cdouble}),a,b)
        end
    end
end

for f in (:eraPas,
          :eraSeps)
    @eval begin
        function ($f)(al::Cdouble,ap::Cdouble,bl::Cdouble,bp::Cdouble)
            ccall(($(Expr(:quote,f)),liberfa),Cdouble,(Cdouble,Cdouble,Cdouble,Cdouble),al,ap,bl,bp)
        end
    end
end

for f in (:eraRxp,
          :eraTrxp)
    @eval begin
        function ($f)(r::Array{Cdouble}, p::Array{Cdouble})
            rp = zeros(3)
            ccall(($(Expr(:quote,f)),liberfa), Void,
                  (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  r,p,rp)
            rp
        end
    end
end

for f in (:eraRx,
          :eraRy,
          :eraRz)
    @eval begin
        function ($f)(a::Cdouble, r::Array{Cdouble})
            ccall(($(Expr(:quote,f)),liberfa),Void,
                  (Cdouble,Ptr{Cdouble}),
                  a,r)
            r
        end
    end
end

for f in (:eraRxpv,
          :eraTrxpv)
    @eval begin
        function ($f)(r::Array{Cdouble}, p::Array{Cdouble})
            rp = zeros((2,3))
            ccall(($(Expr(:quote,f)),liberfa), Void,
                  (Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  r,p,rp)
            rp
        end
    end
end

for f in (:eraXys00a,
          :eraXys00b,
          :eraXys06a)
    @eval begin
        function ($f)(date1::Float64, date2::Float64)
            x = [0.]
            y = [0.]
            s = [0.]
            ccall(($(Expr(:quote,f)),liberfa), Void,
                  (Cdouble,Cdouble,Ptr{Cdouble},Ptr{Cdouble},Ptr{Cdouble}),
                  date1,date2,x,y,s)
            x[1], y[1], s[1]
        end
    end
end

end # module
