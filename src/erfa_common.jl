const DPI = 3.141592653589793
const D2PI = 6.283185307179586
const DR2D = 57.29577951308232
const DD2R = 0.017453292519943295
const DR2AS = 206264.80624709636
const DAS2R = 4.84813681109536e-6
const DS2R = 7.27220521664304e-5
const TURNAS = 1.296e6
const DMAS2R = DAS2R / 1000.0
const DTY = 365.242198781
const DAYSEC = 86400.0
const DJY = 365.25
const DJC = 36525.0
const DJM = 365250.0
const DJ00 = 2.451545e6
const DJM0 = 2.4000005e6
const DJM00 = 51544.5
const DJM77 = 43144.0
const TTMTAI = 32.184
const DAU = 1.4959787e11
const CMPS = 2.99792458e8
const AULT = 499.004782
const DC = DAYSEC / AULT
const ELG = 6.969290134e-10
const ELB = 1.550519768e-8
const TDB0 = -6.55e-5
const SRS = 1.97412574336e-8

@enum Ellipsoid WGS84 = 1 GRS80 = 2 WGS72 = 3

mutable struct ASTROM
    pmt::Cdouble
    eb::NTuple{3,Cdouble}
    eh::NTuple{3,Cdouble}
    em::Cdouble
    v::NTuple{3,Cdouble}
    bm1::Cdouble
    bpn::NTuple{9,Cdouble}
    along::Cdouble
    phi::Cdouble
    xpl::Cdouble
    ypl::Cdouble
    sphi::Cdouble
    cphi::Cdouble
    diurab::Cdouble
    eral::Cdouble
    refa::Cdouble
    refb::Cdouble
end

struct LDBODY
    bm::Cdouble
    dl::Cdouble
    pv::NTuple{6,Cdouble}
end

struct ERFAExcpetion <: Exception
    msg::String
end

Base.showerror(io::IO, ex::ERFAExcpetion) = print(io, ex.msg)
