export ERFAException

"""
    DPI

Pi
"""
const DPI = 3.141592653589793

"""
    D2PI

2Pi
"""
const D2PI = 6.283185307179586

"""
    DR2D

Radians to degrees
"""
const DR2D = 57.29577951308232

"""
    DD2R

Degrees to radians
"""
const DD2R = 0.017453292519943295

"""
    DR2AS

Radians to arcseconds
"""
const DR2AS = 206264.80624709636

"""
    DAS2R

Arcseconds to radians
"""
const DAS2R = 4.84813681109536e-6

"""
    DS2R

Seconds of time to radians
"""
const DS2R = 7.27220521664304e-5

"""
    TURNAS

Arcseconds in a full circle
"""
const TURNAS = 1.296e6

"""
    DMAS2R

Milliarcseconds to radians
"""
const DMAS2R = DAS2R / 1000.0

"""
    DTY

Length of tropical year B1900 (days)
"""
const DTY = 365.242198781

"""
    DAYSEC

Seconds per day
"""
const DAYSEC = 86400.0

"""
    DJY

Days per Julian year
"""
const DJY = 365.25

"""
    DJC

Days per Julian century
"""
const DJC = 36525.0

"""
    DJM

Days per Julian millennium
"""
const DJM = 365250.0

"""
    DJ00

Reference epoch (J2000.0), Julian Date
"""
const DJ00 = 2.451545e6

"""
    DJM0

Julian Date of Modified Julian Date zero
"""
const DJM0 = 2.4000005e6

"""
    DJM00

Reference epoch (J2000.0), Modified Julian Date
"""
const DJM00 = 51544.5

"""
    DJM77

1977 Jan 1.0 as MJD
"""
const DJM77 = 43144.0

"""
    TTMTAI

TT minus TAI (s)
"""
const TTMTAI = 32.184

"""
    DAU

Astronomical unit (m, IAU 2012)
"""
const DAU = 1.4959787e11

"""
    CMPS

Speed of light (m/s)
"""
const CMPS = 2.99792458e8

"""
    AULT

Light time for 1 au (s)
"""
const AULT = 499.004782

"""
    DC

Speed of light (au per day)
"""
const DC = DAYSEC / AULT

"""
    ELG

L_G = 1 - d(TT)/d(TCG)
"""
const ELG = 6.969290134e-10

"""
    L_B

L_B = 1 - d(TDB)/d(TCB)
"""
const ELB = 1.550519768e-8

"""
    TDB0

TDB (s) at TAI 1977/1/1.0
"""
const TDB0 = -6.55e-5

"""
    SRS

Schwarzschild radius of the Sun (au) =
2 * 1.32712440041e20 / (2.99792458e8)^2 / 1.49597870700e11
"""
const SRS = 1.97412574336e-8

"""
    Ellipsoid

Reference ellipsoids

- [`WGS84`](@ref)
- [`GRS80`](@ref)
- [`WGS72`](@ref)
"""
@enum Ellipsoid WGS84=1 GRS80=2 WGS72=3

"""
    WGS84

Enum constant for the WGS84 reference ellipsoid
"""
WGS84

"""
    GRS80

Enum constant for the GRS80 reference ellipsoid
"""
GRS80

"""
    WGS72

Enum constant for the WGS72 reference ellipsoid
"""
WGS72

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

function ASTROM()
    ASTROM(
        0.0,
        zeros(Cdouble, 3),
        zeros(Cdouble, 3),
        0.0,
        zeros(Cdouble, 3),
        0.0,
        zeros(Cdouble, 3, 3),
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )
end

function Base.getproperty(a::ASTROM, field::Symbol)
    val = getfield(a, field)
    if field in (:eb, :eh, :v)
        return collect(val)
    elseif field == :bpn
        return permutedims(reshape(collect(val), 3, 3))
    end
    return val
end

struct LDBODY
    bm::Cdouble
    dl::Cdouble
    pv::NTuple{6,Cdouble}
end

struct ERFAException <: Exception
    msg::String
end

Base.showerror(io::IO, ex::ERFAException) = print(io, ex.msg)

macro checkdims(m::Int, n::Int, arr::Symbol...)
    ex = :()
    for a in arr
        name = string(a)
        expr = quote
            m1, n1 = size($(esc(a)))
            if (m1, n1) != ($m, $n)
                msg = string("`$($name)` must be a $($m)x$($n) matrix but is ", m1, "x", n1, ".")
                throw(ArgumentError(msg))
            end
        end
        push!(ex.args, expr)
    end
    :($ex; nothing)
end

macro checkdims(len::Int, arr::Symbol...)
    ex = :()
    for a in arr
        name = string(a)
        expr = quote
            n = length($(esc(a)))
            if n != $(esc(len))
                throw(ArgumentError("`$($name)` must have $($(esc(len))) elements but has $n."))
            end
        end
        push!(ex.args, expr)
    end
    :($ex; nothing)
end

function array_to_cmatrix(array; n=0)
    any(length.(array) .!= n) && throw(ArgumentError("Expected each vector to have $n elements."))
    return hcat(array...)
end

function array_to_cmatrix(array::Vector{Vector{Int}}; n=0)
    any(length.(array) .!= n) && throw(ArgumentError("Expected each vector to have $n elements."))
    array = map(x->Cint.(x), array)
    return hcat(array...)
end

function cmatrix_to_array(matrix)
    return [matrix[:, i] for i in 1:size(matrix, 2)]
end

function cmatrix_to_array(matrix::Matrix{Cint})
    return [Int.(matrix[:, i]) for i in 1:size(matrix, 2)]
end


