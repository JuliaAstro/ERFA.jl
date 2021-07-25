"""
    tpors(xi, eta, ra, dec)

In the tangent plane projection, given the rectangular coordinates
of a star and its spherical coordinates, determine the spherical
coordinates of the tangent point.

### Given ###

- `xi`, `eta`: rectangular coordinates of star image (Note 2)
- `a`, `b`: star's spherical coordinates (Note 3)

### Returned ###

- `status`: number of solutions:
    - 0 = no solutions returned (Note 5)
    - 1 = only the first solution is useful (Note 6)
    - 2 = both solutions are useful (Note 6)
- `a01`, `b01`: tangent point's spherical coordinates, Soln. 1
- `a02`, `b02`: tangent point's spherical coordinates, Soln. 2

### Notes ###

1. The tangent plane projection is also called the "gnomonic
   projection" and the "central projection".

2. The eta axis points due north in the adopted coordinate system.
   If the spherical coordinates are observed (RA,Dec), the tangent
   plane coordinates (xi,eta) are conventionally called the
   "standard coordinates".  If the spherical coordinates are with
   respect to a right-handed triad, (xi,eta) are also right-handed.
   The units of (xi,eta) are, effectively, radians at the tangent
   point.

3. All angular arguments are in radians.

4. The angles a01 and a02 are returned in the range 0-2pi.  The
   angles b01 and b02 are returned in the range +/-pi, but in the
   usual, non-pole-crossing, case, the range is +/-pi/2.

5. Cases where there is no solution can arise only near the poles.
   For example, it is clearly impossible for a star at the pole
   itself to have a non-zero xi value, and hence it is meaningless
   to ask where the tangent point would have to be to bring about
   this combination of xi and dec.

6. Also near the poles, cases can arise where there are two useful
   solutions.  The return value indicates whether the second of the
   two solutions returned is useful;  1 indicates only one useful
   solution, the usual case.

7. The basis of the algorithm is to solve the spherical triangle PSC,
   where P is the north celestial pole, S is the star and C is the
   tangent point.  The spherical coordinates of the tangent point are
   [a0,b0];  writing rho^2 = (xi^2+eta^2) and r^2 = (1+rho^2), side c
   is then (pi/2-b), side p is sqrt(xi^2+eta^2) and side s (to be
   found) is (pi/2-b0).  Angle C is given by sin(C) = xi/rho and
   cos(C) = eta/rho.  Angle P (to be found) is the longitude
   difference between star and tangent point (a-a0).

8. This function is a member of the following set:

   |   spherical  |   vector   |   solve for |
   |:-------------|:-----------|:------------|
   |   eraTpxes   |  eraTpxev  |      xi,eta |
   |   eraTpsts   |  eraTpstv  |        star |
   | > eraTpors < |  eraTporv  |      origin |

### Called ###

- [`anp`](@ref) normalize angle into range 0 to 2pi

### References ###

- Calabretta M.R. & Greisen, E.W., 2002, "Representations of
  celestial coordinates in FITS", Astron.Astrophys. 395, 1077

- Green, R.M., "Spherical Astronomy", Cambridge University Press,
  1987, Chapter 13.
"""
function tpors(xi, eta, ra, dec)
    az1 = Ref{Cdouble}()
    az2 = Ref{Cdouble}()
    bz1 = Ref{Cdouble}()
    bz2 = Ref{Cdouble}()
    status = ccall((:eraTpors, liberfa), Cint,
                   (Cdouble, Cdouble, Cdouble, Cdouble,
                    Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
                   xi, eta, ra, dec, az1, az2, bz1, bz2)
    return status, az1[], az2[], bz1[], bz2[]
end

"""
    tporv(xi, eta, v)

In the tangent plane projection, given the rectangular coordinates
of a star and its direction cosines, determine the direction
cosines of the tangent point.

### Given ###

- `xi`, `eta`: rectangular coordinates of star image (Note 2)
- `v`: star's direction cosines (Note 3)

### Returned ###

- `int`: number of solutions:
    - 0 = no solutions returned (Note 4)
    - 1 = only the first solution is useful (Note 5)
    - 2 = both solutions are useful (Note 5)
- `v01`: tangent point's direction cosines, Solution 1
- `v02`: tangent point's direction cosines, Solution 2

### Notes ###

1. The tangent plane projection is also called the "gnomonic
   projection" and the "central projection".

2. The eta axis points due north in the adopted coordinate system.
   If the direction cosines represent observed (RA,Dec., the tangent
   plane coordinates (xi,eta. are conventionally called the
   "standard coordinates".  If the direction cosines are with
   respect to a right-handed triad, (xi,eta. are also right-handed.
   The units of (xi,eta. are, effectively, radians at the tangent
   point.

3. The vector v must be of unit length or the result will be wrong.

4. Cases where there is no solution can arise only near the poles.
   For example, it is clearly impossible for a star at the pole
   itself to have a non-zero xi value, and hence it is meaningless
   to ask where the tangent point would have to be.

5. Also near the poles, cases can arise where there are two useful
   solutions.  The return value indicates whether the second of the
   two solutions returned is useful;  1 indicates only one useful
   solution, the usual case.

6. The basis of the algorithm is to solve the spherical triangle
   PSC, where P is the north celestial pole, S is the star and C is
   the tangent point.  Calling the celestial spherical coordinates
   of the star and tangent point (a,b. and (a0,b0) respectively, and
   writing rho^2 = (xi^2+eta^2. and r^2 = (1+rho^2), and
   transforming the vector v into (a,b. in the normal way, side c is
   then (pi/2-b., side p is sqrt(xi^2+eta^2) and side s (to be
   found. is (pi/2-b0), while angle C is given by sin(C) = xi/rho
   and cos(C. = eta/rho;  angle P (to be found) is (a-a0).  After
   solving the spherical triangle, the result (a0,b0. can be
   expressed in vector form as v0.

7. This function is a member of the following set:

     | spherical |    vector      |  solve for |
     |:----------|:---------------|:-----------|
     | eraTpxes  |   eraTpxev     |   xi,eta   |
     | eraTpsts  |   eraTpstv     |    star    |
     | eraTpors  | > eraTporv <   |   origin   |

### References ###

- Calabretta M.R. & Greisen, E.W., 2002, "Representations of
  celestial coordinates in FITS", Astron.Astrophys. 395, 1077

- Green, R.M., "Spherical Astronomy", Cambridge University Press,
  1987, Chapter 13.
"""
function tporv(xi, eta, v)
    vz1 = Array{Cdouble}(undef, 3)
    vz2 = Array{Cdouble}(undef, 3)
    status = ccall((:eraTporv, liberfa), Cint,
                   (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
                   xi, eta, v, vz1, vz2)
    return status, vz1, vz2
end

"""
    tpsts(xi, eta, raz, decz)

In the tangent plane projection, given the star's rectangular
coordinates and the spherical coordinates of the tangent point,
solve for the spherical coordinates of the star.

### Given ###

- `xi`, `eta`: rectangular coordinates of star image (Note 2)
- `a0`, `b0`: tangent point's spherical coordinates

### Returned ###

- `a`, `b`: star's spherical coordinates

### Notes ###

1. The tangent plane projection is also called the "gnomonic
   projection" and the "central projection".

2. The eta axis points due north in the adopted coordinate system.
   If the spherical coordinates are observed (RA,Dec), the tangent
   plane coordinates (xi,eta) are conventionally called the
   "standard coordinates".  If the spherical coordinates are with
   respect to a right-handed triad, (xi,eta) are also right-handed.
   The units of (xi,eta) are, effectively, radians at the tangent
   point.

3. All angular arguments are in radians.

4. This function is a member of the following set:

   |   spherical       | vector           | solve for |
   |:------------------|:-----------------|:----------|
   |  [`tpxes`](@ref)  | [`tpxev`](@ref)  | xi,eta    |
   | *[`tpsts`](@ref)* | [`tpstv`](@ref)  | star      |
   |  [`tpors`](@ref)  | [`tporv`](@ref)  | origin    |

### Called ###

- [`anp`](@ref) normalize angle into range 0 to 2pi

### References ###

- Calabretta M.R. & Greisen, E.W., 2002, "Representations of
  celestial coordinates in FITS", Astron.Astrophys. 395, 1077

- Green, R.M., "Spherical Astronomy", Cambridge University Press,
  1987, Chapter 13.
"""
function tpsts(xi, eta, raz, decz)
    ra = Ref{Cdouble}()
    dec = Ref{Cdouble}()
    ccall((:eraTpsts, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
          xi, eta, raz, decz, ra, dec)
    return ra[], dec[]
end

"""
    tpstv(xi, eta, vz)

In the tangent plane projection, given the star's rectangular
coordinates and the direction cosines of the tangent point, solve
for the direction cosines of the star.

### Given ###

- `xi`, `eta`: rectangular coordinates of star image (Note 2)
- `v0`: tangent point's direction cosines

### Returned ###

- `v`: star's direction cosines

### Notes ###

1. The tangent plane projection is also called the "gnomonic
   projection" and the "central projection".

2. The eta axis points due north in the adopted coordinate system.
   If the direction cosines represent observed (RA,Dec), the tangent
   plane coordinates (xi,eta) are conventionally called the
   "standard coordinates".  If the direction cosines are with
   respect to a right-handed triad, (xi,eta) are also right-handed.
   The units of (xi,eta) are, effectively, radians at the tangent
   point.

3. The method used is to complete the star vector in the (xi,eta)
   based triad and normalize it, then rotate the triad to put the
   tangent point at the pole with the x-axis aligned to zero
   longitude.  Writing (a0,b0) for the celestial spherical
   coordinates of the tangent point, the sequence of rotations is
   (b-pi/2) around the x-axis followed by (-a-pi/2) around the
   z-axis.

4. If vector v0 is not of unit length, the returned vector v will
   be wrong.

5. If vector v0 points at a pole, the returned vector v will be
   based on the arbitrary assumption that the longitude coordinate
   of the tangent point is zero.

6. This function is a member of the following set:

   | spherical       | vector            | solve for |
   |:----------------|:------------------|:----------|
   | [`tpxes`](@ref) |  [`tpxev`](@ref)  | xi,eta    |
   | [`tpsts`](@ref) | *[`tpstv`](@ref)* | star      |
   | [`tpors`](@ref) |  [`tporv`](@ref)  | origin    |

### References ###

- Calabretta M.R. & Greisen, E.W., 2002, "Representations of
  celestial coordinates in FITS", Astron.Astrophys. 395, 1077

- Green, R.M., "Spherical Astronomy", Cambridge University Press,
  1987, Chapter 13.
"""
function tpstv(xi, eta, vz)
    v = Array{Cdouble}(undef, 3)
    ccall((:eraTpstv, liberfa), Cvoid,
          (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}),
          xi, eta, vz, v)
    return v
end

"""
    tpxes(ra, dec, raz, decz)

In the tangent plane projection, given celestial spherical
coordinates for a star and the tangent point, solve for the star's
rectangular coordinates in the tangent plane.

### Given ###

- `a`, `b`: star's spherical coordinates
- `a0`, `b0`: tangent point's spherical coordinates

### Returned ###

- `status`
    - 0 = OK
    - 1 = star too far from axis
    - 2 = antistar on tangent plane
    - 3 = antistar too far from axis
- `xi`, `eta`: rectangular coordinates of star image (Note 2)

### Notes ###

1. The tangent plane projection is also called the "gnomonic
   projection" and the "central projection".

2. The eta axis points due north in the adopted coordinate system.
   If the spherical coordinates are observed (RA,Dec), the tangent
   plane coordinates (xi,eta) are conventionally called the
   "standard coordinates".  For right-handed spherical coordinates,
   (xi,eta) are also right-handed.  The units of (xi,eta) are,
   effectively, radians at the tangent point.

3. All angular arguments are in radians.

4. This function is a member of the following set:

   | spherical         | vector          | solve for |
   |:------------------|:----------------|:----------|
   | *[`tpxes`](@ref)* | [`tpxev`](@ref) | xi,eta    |
   |  [`tpsts`](@ref)  | [`tpstv`](@ref) | star      |
   |  [`tpors`](@ref)  | [`tporv`](@ref) | origin    |

### References ###

- Calabretta M.R. & Greisen, E.W., 2002, "Representations of
  celestial coordinates in FITS", Astron.Astrophys. 395, 1077

- Green, R.M., "Spherical Astronomy", Cambridge University Press,
  1987, Chapter 13.
"""
function tpxes(ra, dec, raz, decz)
    xi = Ref{Cdouble}()
    eta = Ref{Cdouble}()
    status = ccall((:eraTpxes, liberfa), Cint,
                   (Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
                   ra, dec, raz, decz, xi, eta)
    return status, xi[], eta[]
end

"""
    tr(r)

Transpose an r-matrix.

!!! warning "Deprecated"
    Use `r'` instead.

### Given ###

- `r`: R-matrix

### Returned ###

- `rt`: Transpose

### Note ###

   It is permissible for r and rt to be the same array.

### Called ###

- [`cr`](@ref): copy r-matrix

"""
tr

function _tr(r)
    @checkdims 3 3 r
    rt = zeros(Cdouble, 3, 3)
    ccall((:eraTr, liberfa), Cvoid,
          (Ptr{Cdouble}, Ptr{Cdouble}),
          permutedims(r), rt)
    return permutedims(rt)
end

@deprecate tr(r) r'

"""
    trxpv(r, pv)

Multiply a pv-vector by the transpose of an r-matrix.

!!! warning "Deprecated"
    Use `[r' * pv[1], r' * pv[2]]` instead.

### Given ###

- `r`: R-matrix
- `pv`: Pv-vector

### Returned ###

- `trpv`: R * pv

### Note ###

   It is permissible for pv and trpv to be the same array.

### Called ###

- [`tr`](@ref): transpose r-matrix
- [`rxpv`](@ref): product of r-matrix and pv-vector

"""
trxpv

function _trxpv(r, pv)
    @checkdims 3 3 r
    _pv = array_to_cmatrix(pv; n=3)
    rp = zeros(Cdouble, 3, 2)
    ccall((:eraTrxpv, liberfa), Cvoid,
            (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
            permutedims(r), _pv, rp)
    return cmatrix_to_array(rp)
end

@deprecate trxpv(r, pv) [r' * pv[1], r' * pv[2]]

"""
    trxp(r, p)

Multiply a p-vector by the transpose of an r-matrix.

!!! warning "Deprecated"
    Use `r' * p` instead.

### Given ###

- `r`: R-matrix
- `p`: P-vector

### Returned ###

- `trp`: R * p

### Note ###

   It is permissible for p and trp to be the same array.

### Called ###

- [`tr`](@ref): transpose r-matrix
- [`rxp`](@ref): product of r-matrix and p-vector

"""
trxp

function _trxp(r, p)
    @checkdims 3 3 r
    @checkdims 3 p
    rp = zeros(Cdouble, 3)
    ccall((:eraTrxp, liberfa), Cvoid,
            (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
            permutedims(r), p, rp)
    return rp
end

@deprecate trxp(r, p) r' * p

"""
    taiut1(tai1, tai2, dta)

Time scale transformation:  International Atomic Time, TAI, to
Universal Time, UT1.

### Given ###

- `tai1`, `tai2`: TAI as a 2-part Julian Date
- `dta`: UT1-TAI in seconds

### Returned ###

- `ut11`, `ut12`: UT1 as a 2-part Julian Date

### Notes ###

1. tai1+tai2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tai1 is the Julian
   Day Number and tai2 is the fraction of a day.  The returned
   UT11,UT12 follow suit.

2. The argument dta, i.e. UT1-TAI, is an observed quantity, and is
   available from IERS tabulations.

### Reference ###

- Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992)

"""
taiut1

"""
    tdbtt(tdb1, tdb2, dtr)

Time scale transformation:  Barycentric Dynamical Time, TDB, to
Terrestrial Time, TT.

### Given ###

- `tdb1`, `tdb2`: TDB as a 2-part Julian Date
- `dtr`: TDB-TT in seconds

### Returned ###

- `tt1`, `tt2`: TT as a 2-part Julian Date

### Notes ###

1. tdb1+tdb2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tdb1 is the Julian
   Day Number and tdb2 is the fraction of a day.  The returned
   tt1,tt2 follow suit.

2. The argument dtr represents the quasi-periodic component of the
   GR transformation between TT and TCB.  It is dependent upon the
   adopted solar-system ephemeris, and can be obtained by numerical
   integration, by interrogating a precomputed time ephemeris or by
   evaluating a model such as that implemented in the ERFA function
   [`dtdb`](@ref).   The quantity is dominated by an annual term of 1.7 ms
   amplitude.

3. TDB is essentially the same as Teph, the time argument for the
   JPL solar system ephemerides.

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- IAU 2006 Resolution 3

"""
tdbtt

"""
    tttdb(tt1, tt2, dtr)

Time scale transformation:  Terrestrial Time, TT, to Barycentric
Dynamical Time, TDB.

### Given ###

- `tt1`, `tt2`: TT as a 2-part Julian Date
- `dtr`: TDB-TT in seconds

### Returned ###

- `tdb1`, `tdb2`: TDB as a 2-part Julian Date

### Notes ###

1. tt1+tt2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where tt1 is the Julian Day Number
   and tt2 is the fraction of a day.  The returned tdb1,tdb2 follow
   suit.

2. The argument dtr represents the quasi-periodic component of the
   GR transformation between TT and TCB.  It is dependent upon the
   adopted solar-system ephemeris, and can be obtained by numerical
   integration, by interrogating a precomputed time ephemeris or by
   evaluating a model such as that implemented in the ERFA function
   [`dtdb`](@ref).   The quantity is dominated by an annual term of 1.7 ms
   amplitude.

3. TDB is essentially the same as Teph, the time argument for the JPL
   solar system ephemerides.

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- IAU 2006 Resolution 3

"""
tttdb

"""
    ttut1(tt1, tt2, dt)

Time scale transformation:  Terrestrial Time, TT, to Universal Time,
UT1.

### Given ###

- `tt1`, `tt2`: TT as a 2-part Julian Date
- `dt`: TT-UT1 in seconds

### Returned ###

- `ut11`, `ut12`: UT1 as a 2-part Julian Date

### Notes ###

1. tt1+tt2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where tt1 is the Julian Day Number
   and tt2 is the fraction of a day.  The returned ut11,ut12 follow
   suit.

2. The argument dt is classical Delta T.

### Reference ###

- Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992)

"""
ttut1

for name in ("taiut1",
             "tdbtt",
             "tttdb",
             "ttut1")
    f = Symbol(name)
    fc = "era" * uppercasefirst(name)
    @eval begin
        function ($f)(a, b, c)
            r1 = Ref{Cdouble}()
            r2 = Ref{Cdouble}()
            i = ccall(($fc, liberfa), Cint,
                      (Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
                      a, b, c, r1, r2)
            @assert i == 0
            r1[], r2[]
        end
    end
end

"""
    taitt(tai1, tai2)

Time scale transformation:  International Atomic Time, TAI, to
Terrestrial Time, TT.

### Given ###

- `tai1`, `tai2`: TAI as a 2-part Julian Date

### Returned ###

- `tt1`, `tt2`: TT as a 2-part Julian Date

### Note ###

   tai1+tai2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tai1 is the Julian
   Day Number and tai2 is the fraction of a day.  The returned
   tt1,tt2 follow suit.

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992)

"""
taitt

"""
    taiutc(tai1, tai2)

Time scale transformation:  International Atomic Time, TAI, to
Coordinated Universal Time, UTC.

### Given ###

- `tai1`, `tai2`: TAI as a 2-part Julian Date (Note 1)

### Returned ###

- `utc1`, `utc2`: UTC as a 2-part quasi Julian Date (Notes 1-3)

### Notes ###

1. tai1+tai2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tai1 is the Julian
   Day Number and tai2 is the fraction of a day.  The returned utc1
   and utc2 form an analogous pair, except that a special convention
   is used, to deal with the problem of leap seconds - see the next
   note.

2. JD cannot unambiguously represent UTC during a leap second unless
   special measures are taken.  The convention in the present
   function is that the JD day represents UTC days whether the
   length is 86399, 86400 or 86401 SI seconds.  In the 1960-1972 era
   there were smaller jumps (in either direction) each time the
   linear UTC(TAI) expression was changed, and these "mini-leaps"
   are also included in the ERFA convention.

3. The function [`d2dtf`](@ref) can be used to transform the UTC quasi-JD
   into calendar date and clock time, including UTC leap second
   handling.

4. The warning status "dubious year" flags UTCs that predate the
   introduction of the time scale or that are too far in the future
   to be trusted.  See [`dat`](@ref) for further details.

### Called ###

- [`utctai`](@ref): UTC to TAI

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992)

"""
taiutc

"""
    tcbtdb(tcb1, tcb2)

Time scale transformation:  Barycentric Coordinate Time, TCB, to
Barycentric Dynamical Time, TDB.

### Given ###

- `tcb1`, `tcb2`: TCB as a 2-part Julian Date

### Returned ###

- `tdb1`, `tdb2`: TDB as a 2-part Julian Date

### Notes ###

1. tcb1+tcb2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tcb1 is the Julian
   Day Number and tcb2 is the fraction of a day.  The returned
   tdb1,tdb2 follow suit.

2. The 2006 IAU General Assembly introduced a conventional linear
   transformation between TDB and TCB.  This transformation
   compensates for the drift between TCB and terrestrial time TT,
   and keeps TDB approximately centered on TT.  Because the
   relationship between TT and TCB depends on the adopted solar
   system ephemeris, the degree of alignment between TDB and TT over
   long intervals will vary according to which ephemeris is used.
   Former definitions of TDB attempted to avoid this problem by
   stipulating that TDB and TT should differ only by periodic
   effects.  This is a good description of the nature of the
   relationship but eluded precise mathematical formulation.  The
   conventional linear relationship adopted in 2006 sidestepped
   these difficulties whilst delivering a TDB that in practice was
   consistent with values before that date.

3. TDB is essentially the same as Teph, the time argument for the
   JPL solar system ephemerides.

### Reference ###

- IAU 2006 Resolution B3

"""
tcbtdb

"""
    tcgtt(tcg1, tcg2)

Time scale transformation:  Geocentric Coordinate Time, TCG, to
Terrestrial Time, TT.

### Given ###

- `tcg1`, `tcg2`: TCG as a 2-part Julian Date

### Returned ###

- `tt1`, `tt2`: TT as a 2-part Julian Date

### Note ###

   tcg1+tcg2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tcg1 is the Julian
   Day Number and tcg22 is the fraction of a day.  The returned
   tt1,tt2 follow suit.

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),.
    IERS Technical Note No. 32, BKG (2004)

- IAU 2000 Resolution B1.9

"""
tcgtt

"""
    tdbtcb(tdb1, tdb2)

Time scale transformation:  Barycentric Dynamical Time, TDB, to
Barycentric Coordinate Time, TCB.

### Given ###

- `tdb1`, `tdb2`: TDB as a 2-part Julian Date

### Returned ###

- `tcb1`, `tcb2`: TCB as a 2-part Julian Date

### Notes ###

1. tdb1+tdb2 is Julian Date, apportioned in any convenient way
   between the two arguments, for example where tdb1 is the Julian
   Day Number and tdb2 is the fraction of a day.  The returned
   tcb1,tcb2 follow suit.

2. The 2006 IAU General Assembly introduced a conventional linear
   transformation between TDB and TCB.  This transformation
   compensates for the drift between TCB and terrestrial time TT,
   and keeps TDB approximately centered on TT.  Because the
   relationship between TT and TCB depends on the adopted solar
   system ephemeris, the degree of alignment between TDB and TT over
   long intervals will vary according to which ephemeris is used.
   Former definitions of TDB attempted to avoid this problem by
   stipulating that TDB and TT should differ only by periodic
   effects.  This is a good description of the nature of the
   relationship but eluded precise mathematical formulation.  The
   conventional linear relationship adopted in 2006 sidestepped
   these difficulties whilst delivering a TDB that in practice was
   consistent with values before that date.

3. TDB is essentially the same as Teph, the time argument for the
   JPL solar system ephemerides.

### Reference ###

- IAU 2006 Resolution B3

"""
tdbtcb

"""
    tttai(tt1, tt2)

Time scale transformation:  Terrestrial Time, TT, to International
Atomic Time, TAI.

### Given ###

- `tt1`, `tt2`: TT as a 2-part Julian Date

### Returned ###

- `tai1`, `tai2`: TAI as a 2-part Julian Date

### Note ###

   tt1+tt2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where tt1 is the Julian Day Number
   and tt2 is the fraction of a day.  The returned tai1,tai2 follow
   suit.

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Explanatory Supplement to the Astronomical Almanac,
    P. Kenneth Seidelmann (ed), University Science Books (1992)

"""
tttai

"""
    tttcg(tt1, tt2)

Time scale transformation:  Terrestrial Time, TT, to Geocentric
Coordinate Time, TCG.

### Given ###

- `tt1`, `tt2`: TT as a 2-part Julian Date

### Returned ###

- `tcg1`, `tcg2`: TCG as a 2-part Julian Date

### Note ###

   tt1+tt2 is Julian Date, apportioned in any convenient way between
   the two arguments, for example where tt1 is the Julian Day Number
   and tt2 is the fraction of a day.  The returned tcg1,tcg2 follow
   suit.

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- IAU 2000 Resolution B1.9

"""
tttcg

for name in ("taitt",
             "taiutc",
             "tcbtdb",
             "tcgtt",
             "tdbtcb",
             "tttai",
             "tttcg")
    f = Symbol(name)
    fc = "era" * uppercasefirst(name)
    @eval begin
        function ($f)(a, b)
            r1 = Ref{Cdouble}()
            r2 = Ref{Cdouble}()
            i = ccall(($fc, liberfa), Cint,
                      (Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
                      a, b, r1, r2)
            @assert i == 0
            return r1[], r2[]
        end
    end
end

"""
    tf2a(s, ihour, imin, sec)

Convert hours, minutes, seconds to radians.

### Given ###

- `s`: Sign:  '-' = negative, otherwise positive
- `ihour`: Hours
- `imin`: Minutes
- `sec`: Seconds

### Returned ###

- `rad`: Angle in radians

### Notes ###

1.  The result is computed even if any of the range checks fail.

2.  Negative ihour, imin and/or sec produce a warning status, but
    the absolute value is used in the conversion.

3.  If there are multiple errors, the status value reflects only the
    first, the smallest taking precedence.

"""
tf2a

"""
    tf2d(s, ihour, imin, sec)

Convert hours, minutes, seconds to days.

### Given ###

- `s`: Sign:  '-' = negative, otherwise positive
- `ihour`: Hours
- `imin`: Minutes
- `sec`: Seconds

### Returned ###

- `days`: Interval in days

### Notes ###

1.  The result is computed even if any of the range checks fail.

2.  Negative ihour, imin and/or sec produce a warning status, but
    the absolute value is used in the conversion.

3.  If there are multiple errors, the status value reflects only the
    first, the smallest taking precedence.

"""
tf2d

for name in ("tf2a",
             "tf2d")
    f = Symbol(name)
    fc = "era" * uppercasefirst(name)
    @eval begin
        function ($f)(s, ideg, iamin, asec)
            rad = Ref{Cdouble}()
            i = ccall(($fc, liberfa), Cint,
                      (Cchar, Cint, Cint, Cdouble, Ref{Cdouble}),
                       s, ideg, iamin, asec, rad)
            if i == 1
                throw(ERFAException("ihour outside range 0-23"))
            elseif i == 2
                throw(ERFAException("imin outside range 0-59"))
            elseif i == 3
                throw(ERFAException("sec outside range 0-59.999..."))
            end
            return rad[]
        end
    end
end

