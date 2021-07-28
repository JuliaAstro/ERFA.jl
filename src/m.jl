"""
    moon98(date1, date2)

Approximate geocentric position and velocity of the Moon.

n.b. Not IAU-endorsed and without canonical status.

### Given ###

- `date1`: TT date part A (Notes 1,4)
- `date2`: TT date part B (Notes 1,4)

### Returned ###

- `pv`: Moon p,v, GCRS (AU, AU/d, Note 5)

### Notes ###

1. The TT date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TT)=2450123.7 could be expressed in any of these ways, among
   others:

   | `date1`   |     `date2` | Method      |
   |:----------|:------------|:------------|
   | 2450123.7 |         0.0 | JD          |
   | 2451545.0 |     -1421.3 | J2000       |
   | 2400000.5 |     50123.2 | MJD         |
   | 2450123.5 |         0.2 | date & time |

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 method is best matched to the way the
   argument is handled internally and will deliver the optimum
   resolution.  The MJD method and the date & time methods are both
   good compromises between resolution and convenience.  The limited
   accuracy of the present algorithm is such that any of the methods
   is satisfactory.

2. This function is a full implementation of the algorithm
   published by Meeus (see reference) except that the light-time
   correction to the Moon's mean longitude has been omitted.

3. Comparisons with ELP/MPP02 over the interval 1950-2100 gave RMS
   errors of 2.9 arcsec in geocentric direction, 6.1 km in position
   and 36 mm/s in velocity.  The worst case errors were 18.3 arcsec
   in geocentric direction, 31.7 km in position and 172 mm/s in
   velocity.

4. The original algorithm is expressed in terms of "dynamical time",
   which can either be TDB or TT without any significant change in
   accuracy.  UT cannot be used without incurring significant errors
   (30 arcsec in the present era) due to the Moon's 0.5 arcsec/sec
   movement.

5. The result is with respect to the GCRS (the same as J2000.0 mean
   equator and equinox to within 23 mas).

6. Velocity is obtained by a complete analytical differentiation
   of the Meeus model.

7. The Meeus algorithm generates position and velocity in mean
   ecliptic coordinates of date, which the present function then
   rotates into GCRS.  Because the ecliptic system is precessing,
   there is a coupling between this spin (about 1.4 degrees per
   century) and the Moon position that produces a small velocity
   contribution.  In the present function this effect is neglected
   as it corresponds to a maximum difference of less than 3 mm/s and
   increases the RMS error by only 0.4%.

### References ###

- Meeus, J., Astronomical Algorithms, 2nd edition, Willmann-Bell,
  1998, p337.

- Simon, J.L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
  Francou, G. & Laskar, J., Astron.Astrophys., 1994, 282, 663

### Called ###

- [`s2pv`](@ref): spherical coordinates to pv-vector
- [`pfw06`](@ref): bias-precession F-W angles, IAU 2006
- [`ir`](@ref): initialize r-matrix to identity
- [`rz`](@ref): rotate around Z-axis
- [`rx`](@ref): rotate around X-axis
- [`rxpv`](@ref): product of r-matrix and pv-vector
"""
function moon98(date1, date2)
    pv = Array{Cdouble}(undef, 3, 2)
    ccall((:eraMoon98, liberfa), Cvoid,
          (Cdouble, Cdouble, Ptr{Cdouble}),
          date1, date2, pv)
    return cmatrix_to_array(pv)
end

