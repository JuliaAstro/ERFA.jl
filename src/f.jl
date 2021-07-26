"""
    fk425(r1950, d1950, dr1950, dd1950, p1950, v1950)

Convert B1950.0 FK4 star catalog data to J2000.0 FK5.

This function converts a star's catalog data from the old FK4
Bessel-Newcomb) system to the later IAU 1976 FK5 (Fricke) system.

### Given (all B1950.0, FK4) ###

- `r1950,d1950`: B1950.0 RA,Dec (rad)
- `dr1950,dd1950`: B1950.0 proper motions (rad/trop.yr)
- `p1950`: Parallax (arcsec)
- `v1950`: Radial velocity (km/s, +ve = moving away)

### Returned (all J2000.0, FK5) ###

- `r2000,d2000`: J2000.0 RA,Dec (rad)
- `dr2000,dd2000`: J2000.0 proper motions (rad/Jul.yr)
- `p2000`: Parallax (arcsec)
- `v2000`: Radial velocity (km/s, +ve = moving away)

### Notes ###

1. The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
   and are per year rather than per century.

2. The conversion is somewhat complicated, for several reasons:

   - Change of standard epoch from B1950.0 to J2000.0.

   - An intermediate transition date of 1984 January 1.0 TT.

   - A change of precession model.

   - Change of time unit for proper motion (tropical to Julian).

   - FK4 positions include the E-terms of aberration, to simplify
     the hand computation of annual aberration.  FK5 positions
     assume a rigorous aberration computation based on the Earth's
     barycentric velocity.

   - The E-terms also affect proper motions, and in particular cause
     objects at large distances to exhibit fictitious proper
     motions.

   The algorithm is based on Smith et al. (1989) and Yallop et al.
   (1989), which presented a matrix method due to Standish (1982) as
   developed by Aoki et al. (1983), using Kinoshita's development of
   Andoyer's post-Newcomb precession.  The numerical constants from
   Seidelmann (1992) are used canonically.

3. Conversion from B1950.0 FK4 to J2000.0 FK5 only is provided for.
   Conversions for different epochs and equinoxes would require
   additional treatment for precession, proper motion and E-terms.

4. In the FK4 catalog the proper motions of stars within 10 degrees
   of the poles do not embody differential E-terms effects and
   should, strictly speaking, be handled in a different manner from
   stars outside these regions.  However, given the general lack of
   homogeneity of the star data available for routine astrometry,
   the difficulties of handling positions that may have been
   determined from astrometric fields spanning the polar and non-
   polar regions, the likelihood that the differential E-terms
   effect was not taken into account when allowing for proper motion
   in past astrometry, and the undesirability of a discontinuity in
   the algorithm, the decision has been made in this ERFA algorithm
   to include the effects of differential E-terms on the proper
   motions for all stars, whether polar or not.  At epoch J2000.0,
   and measuring "on the sky" rather than in terms of RA change, the
   errors resulting from this simplification are less than
   1 milliarcsecond in position and 1 milliarcsecond per century in
   proper motion.

### Called ###

- [`anp`](@ref): normalize angle into range 0 to 2pi
- [`pv2s`](@ref): pv-vector to spherical coordinates
- [`pdp`](@ref): scalar product of two p-vectors
- [`pvmpv`](@ref): pv-vector minus pv_vector
- [`pvppv`](@ref): pv-vector plus pv_vector
- [`s2pv`](@ref): spherical coordinates to pv-vector
- [`sxp`](@ref): multiply p-vector by scalar

### References ###

- Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0
  FK4-based positions of stars to epoch J2000.0 positions in
  accordance with the new IAU resolutions".  Astron.Astrophys.
  128, 263-267.

- Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
  Astronomical Almanac", ISBN 0-935702-68-7.

- Smith, C.A. et al., 1989, "The transformation of astrometric
  catalog systems to the equinox J2000.0".  Astron.J. 97, 265.

- Standish, E.M., 1982, "Conversion of positions and proper motions
  from B1950.0 to the IAU system at J2000.0".  Astron.Astrophys.,
  115, 1, 20-22.

- Yallop, B.D. et al., 1989, "Transformation of mean star places
  from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space".
  Astron.J. 97, 274.
"""
function fk425(r1950, d1950, dr1950, dd1950, p1950, v1950)
    r2000 = Ref{Cdouble}()
    d2000 = Ref{Cdouble}()
    dr2000 = Ref{Cdouble}()
    dd2000 = Ref{Cdouble}()
    p2000 = Ref{Cdouble}()
    v2000 = Ref{Cdouble}()
    ccall((:eraFk425, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
           Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble},
           Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble},),
          r1950, d1950, dr1950, dd1950, p1950, v1950,
          r2000, d2000, dr2000, dd2000, p2000, v2000)
    return r2000[], d2000[], dr2000[], dd2000[], p2000[], v2000[]
end

"""
    fk45z(r1950, d1950, bepoch)

Convert a B1950.0 FK4 star position to J2000.0 FK5, assuming zero
proper motion in the FK5 system.

This function converts a star's catalog data from the old FK4
(Bessel-Newcomb) system to the later IAU 1976 FK5 (Fricke) system,
in such a way that the FK5 proper motion is zero.  Because such a
star has, in general, a non-zero proper motion in the FK4 system,
the function requires the epoch at which the position in the FK4
system was determined.

### Given ###

- `r1950`, `d1950`: B1950.0 FK4 RA,Dec at epoch (rad)
- `bepoch`: Besselian epoch (e.g. 1979.3)

### Returned ###

- `r2000`, `d2000`: J2000.0 FK5 RA,Dec (rad)

### Notes ###

1. The epoch bepoch is strictly speaking Beselian, but if a
   Julian epoch is supplied the result will be affected only to a
   negligible extent.

2. The method is from Appendix 2 of Aoki et al. (1983), but using
   the constants of Seidelmann (1992).  See the function eraFk425
   for a general introduction to the FK4 to FK5 conversion.

3. Conversion from equinox B1950.0 FK4 to equinox J2000.0 FK5 only
   is provided for.  Conversions for different starting and/or
   ending epochs would require additional treatment for precession,
   proper motion and E-terms.

4. In the FK4 catalog the proper motions of stars within 10 degrees
   of the poles do not embody differential E-terms effects and
   should, strictly speaking, be handled in a different manner from
   stars outside these regions.  However, given the general lack of
   homogeneity of the star data available for routine astrometry,
   the difficulties of handling positions that may have been
   determined from astrometric fields spanning the polar and non-
   polar regions, the likelihood that the differential E-terms
   effect was not taken into account when allowing for proper motion
   in past astrometry, and the undesirability of a discontinuity in
   the algorithm, the decision has been made in this ERFA algorithm
   to include the effects of differential E-terms on the proper
   motions for all stars, whether polar or not.  At epoch 2000.0,
   and measuring "on the sky" rather than in terms of RA change, the
   errors resulting from this simplification are less than
   1 milliarcsecond in position and 1 milliarcsecond per century in
   proper motion.

### References ###

- Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0
  FK4-based positions of stars to epoch J2000.0 positions in
  accordance with the new IAU resolutions".  Astron.Astrophys.
  128, 263-267.

- Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
  Astronomical Almanac", ISBN 0-935702-68-7.

### Called ###

- [`anp`](@ref): normalize angle into range 0 to 2pi
- [`c2s`](@ref): p-vector to spherical
- [`epb2jd`](@ref): Besselian epoch to Julian date
- [`epj`](@ref): Julian date to Julian epoch
- [`pdp`](@ref): scalar product of two p-vectors
- [`pmp`](@ref): p-vector minus p-vector
- [`ppsp`](@ref): p-vector plus scaled p-vector
- [`pvu`](@ref): update a pv-vector
- [`s2c`](@ref): spherical to p-vectors
"""
function fk45z(r1950, d1950, bepoch)
    r2000 = Ref{Cdouble}()
    d2000 = Ref{Cdouble}()
    ccall((:eraFk45z, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
          r1950, d1950, bepoch, r2000, d2000)
    return r2000[], d2000[]
end

"""
    function fk524(r2000, d2000, dr2000, dd2000, p2000, v2000)

Convert J2000.0 FK5 star catalog data to B1950.0 FK4.

### Given (all J2000.0, FK5) ###

- `r2000`, `d2000`: J2000.0 RA,Dec (rad)
- `dr2000`, `dd2000`: J2000.0 proper motions (rad/Jul.yr)
- `p2000`: parallax (arcsec)
- `v2000`: radial velocity (km/s, +ve = moving away)

### Returned (all B1950.0, FK4) ###

- `r1950`, `d1950`: B1950.0 RA,Dec (rad)
- `dr1950`, `dd1950`: B1950.0 proper motions (rad/trop.yr)
- `p1950`: parallax (arcsec)
- `v1950`: radial velocity (km/s, +ve = moving away)

### Notes ###

1. The proper motions in RA are dRA/dt rather than cos(Dec)*dRA/dt,
   and are per year rather than per century.

2. The conversion is somewhat complicated, for several reasons:

   - Change of standard epoch from J2000.0 to B1950.0.

   - An intermediate transition date of 1984 January 1.0 TT.

   - A change of precession model.

   - Change of time unit for proper motion (Julian to tropical).

   - FK4 positions include the E-terms of aberration, to simplify
     the hand computation of annual aberration.  FK5 positions
     assume a rigorous aberration computation based on the Earth's
     barycentric velocity.

   - The E-terms also affect proper motions, and in particular cause
     objects at large distances to exhibit fictitious proper
     motions.

   The algorithm is based on Smith et al. (1989) and Yallop et al.
   (1989), which presented a matrix method due to Standish (1982) as
   developed by Aoki et al. (1983), using Kinoshita's development of
   Andoyer's post-Newcomb precession.  The numerical constants from
   Seidelmann (1992) are used canonically.

4. In the FK4 catalog the proper motions of stars within 10 degrees
   of the poles do not embody differential E-terms effects and
   should, strictly speaking, be handled in a different manner from
   stars outside these regions.  However, given the general lack of
   homogeneity of the star data available for routine astrometry,
   the difficulties of handling positions that may have been
   determined from astrometric fields spanning the polar and non-
   polar regions, the likelihood that the differential E-terms
   effect was not taken into account when allowing for proper motion
   in past astrometry, and the undesirability of a discontinuity in
   the algorithm, the decision has been made in this ERFA algorithm
   to include the effects of differential E-terms on the proper
   motions for all stars, whether polar or not.  At epoch J2000.0,
   and measuring "on the sky" rather than in terms of RA change, the
   errors resulting from this simplification are less than
   1 milliarcsecond in position and 1 milliarcsecond per century in
   proper motion.

### Called ###

- [`anp`](@ref): normalize angle into range 0 to 2pi
- [`pdp`](@ref): scalar product of two p-vectors
- [`pm`](@ref): modulus of p-vector
- [`pmp`](@ref): p-vector minus p-vector
- [`ppp`](@ref): p-vector plus p-vector
- [`pv2s`](@ref): pv-vector to spherical coordinates
- [`s2pv`](@ref): spherical coordinates to pv-vector
- [`sxp`](@ref): multiply p-vector by scalar

### References ###

- Aoki, S. et al., 1983, "Conversion matrix of epoch B1950.0
  FK4-based positions of stars to epoch J2000.0 positions in
  accordance with the new IAU resolutions".  Astron.Astrophys.
  128, 263-267.

- Seidelmann, P.K. (ed), 1992, "Explanatory Supplement to the
  Astronomical Almanac", ISBN 0-935702-68-7.

- Smith, C.A. et al., 1989, "The transformation of astrometric
  catalog systems to the equinox J2000.0".  Astron.J. 97, 265.

- Standish, E.M., 1982, "Conversion of positions and proper motions
  from B1950.0 to the IAU system at J2000.0".  Astron.Astrophys.,
  115, 1, 20-22.

- Yallop, B.D. et al., 1989, "Transformation of mean star places
  from FK4 B1950.0 to FK5 J2000.0 using matrices in 6-space".
  Astron.J. 97, 274.
"""
function fk524(r2000, d2000, dr2000, dd2000, p2000, v2000)
    r1950 = Ref{Cdouble}()
    d1950 = Ref{Cdouble}()
    dr1950 = Ref{Cdouble}()
    dd1950 = Ref{Cdouble}()
    p1950 = Ref{Cdouble}()
    v1950 = Ref{Cdouble}()
    ccall((:eraFk524, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
           Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble},
           Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
          r2000, d2000, dr2000, dd2000, p2000, v2000,
          r1950, d1950, dr1950, dd1950, p1950, v1950)
    return r1950[], d1950[], dr1950[], dd1950[], p1950[], v1950[]
end

"""
    fk54z(r2000, d2000, bepoch)

Convert a J2000.0 FK5 star position to B1950.0 FK4, assuming zero
proper motion in FK5 and parallax.

### Given ###

- `r2000`, `d2000`: J2000.0 FK5 RA,Dec (rad)
- `bepoch`: Besselian epoch (e.g. 1950.0)

### Returned ###

- `r1950`, `d1950`: B1950.0 FK4 RA,Dec (rad) at epoch BEPOCH
- `dr1950`, `dd1950`: B1950.0 FK4 proper motions (rad/trop.yr)

### Notes ###

1. In contrast to the eraFk524 function, here the FK5 proper
   motions, the parallax and the radial velocity are presumed zero.

2. This function converts a star position from the IAU 1976 FK5
   (Fricke) system to the former FK4 (Bessel-Newcomb) system, for
   cases such as distant radio sources where it is presumed there is
   zero parallax and no proper motion.  Because of the E-terms of
   aberration, such objects have (in general) non-zero proper motion
   in FK4, and the present function returns those fictitious proper
   motions.

3. Conversion from B1950.0 FK4 to J2000.0 FK5 only is provided for.
   Conversions involving other equinoxes would require additional
   treatment for precession.

4. The position returned by this function is in the B1950.0 FK4
   reference system but at Besselian epoch BEPOCH.  For comparison
   with catalogs the BEPOCH argument will frequently be 1950.0. (In
   this context the distinction between Besselian and Julian epoch
   is insignificant.)

5. The RA component of the returned (fictitious) proper motion is
   dRA/dt rather than cos(Dec)*dRA/dt.

### Called ###

- [`anp`](@ref): normalize angle into range 0 to 2pi
- [`c2s`](@ref): p-vector to spherical
- [`fk524`](@ref): FK4 to FK5
- [`s2c`](@ref): spherical to p-vector
"""
function fk54z(r2000, d2000, bepoch)
    r1950 = Ref{Cdouble}()
    d1950 = Ref{Cdouble}()
    dr1950 = Ref{Cdouble}()
    dd1950 = Ref{Cdouble}()
    ccall((:eraFk54z, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble,
           Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
          r2000, d2000, bepoch, r1950, d1950, dr1950, dd1950)
    return r1950[], d1950[], dr1950[], dd1950[]
end

"""
    fk5hip()

FK5 to Hipparcos rotation and spin.

### Returned ###

- `r5h`: R-matrix: FK5 rotation wrt Hipparcos (Note 2)
- `s5h`: R-vector: FK5 spin wrt Hipparcos (Note 3)

### Notes ###

1. This function models the FK5 to Hipparcos transformation as a
   pure rotation and spin;  zonal errors in the FK5 catalogue are
   not taken into account.

2. The r-matrix r5h operates in the sense:

         P_Hipparcos = r5h x P_FK5

   where P_FK5 is a p-vector in the FK5 frame, and P_Hipparcos is
   the equivalent Hipparcos p-vector.

3. The r-vector s5h represents the time derivative of the FK5 to
   Hipparcos rotation.  The units are radians per year (Julian,
   TDB).

### Called ###

- [`rv2m`](@ref): r-vector to r-matrix

### Reference ###

- F.Mignard & M.Froeschle, Astron. Astrophys. 354, 732-739 (2000).

"""
function fk5hip()
    r5h = zeros(Cdouble, 3, 3)
    s5h = zeros(Cdouble, 3)
    ccall((:eraFk5hip, liberfa), Cvoid,
          (Ptr{Cdouble}, Ptr{Cdouble}),
          r5h, s5h)
    return r5h, s5h
end

"""
    fk5hz(r5, d5, date1, date2)

Transform an FK5 (J2000.0) star position into the system of the
Hipparcos catalogue, assuming zero Hipparcos proper motion.

### Given ###

- `r5`: FK5 RA (radians), equinox J2000.0, at date
- `d5`: FK5 Dec (radians), equinox J2000.0, at date
- `date1`, `date2`: TDB date (Notes 1,2)

### Returned ###

- `rh`: Hipparcos RA (radians)
- `dh`: Hipparcos Dec (radians)

### Notes ###

1. This function converts a star position from the FK5 system to
   the Hipparcos system, in such a way that the Hipparcos proper
   motion is zero.  Because such a star has, in general, a non-zero
   proper motion in the FK5 system, the function requires the date
   at which the position in the FK5 system was determined.

2. The TT date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TT)=2450123.7 could be expressed in any of these ways,
   among others:

   | `date1`   |     `date2` | Method      |
   |:----------|:------------|:------------|
   | 2450123.7 |         0.0 | JD          |
   | 2451545.0 |     -1421.3 | J2000       |
   | 2400000.5 |     50123.2 | MJD         |
   | 2450123.5 |         0.2 | date & time |

   The JD method is the most natural and convenient to use in
   cases where the loss of several decimal digits of resolution
   is acceptable.  The J2000 method is best matched to the way
   the argument is handled internally and will deliver the
   optimum resolution.  The MJD method and the date & time methods
   are both good compromises between resolution and convenience.

3. The FK5 to Hipparcos transformation is modeled as a pure
   rotation and spin;  zonal errors in the FK5 catalogue are not
   taken into account.

4. The position returned by this function is in the Hipparcos
   reference system but at date date1+date2.

5. See also [`fk52h`](@ref), [`h2fk5`](@ref), [`hfk5z`](@ref).

### Called ###

- [`s2c`](@ref): spherical coordinates to unit vector
- [`fk5hip`](@ref): FK5 to Hipparcos rotation and spin
- [`sxp`](@ref): multiply p-vector by scalar
- [`rv2m`](@ref): r-vector to r-matrix
- [`trxp`](@ref): product of transpose of r-matrix and p-vector
- [`pxp`](@ref): vector product of two p-vectors
- [`c2s`](@ref): p-vector to spherical
- [`anp`](@ref): normalize angle into range 0 to 2pi

### Reference ###

- F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.

"""
function fk5hz(r5, d5, date1, date2)
    rh = Ref{Cdouble}()
    dh = Ref{Cdouble}()
    ccall((:eraFk5hz, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
          r5, d5, date1, date2, rh, dh)
    return rh[], dh[]
end

"""
    fw2xy(gamb, phib, psi, eps)

CIP X,Y given Fukushima-Williams bias-precession-nutation angles.

### Given ###

- `gamb`: F-W angle gamma_bar (radians)
- `phib`: F-W angle phi_bar (radians)
- `psi`: F-W angle psi (radians)
- `eps`: F-W angle epsilon (radians)

### Returned ###

- `x`, `y`: CIP unit vector X,Y

### Notes ###

1. Naming the following points:

   - e = J2000.0 ecliptic pole,
   - p = GCRS pole
   - E = ecliptic pole of date,
   - and P = CIP,

   the four Fukushima-Williams angles are as follows:

   - gamb = gamma = epE
   - phib = phi = pE
   - psi = psi = pEP
   - eps = epsilon = EP

2. The matrix representing the combined effects of frame bias,
   precession and nutation is:

   `NxPxB = R_1(-epsA).R_3(-psi).R_1(phib).R_3(gamb)`

   The returned values x,y are elements `[3, 1]` and `[2, 1]` of the
   matrix.  Near J2000.0, they are essentially angles in radians.

### Called ###

- [`fw2m`](@ref): F-W angles to r-matrix
- [`bpn2xy`](@ref): extract CIP X,Y coordinates from NPB matrix

### Reference ###

- Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351

"""
function fw2xy(gamb, phib, psi, eps)
    x = Ref{Cdouble}()
    y = Ref{Cdouble}()
    ccall((:eraFw2xy, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
          gamb, phib, psi, eps, x, y)
    return x[], y[]
end

"""
    fad03(t)

Fundamental argument, IERS Conventions (2003):

mean elongation of the Moon from the Sun.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- `D`, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   is from Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

"""
fad03

"""
    fae03(t)

Fundamental argument, IERS Conventions (2003): Mean longitude of Earth.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- Mean longitude of Earth, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

- Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111

"""
fae03

"""
    faf03(t)

Fundamental argument, IERS Conventions (2003): Mean longitude of the Moon minus
mean longitude of the ascending node.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- `F`, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   is from Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

"""
faf03

"""
    faju03(t)

Fundamental argument, IERS Conventions (2003): Mean longitude of Jupiter.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- Mean longitude of Jupiter, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

- Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111

"""
faju03

"""
    fal03(t)

Fundamental argument, IERS Conventions (2003):

mean anomaly of the Moon.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- `l`, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   is from Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

"""
fal03

"""
    falp03(t)

Fundamental argument, IERS Conventions (2003):

mean anomaly of the Sun.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- `l'`, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   is from Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

"""
falp03

"""
    fama03(t)

Fundamental argument, IERS Conventions (2003): Mean longitude of Mars.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- Mean longitude of Mars, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

- Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111

"""
fama03

"""
    fame03(t)

Fundamental argument, IERS Conventions (2003): Mean longitude of Mercury.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- Mean longitude of Mercury, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

- Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111

"""
fame03

"""
    fane03(t)

Fundamental argument, IERS Conventions (2003): Mean longitude of Neptune.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- Mean longitude of Neptune, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   is adapted from Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

"""
fane03

"""
    faom03(t)

Fundamental argument, IERS Conventions (2003): Mean longitude of the Moon's ascending node.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- `Omega`, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   is from Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

"""
faom03

"""
    fapa03(t)

Fundamental argument, IERS Conventions (2003):

general accumulated precession in longitude.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- General precession in longitude, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003).  It
   is taken from Kinoshita & Souchay (1990) and comes originally
   from Lieske et al. (1977).

### References ###

- Kinoshita, H. and Souchay J. 1990, Celest.Mech. and Dyn.Astron.
    48, 187

- Lieske, J.H., Lederle, T., Fricke, W. & Morando, B. 1977,
    Astron.Astrophys. 58, 1-16

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

"""
fapa03

"""
    fasa03(t)

Fundamental argument, IERS Conventions (2003): Mean longitude of Saturn.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- Mean longitude of Saturn, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

- Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111

"""
fasa03

"""
    faur03(t)

Fundamental argument, IERS Conventions (2003): Mean longitude of Uranus.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned  ###

- Mean longitude of Uranus, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   is adapted from Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

"""
faur03

"""
    fave03(t)

Fundamental argument, IERS Conventions (2003): Mean longitude of Venus.

### Given ###

- `t`: TDB, Julian centuries since J2000.0 (Note 1)

### Returned ###

- Mean longitude of Venus, radians (Note 2)

### Notes ###

1. Though t is strictly TDB, it is usually more convenient to use
   TT, which makes no significant difference.

2. The expression used is as adopted in IERS Conventions (2003) and
   comes from Souchay et al. (1999) after Simon et al. (1994).

### References ###

- McCarthy, D. D., Petit, G. (eds.), IERS Conventions (2003),
    IERS Technical Note No. 32, BKG (2004)

- Simon, J.-L., Bretagnon, P., Chapront, J., Chapront-Touze, M.,
    Francou, G., Laskar, J. 1994, Astron.Astrophys. 282, 663-683

- Souchay, J., Loysel, B., Kinoshita, H., Folgueira, M. 1999,
    Astron.Astrophys.Supp.Ser. 135, 111

"""
fave03

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
    fc = "era" * uppercasefirst(name)
    @eval ($f)(d) = ccall(($fc, liberfa), Cdouble, (Cdouble,), d)
end

"""
    fk52h(ra, dec, dra, ddec, px, rv)

Transform FK5 (J2000.0) star data into the Hipparcos system.

### Given (all FK5, equinox J2000.0, epoch J2000.0) ###

- `r5`: RA (radians)
- `d5`: Dec (radians)
- `dr5`: Proper motion in RA (dRA/dt, rad/Jyear)
- `dd5`: Proper motion in Dec (dDec/dt, rad/Jyear)
- `px5`: Parallax (arcsec)
- `rv5`: Radial velocity (km/s, positive = receding)

### Returned (all Hipparcos, epoch J2000.0) ###

- `rh`: RA (radians)
- `dh`: Dec (radians)
- `drh`: proper motion in RA (dRA/dt, rad/Jyear)
- `ddh`: proper motion in Dec (dDec/dt, rad/Jyear)
- `pxh`: parallax (arcsec)
- `rvh`: radial velocity (km/s, positive = receding)

### Notes ###

1. This function transforms FK5 star positions and proper motions
   into the system of the Hipparcos catalog.

2. The proper motions in RA are dRA/dt rather than
   cos(Dec)*dRA/dt, and are per year rather than per century.

3. The FK5 to Hipparcos transformation is modeled as a pure
   rotation and spin;  zonal errors in the FK5 catalog are not
   taken into account.

4. See also [`h2fk5`](@ref), [`fk5hz`](@ref), [`hfk5z`](@ref).

### Called ###

- [`starpv`](@ref): star catalog data to space motion pv-vector
- [`fk5hip`](@ref): FK5 to Hipparcos rotation and spin
- [`rxp`](@ref): product of r-matrix and p-vector
- [`pxp`](@ref): vector product of two p-vectors
- [`ppp`](@ref): p-vector plus p-vector
- [`pvstar`](@ref): space motion pv-vector to star catalog data

### Reference ###

- F.Mignard & M.Froeschle, Astron. Astrophys. 354, 732-739 (2000).

"""
function fk52h(ra, dec, dra, ddec, px, rv)
    r = Ref{Cdouble}()
    d = Ref{Cdouble}()
    dr = Ref{Cdouble}()
    dd = Ref{Cdouble}()
    p = Ref{Cdouble}()
    v = Ref{Cdouble}()
    ccall((:eraFk52h, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble},
           Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
          ra, dec, dra, ddec, px, rv, r, d, dr, dd, p, v)
    return r[], d[], dr[], dd[], p[], v[]
end

"""
    fw2m(x, y, s, t)

Form rotation matrix given the Fukushima-Williams angles.

### Given ###

- `gamb`: F-W angle gamma_bar (radians)
- `phib`: F-W angle phi_bar (radians)
- `psi`: F-W angle psi (radians)
- `eps`: F-W angle epsilon (radians)

### Returned ###

- `r`: Rotation matrix

### Notes ###

1. Naming the following points:

   - e = J2000.0 ecliptic pole,
   - p = GCRS pole,
   - E = ecliptic pole of date,
   - and P = CIP,

   the four Fukushima-Williams angles are as follows:

    - gamb = gamma = epE
    - phib = phi = pE
    - psi = psi = pEP
    - eps = epsilon = EP

2. The matrix representing the combined effects of frame bias,
   precession and nutation is:

   `NxPxB = R_1(-eps).R_3(-psi).R_1(phib).R_3(gamb)`

3. Three different matrices can be constructed, depending on the
   supplied angles:

   -  To obtain the nutation x precession x frame bias matrix,
      generate the four precession angles, generate the nutation
      components and add them to the `psi_bar` and `epsilon_A` angles,
      and call the present function.

   -  To obtain the precession x frame bias matrix, generate the
      four precession angles and call the present function.

   -  To obtain the frame bias matrix, generate the four precession
      angles for date J2000.0 and call the present function.

   The nutation-only and precession-only matrices can if necessary
   be obtained by combining these three appropriately.

### Called ###

- [`ir`](@ref): initialize r-matrix to identity
- [`rz`](@ref): rotate around Z-axis
- [`rx`](@ref): rotate around X-axis

### Reference ###

- Hilton, J. et al., 2006, Celest.Mech.Dyn.Astron. 94, 351

"""
function fw2m(x, y, s, t)
    r = zeros(Cdouble, 3, 3)
    ccall((:eraFw2m, liberfa), Cvoid,
            (Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Cdouble}),
            x, y, s, t, r)
    return permutedims(r)
end
