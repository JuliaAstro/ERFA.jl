"""
    hfk5z(rh, dh, date1, date2)

Transform a Hipparcos star position into FK5 J2000.0, assuming
zero Hipparcos proper motion.

### Given ###

- `rh`: Hipparcos RA (radians)
- `dh`: Hipparcos Dec (radians)
- `date1`, `date2`: TDB date (Note 1)

### Returned (all FK5, equinox J2000.0, date date1+date2) ###

- `r5`: RA (radians)
- `d5`: Dec (radians)
- `dr5`: FK5 RA proper motion (rad/year, Note 4)
- `dd5`: Dec proper motion (rad/year, Note 4)

### Notes ###

1. The TT date date1+date2 is a Julian Date, apportioned in any
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

2. The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

3. The FK5 to Hipparcos transformation is modeled as a pure rotation
   and spin;  zonal errors in the FK5 catalogue are not taken into
   account.

4. It was the intention that Hipparcos should be a close
   approximation to an inertial frame, so that distant objects have
   zero proper motion;  such objects have (in general) non-zero
   proper motion in FK5, and this function returns those fictitious
   proper motions.

5. The position returned by this function is in the FK5 J2000.0
   reference system but at date date1+date2.

6. See also [`fk52h`](@ref), [`h2fk5`](@ref), [`fk5hz`](@ref).

### Called ###

- [`s2c`](@ref): spherical coordinates to unit vector
- [`fk5hip`](@ref): FK5 to Hipparcos rotation and spin
- [`rxp`](@ref): product of r-matrix and p-vector
- [`sxp`](@ref): multiply p-vector by scalar
- [`rxr`](@ref): product of two r-matrices
- [`trxp`](@ref): product of transpose of r-matrix and p-vector
- [`pxp`](@ref): vector product of two p-vectors
- [`pv2s`](@ref): pv-vector to spherical
- [`anp`](@ref): normalize angle into range 0 to 2pi

### Reference ###

- F.Mignard & M.Froeschle, 2000, Astron.Astrophys. 354, 732-739.

"""
function hfk5z(rh, dh, date1, date2)
    r5 = Ref{Cdouble}()
    d5 = Ref{Cdouble}()
    dr5 = Ref{Cdouble}()
    dd5 = Ref{Cdouble}()
    ccall((:eraHfk5z, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
          rh, dh, date1, date2, r5, d5, dr5, dd5)
    return r5[], d5[], dr5[], dd5[]
end

"""
    h2fk5(ra, dec, dra, ddec, px, rv)

Transform Hipparcos star data into the FK5 (J2000.0) system.

### Given (all Hipparcos, epoch J2000.0) ###

- `rh`: RA (radians)
- `dh`: Dec (radians)
- `drh`: Proper motion in RA (dRA/dt, rad/Jyear)
- `ddh`: Proper motion in Dec (dDec/dt, rad/Jyear)
- `pxh`: Parallax (arcsec)
- `rvh`: Radial velocity (km/s, positive = receding)

### Returned (all FK5, equinox J2000.0, epoch J2000.0) ###

- `r5`: RA (radians)
- `d5`: Dec (radians)
- `dr5`: Proper motion in RA (dRA/dt, rad/Jyear)
- `dd5`: Proper motion in Dec (dDec/dt, rad/Jyear)
- `px5`: Parallax (arcsec)
- `rv5`: Radial velocity (km/s, positive = receding)

### Notes ###

1. This function transforms Hipparcos star positions and proper
   motions into FK5 J2000.0.

2. The proper motions in RA are dRA/dt rather than
   cos(Dec)*dRA/dt, and are per year rather than per century.

3. The FK5 to Hipparcos transformation is modeled as a pure
   rotation and spin;  zonal errors in the FK5 catalog are not
   taken into account.

4. See also [`fk52h`](@ref), [`fk5hz`](@ref), [`hfk5z`](@ref).

### Called ###

- [`starpv`](@ref): star catalog data to space motion pv-vector
- [`fk5hip`](@ref): FK5 to Hipparcos rotation and spin
- [`rv2m`](@ref): r-vector to r-matrix
- [`rxp`](@ref): product of r-matrix and p-vector
- [`trxp`](@ref): product of transpose of r-matrix and p-vector
- [`pxp`](@ref): vector product of two p-vectors
- [`pmp`](@ref): p-vector minus p-vector
- [`pvstar`](@ref): space motion pv-vector to star catalog data

### Reference ###

- F.Mignard & M.Froeschle, Astron. Astrophys. 354, 732-739 (2000).

"""
function h2fk5(ra, dec, dra, ddec, px, rv)
    r = Ref{Cdouble}()
    d = Ref{Cdouble}()
    dr = Ref{Cdouble}()
    dd = Ref{Cdouble}()
    p = Ref{Cdouble}()
    v = Ref{Cdouble}()
    ccall((:eraH2fk5, liberfa), Cvoid,
            (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
            ra, dec, dra, ddec, px, rv, r, d, dr, dd, p, v)
    return r[], d[], dr[], dd[], p[], v[]
end

