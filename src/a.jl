"""
    ab(pnat, v, s, bm1)

Apply aberration to transform natural direction into proper
direction.

### Given ###

- `pnat`: Natural direction to the source (unit vector)
- `v`: Observer barycentric velocity in units of c
- `s`: Distance between the Sun and the observer (au)
- `bm1`: ``\\sqrt{1-|v|^2}`` reciprocal of Lorenz factor

### Returned ###

- `ppr`: Proper direction to source (unit vector)

### Notes ###

1. The algorithm is based on Expr. (7.40) in the Explanatory
   Supplement (Urban & Seidelmann 2013), but with the following
   changes:

   -  Rigorous rather than approximate normalization is applied.

   -  The gravitational potential term from Expr. (7) in
      Klioner (2003) is added, taking into account only the Sun's
      contribution.  This has a maximum effect of about
      0.4 microarcsecond.

2. In almost all cases, the maximum accuracy will be limited by the
   supplied velocity.  For example, if the ERFA [`epv00`](@ref) function is
   used, errors of up to 5 microarcseconds could occur.

### References ###

- Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
    the Astronomical Almanac, 3rd ed., University Science Books
    (2013).

- Klioner, Sergei A., "A practical relativistic model for micro-
    arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).

### Called ###

- [`pdp`](@ref): scalar product of two p-vectors

"""
function ab(pnat, v, s, bm1)
    @checkdims 3 pnat v
    ppr = zeros(Cdouble, 3)
    ccall((:eraAb, liberfa), Cvoid,
          (Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Ptr{Cdouble}),
          pnat, v, s, bm1, ppr)
    ppr
end


"""
    ae2hd(az, el, phi)

Horizon to equatorial coordinates:  transform azimuth and altitude
to hour angle and declination.

### Given ###

- `az`: Azimuth
- `el`: Altitude (informally, elevation)
- `phi`: Site latitude

### Returned ###

- `ha`: Hour angle (local)
- `dec`: Declination

### Notes ###

1. All the arguments are angles in radians.

2. The sign convention for azimuth is north zero, east +pi/2.

3. HA is returned in the range +/-pi.  Declination is returned in
   the range +/-pi/2.

4. The latitude phi is pi/2 minus the angle between the Earth's
   rotation axis and the adopted zenith.  In many applications it
   will be sufficient to use the published geodetic latitude of the
   site.  In very precise (sub-arcsecond. applications, phi can be
   corrected for polar motion.

5. The azimuth az must be with respect to the rotational north pole,
   as opposed to the ITRS pole, and an azimuth with respect to north
   on a map of the Earth's surface will need to be adjusted for
   polar motion if sub-arcsecond accuracy is required.

6. Should the user wish to work with respect to the astronomical
   zenith rather than the geodetic zenith, phi will need to be
   adjusted for deflection of the vertical (often tens of
   arcseconds), and the zero point of ha will also be affected.

7. The transformation is the same as Ve = Ry(phi-pi/2)*Rz(pi)*Vh,
   where Ve and Vh are lefthanded unit vectors in the (ha,dec. and
   (az,el. systems respectively and Rz and Ry are rotations about
   first the z-axis and then the y-axis.  (n.b. Rz(pi. simply
   reverses the signs of the x and y components.. For efficiency,
   the algorithm is written out rather than calling other utility
   functions.  For applications that require even greater
   efficiency, additional savings are possible if constant terms
   such as functions of latitude are computed once and for all.

8. Again for efficiency, no range checking of arguments is carried
   out.
"""
function ae2hd(az, el, phi)
    h = Ref{Cdouble}()
    d = Ref{Cdouble}()
    ccall((:eraAe2hd, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
          az, el, phi, h, d)
    return h[], d[]
end

"""
    apcg(date1, date2, ebpv, ehp)

For a geocentric observer, prepare star-independent astrometry
parameters for transformations between ICRS and GCRS coordinates.
The Earth ephemeris is supplied by the caller.

The parameters produced by this function are required in the
parallax, light deflection and aberration parts of the astrometric
transformation chain.

### Given ###

- `date1`: TDB as a 2-part...
- `date2`: ...Julian Date (Note 1)
- `ebpv`: Earth barycentric pos/vel (au, au/day)
- `ehp`: Earth heliocentric position (au)

### Returned ###

- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: unchanged
    - `xpl`: unchanged
    - `ypl`: unchanged
    - `sphi`: unchanged
    - `cphi`: unchanged
    - `diurab`: unchanged
    - `l`: unchanged
    - `refa`: unchanged
    - `refb`: unchanged

### Notes ###

1. The TDB date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TDB)=2450123.7 could be expressed in any of these ways, among
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2. All the vectors are with respect to BCRS axes.

3. This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion,
   parallax, light deflection, and aberration.  From GCRS to CIRS
   comprises frame bias and precession-nutation.  From CIRS to
   observed takes account of Earth rotation, polar motion, diurnal
   aberration and parallax (unless subsumed into the ICRS <-> GCRS
   transformation), and atmospheric refraction.

4. The context structure astrom produced by this function is used by
   `atciq*` and `aticq*`.

### Called ###

- [`apcs`](@ref): astrometry parameters, ICRS-GCRS, space observer

"""
function apcg(date1, date2, ebpv, ehp)
    @checkdims 3 ehp
    _ebpv = array_to_cmatrix(ebpv; n=3)
    astrom = ASTROM()
    ccall((:eraApcg, liberfa), Cvoid,
          (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ref{ASTROM}),
          date1, date2, _ebpv, ehp, astrom)
    return astrom
end

"""
    apcg13(date1, date2)

For a geocentric observer, prepare star-independent astrometry
parameters for transformations between ICRS and GCRS coordinates.
The caller supplies the date, and ERFA models are used to predict
the Earth ephemeris.

The parameters produced by this function are required in the
parallax, light deflection and aberration parts of the astrometric
transformation chain.

### Given ###

- `date1`: TDB as a 2-part...
- `date2`: ...Julian Date (Note 1)

### Returned ###

- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: unchanged
    - `xpl`: unchanged
    - `ypl`: unchanged
    - `sphi`: unchanged
    - `cphi`: unchanged
    - `diurab`: unchanged
    - `l`: unchanged
    - `refa`: unchanged
    - `refb`: unchanged

### Notes ###

1. The TDB date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TDB)=2450123.7 could be expressed in any of these ways, among
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2. All the vectors are with respect to BCRS axes.

3. In cases where the caller wishes to supply his own Earth
   ephemeris, the function [`apcg`](@ref) can be used instead of the present
   function.

4. This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion,
   parallax, light deflection, and aberration.  From GCRS to CIRS
   comprises frame bias and precession-nutation.  From CIRS to
   observed takes account of Earth rotation, polar motion, diurnal
   aberration and parallax (unless subsumed into the ICRS <-> GCRS
   transformation), and atmospheric refraction.

5. The context structure astrom produced by this function is used by
   `atciq*` and `aticq*`.

### Called ###

- [`epv00`](@ref): Earth position and velocity
- [`apcg`](@ref): astrometry parameters, ICRS-GCRS, geocenter

"""
function apcg13(date1, date2)
    astrom = ASTROM()
    ccall((:eraApcg13, liberfa), Cvoid,
          (Cdouble, Cdouble, Ref{ASTROM}),
          date1, date2, astrom)
    return astrom
end

"""
    apci(date1, date2, ebpv, ehp, x, y, s)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and geocentric CIRS
coordinates.  The Earth ephemeris and CIP/CIO are supplied by the
caller.

The parameters produced by this function are required in the
parallax, light deflection, aberration, and bias-precession-nutation
parts of the astrometric transformation chain.

### Given ###

- `date1`: TDB as a 2-part...
- `date2`: ...Julian Date (Note 1)
- `ebpv`: Earth barycentric position/velocity (au, au/day)
- `ehp`: Earth heliocentric position (au)
- `x`, `y`: CIP X,Y (components of unit vector)
- `s`: The CIO locator s (radians)

### Returned ###

- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: unchanged
    - `xpl`: unchanged
    - `ypl`: unchanged
    - `sphi`: unchanged
    - `cphi`: unchanged
    - `diurab`: unchanged
    - `l`: unchanged
    - `refa`: unchanged
    - `refb`: unchanged

### Notes ###

1. The TDB date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TDB)=2450123.7 could be expressed in any of these ways, among
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2. All the vectors are with respect to BCRS axes.

3. In cases where the caller does not wish to provide the Earth
   ephemeris and CIP/CIO, the function [`apci13`](@ref) can be used instead
   of the present function.  This computes the required quantities
   using other ERFA functions.

4. This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion,
   parallax, light deflection, and aberration.  From GCRS to CIRS
   comprises frame bias and precession-nutation.  From CIRS to
   observed takes account of Earth rotation, polar motion, diurnal
   aberration and parallax (unless subsumed into the ICRS <-> GCRS
   transformation), and atmospheric refraction.

5. The context structure astrom produced by this function is used by
   `atciq*` and `aticq*`.

### Called ###

- [`apcg`](@ref): astrometry parameters, ICRS-GCRS, geocenter
- [`c2ixys`](@ref): celestial-to-intermediate matrix, given X,Y and s

"""
function apci(date1, date2, ebpv, ehp, x, y, s)
    @checkdims 3 ehp
    _ebpv = array_to_cmatrix(ebpv; n=3)
    astrom = ASTROM()
    ccall((:eraApci, liberfa), Cvoid,
          (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Ref{ASTROM}),
          date1, date2, _ebpv, ehp, x, y, s, astrom)
    return astrom
end

"""
    apci13(date1, date2)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and geocentric CIRS
coordinates.  The caller supplies the date, and ERFA models are used
to predict the Earth ephemeris and CIP/CIO.

The parameters produced by this function are required in the
parallax, light deflection, aberration, and bias-precession-nutation
parts of the astrometric transformation chain.

### Given ###

- `date1`: TDB as a 2-part...
- `date2`: ...Julian Date (Note 1)

### Returned ###

- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: unchanged
    - `xpl`: unchanged
    - `ypl`: unchanged
    - `sphi`: unchanged
    - `cphi`: unchanged
    - `diurab`: unchanged
    - `l`: unchanged
    - `refa`: unchanged
    - `refb`: unchanged
- `eo`: Equation of the origins (ERA-GST)

### Notes ###

1. The TDB date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TDB)=2450123.7 could be expressed in any of these ways, among
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2. All the vectors are with respect to BCRS axes.

3. In cases where the caller wishes to supply his own Earth
   ephemeris and CIP/CIO, the function [`apci`](@ref) can be used instead
   of the present function.

4. This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion,
   parallax, light deflection, and aberration.  From GCRS to CIRS
   comprises frame bias and precession-nutation.  From CIRS to
   observed takes account of Earth rotation, polar motion, diurnal
   aberration and parallax (unless subsumed into the ICRS <-> GCRS
   transformation), and atmospheric refraction.

5. The context structure astrom produced by this function is used by
   `atciq*` and `aticq*`.

### Called ###

- [`epv00`](@ref): Earth position and velocity
- [`pnm06a`](@ref): classical NPB matrix, IAU 2006/2000A
- [`bpn2xy`](@ref): extract CIP X,Y coordinates from NPB matrix
- [`s06`](@ref): the CIO locator s, given X,Y, IAU 2006
- [`apci`](@ref): astrometry parameters, ICRS-CIRS
- [`eors`](@ref): equation of the origins, given NPB matrix and s

"""
function apci13(date1, date2)
    astrom = ASTROM()
    eo = Ref{Cdouble}()
    ccall((:eraApci13, liberfa), Cvoid,
          (Cdouble, Cdouble, Ref{ASTROM}, Ref{Cdouble}),
          date1, date2, astrom, eo)
    return astrom, eo[]
end

"""
    apco(date1, date2, ebpv, ehp, x, y, s, theta, elong, phi, hm, xp, yp, sp, refa, refb)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and observed
coordinates.  The caller supplies the Earth ephemeris, the Earth
rotation information and the refraction constants as well as the
site coordinates.

### Given ###

- `date1`: TDB as a 2-part...
- `date2`: ...Julian Date (Note 1)
- `ebpv`: Earth barycentric PV (au, au/day, Note 2)
- `ehp`: Earth heliocentric P (au, Note 2)
- `x`, `y`: CIP X,Y (components of unit vector)
- `s`: The CIO locator s (radians)
- `theta`: Earth rotation angle (radians)
- `elong`: Longitude (radians, east +ve, Note 3)
- `phi`: Latitude (geodetic, radians, Note 3)
- `hm`: Height above ellipsoid (m, geodetic, Note 3)
- `xp`, `yp`: Polar motion coordinates (radians, Note 4)
- `sp`: The TIO locator s' (radians, Note 4)
- `refa`: Refraction constant A (radians, Note 5)
- `refb`: Refraction constant B (radians, Note 5)

### Returned ###

- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: Longitude + s' (radians)
    - `xp1`: Polar motion xp wrt local meridian (radians)
    - `yp1`: Polar motion yp wrt local meridian (radians)
    - `sphi`: Sine of geodetic latitude
    - `cphi`: Cosine of geodetic latitude
    - `diurab`: Magnitude of diurnal aberration vector
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: Refraction constant A (radians)
    - `refb`: Refraction constant B (radians)

### Notes ###

1. The TDB date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TDB)=2450123.7 could be expressed in any of these ways, among
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2. The vectors eb, eh, and all the astrom vectors, are with respect
   to BCRS axes.

3. The geographical coordinates are with respect to the [`WGS84`](@ref)
   reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN
   CONVENTION:  the longitude required by the present function is
   right-handed, i.e. east-positive, in accordance with geographical
   convention.

4. xp and yp are the coordinates (in radians) of the Celestial
   Intermediate Pole with respect to the International Terrestrial
   Reference System (see IERS Conventions), measured along the
   meridians 0 and 90 deg west respectively.  sp is the TIO locator
   s', in radians, which positions the Terrestrial Intermediate
   Origin on the equator.  For many applications, xp, yp and
   (especially) sp can be set to zero.

   Internally, the polar motion is stored in a form rotated onto the
   local meridian.

5. The refraction constants refa and refb are for use in a
   ``dZ = A*\\tan(Z)+B*\\tan^3(Z)`` model, where Z is the observed
   (i.e. refracted) zenith distance and dZ is the amount of
   refraction.

6. It is advisable to take great care with units, as even unlikely
   values of the input parameters are accepted and processed in
   accordance with the models used.

7. In cases where the caller does not wish to provide the Earth
   Ephemeris, the Earth rotation information and refraction
   constants, the function [`apco13`](@ref) can be used instead of the
   present function.  This starts from UTC and weather readings etc.
   and computes suitable values using other ERFA functions.

8. This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion,
   parallax, light deflection, and aberration.  From GCRS to CIRS
   comprises frame bias and precession-nutation.  From CIRS to
   observed takes account of Earth rotation, polar motion, diurnal
   aberration and parallax (unless subsumed into the ICRS <-> GCRS
   transformation), and atmospheric refraction.

9. The context structure astrom produced by this function is used by
   [`atioq`](@ref), [`atoiq`](@ref), [`atciq`](@ref) and [`aticq`](@ref).

### Called ###

- [`aper`](@ref): astrometry parameters: update ERA
- [`c2ixys`](@ref): celestial-to-intermediate matrix, given X,Y and s
- [`pvtob`](@ref): position/velocity of terrestrial station
- [`trxpv`](@ref): product of transpose of r-matrix and pv-vector
- [`apcs`](@ref): astrometry parameters, ICRS-GCRS, space observer
- [`cr`](@ref): copy r-matrix

"""
function apco(date1, date2, ebpv, ehp, x, y, s, theta, elong, phi, hm, xp, yp, sp, refa, refb)
    @checkdims 3 ehp
    _ebpv = array_to_cmatrix(ebpv; n=3)
    astrom = ASTROM()
    ccall((:eraApco, liberfa), Cvoid,
          (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cdouble,
           Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ref{ASTROM}),
          date1, date2, _ebpv, ehp, x, y, s, theta, elong, phi, hm, xp, yp, sp, refa, refb, astrom)
    return astrom
end

"""
    apco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between ICRS and observed
coordinates.  The caller supplies UTC, site coordinates, ambient air
conditions and observing wavelength, and ERFA models are used to
obtain the Earth ephemeris, CIP/CIO and refraction constants.

The parameters produced by this function are required in the
parallax, light deflection, aberration, and bias-precession-nutation
parts of the ICRS/CIRS transformations.

### Given ###

- `utc1`: UTC as a 2-part...
- `utc2`: ...quasi Julian Date (Notes 1,2)
- `dut1`: UT1-UTC (seconds, Note 3)
- `elong`: Longitude (radians, east +ve, Note 4)
- `phi`: Latitude (geodetic, radians, Note 4)
- `hm`: Height above ellipsoid (m, geodetic, Notes 4,6)
- `xp`, `yp`: Polar motion coordinates (radians, Note 5)
- `phpa`: Pressure at the observer (hPa = mB, Note 6)
- `tc`: Ambient temperature at the observer (deg C)
- `rh`: Relative humidity at the observer (range 0-1)
- `wl`: Wavelength (micrometers, Note 7)

### Returned ###

- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: Longitude + s' (radians)
    - `xp1`: Polar motion xp wrt local meridian (radians)
    - `yp1`: Polar motion yp wrt local meridian (radians)
    - `sphi`: Sine of geodetic latitude
    - `cphi`: Cosine of geodetic latitude
    - `diurab`: Magnitude of diurnal aberration vector
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: Refraction constant A (radians)
    - `refb`: Refraction constant B (radians)
- `eo`: Equation of the origins (ERA-GST)

### Notes ###

1.  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function [`dtf2d`](@ref) to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

2.  The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the
    future to be trusted.  See [`dat`](@ref) for further details.

3.  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

4.  The geographical coordinates are with respect to the [`WGS84`](@ref)
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

5.  The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many
    applications, xp and yp can be set to zero.

    Internally, the polar motion is stored in a form rotated onto
    the local meridian.

6.  If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB), is
    available, an adequate estimate of hm can be obtained from the
    expression

    ```
    hm = -29.3 * tsl * log ( phpa / 1013.25 );
    ```

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    ```
    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    ```

    Note, however, that the refraction is nearly proportional to
    the pressure and that an accurate phpa value is important for
    precise work.

7.  The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

8.  It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

9.  In cases where the caller wishes to supply his own Earth
    ephemeris, Earth rotation information and refraction constants,
    the function [`apco`](@ref) can be used instead of the present function.

10. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

    Those with names ending in "13" use contemporary ERFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

11. The context structure astrom produced by this function is used
    by [`atioq`](@ref), [`atoiq`](@ref), `atciq*` and `aticq*`.

### Called ###

- [`utctai`](@ref): UTC to TAI
- [`taitt`](@ref): TAI to TT
- [`utcut1`](@ref): UTC to UT1
- [`epv00`](@ref): Earth position and velocity
- [`pnm06a`](@ref): classical NPB matrix, IAU 2006/2000A
- [`bpn2xy`](@ref): extract CIP X,Y coordinates from NPB matrix
- [`s06`](@ref): the CIO locator s, given X,Y, IAU 2006
- [`era00`](@ref): Earth rotation angle, IAU 2000
- [`sp00`](@ref): the TIO locator s', IERS 2000
- [`refco`](@ref): refraction constants for given ambient conditions
- [`apco`](@ref): astrometry parameters, ICRS-observed
- [`eors`](@ref): equation of the origins, given NPB matrix and s

"""
function apco13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)
    astrom = ASTROM()
    eo = Ref{Cdouble}()
    i = ccall((:eraApco13, liberfa), Cint,
              (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
               Cdouble, Cdouble, Cdouble, Cdouble, Ref{ASTROM}, Ref{Cdouble}),
              utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl, astrom, eo)
    if i == -1
        throw(ERFAException("unacceptable date"))
    elseif i == +1
        @warn "dubious year"
    end
    return astrom, eo[]
end

"""
    apcs(date1, date2, pv, ebpv, ehp)

For an observer whose geocentric position and velocity are known,
prepare star-independent astrometry parameters for transformations
between ICRS and GCRS.  The Earth ephemeris is supplied by the
caller.

The parameters produced by this function are required in the space
motion, parallax, light deflection and aberration parts of the
astrometric transformation chain.

### Given ###

- `date1`: TDB as a 2-part...
- `date2`: ...Julian Date (Note 1)
- `pv`: Observer's geocentric pos/vel (m, m/s)
- `ebpv`: Earth barycentric PV (au, au/day)
- `ehp`: Earth heliocentric P (au)

### Returned ###

- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: unchanged
    - `xpl`: unchanged
    - `ypl`: unchanged
    - `sphi`: unchanged
    - `cphi`: unchanged
    - `diurab`: unchanged
    - `l`: unchanged
    - `refa`: unchanged
    - `refb`: unchanged

### Notes ###

1. The TDB date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TDB)=2450123.7 could be expressed in any of these ways, among
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2. All the vectors are with respect to BCRS axes.

3. Providing separate arguments for (i) the observer's geocentric
   position and velocity and (ii) the Earth ephemeris is done for
   convenience in the geocentric, terrestrial and Earth orbit cases.
   For deep space applications it maybe more convenient to specify
   zero geocentric position and velocity and to supply the
   observer's position and velocity information directly instead of
   with respect to the Earth.  However, note the different units:
   m and m/s for the geocentric vectors, au and au/day for the
   heliocentric and barycentric vectors.

4. In cases where the caller does not wish to provide the Earth
   ephemeris, the function [`apcs13`](@ref) can be used instead of the
   present function.  This computes the Earth ephemeris using the
   ERFA function [`epv00`](@ref).

5. This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion,
   parallax, light deflection, and aberration.  From GCRS to CIRS
   comprises frame bias and precession-nutation.  From CIRS to
   observed takes account of Earth rotation, polar motion, diurnal
   aberration and parallax (unless subsumed into the ICRS <-> GCRS
   transformation), and atmospheric refraction.

6. The context structure astrom produced by this function is used by
   `atciq*` and `aticq*`.

### Called ###

- [`erfa_cp`](@ref): copy p-vector
- [`pm`](@ref): modulus of p-vector
- [`pn`](@ref): decompose p-vector into modulus and direction
- [`ir`](@ref): initialize r-matrix to identity

"""
function apcs(date1, date2, pv, ebpv, ehp)
    @checkdims 3 ehp
    _pv = array_to_cmatrix(pv; n=3)
    _ebpv = array_to_cmatrix(ebpv; n=3)
    astrom = ASTROM()
    ccall((:eraApcs, liberfa), Cvoid,
          (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{ASTROM}),
          date1, date2, _pv, _ebpv, ehp, astrom)
    return astrom
end

"""
    apcs13(date1, date2, pv)

For an observer whose geocentric position and velocity are known,
prepare star-independent astrometry parameters for transformations
between ICRS and GCRS.  The Earth ephemeris is from ERFA models.

The parameters produced by this function are required in the space
motion, parallax, light deflection and aberration parts of the
astrometric transformation chain.

### Given ###

- `date1`: TDB as a 2-part...
- `date2`: ...Julian Date (Note 1)
- `pv`: Observer's geocentric pos/vel (Note 3)

### Returned ###

- EraASTROM*   star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: unchanged
    - `xpl`: unchanged
    - `ypl`: unchanged
    - `sphi`: unchanged
    - `cphi`: unchanged
    - `diurab`: unchanged
    - `l`: unchanged
    - `refa`: unchanged
    - `refb`: unchanged

### Notes ###

1. The TDB date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TDB)=2450123.7 could be expressed in any of these ways, among
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2. All the vectors are with respect to BCRS axes.

3. The observer's position and velocity pv are geocentric but with
   respect to BCRS axes, and in units of m and m/s.  No assumptions
   are made about proximity to the Earth, and the function can be
   used for deep space applications as well as Earth orbit and
   terrestrial.

4. In cases where the caller wishes to supply his own Earth
   ephemeris, the function [`apcs`](@ref) can be used instead of the present
   function.

5. This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion,
   parallax, light deflection, and aberration.  From GCRS to CIRS
   comprises frame bias and precession-nutation.  From CIRS to
   observed takes account of Earth rotation, polar motion, diurnal
   aberration and parallax (unless subsumed into the ICRS <-> GCRS
   transformation), and atmospheric refraction.

6. The context structure astrom produced by this function is used by
   `atciq*` and `aticq*`.

### Called ###

- [`epv00`](@ref): Earth position and velocity
- [`apcs`](@ref): astrometry parameters, ICRS-GCRS, space observer

"""
function apcs13(date1, date2, pv)
    _pv = array_to_cmatrix(pv; n=3)
    astrom = ASTROM()
    ccall((:eraApcs13, liberfa), Cvoid,
          (Cdouble, Cdouble, Ptr{Cdouble}, Ref{ASTROM}),
          date1, date2, _pv, astrom)
    return astrom
end

"""
    aper(theta, astrom)

In the star-independent astrometry parameters, update only the
Earth rotation angle, supplied by the caller explicitly.

### Given ###

- `theta`: Earth rotation angle (radians, Note 2)
- `astrom`: Star-independent astrometry parameters:
    - `pmt`: unchanged
    - `eb`: unchanged
    - `eh`: unchanged
    - `em`: unchanged
    - `v`: unchanged
    - `bm1`: unchanged
    - `bpn`: unchanged
    - `along`: Longitude + s' (radians)
    - `xpl`: unchanged
    - `ypl`: unchanged
    - `sphi`: unchanged
    - `cphi`: unchanged
    - `diurab`: unchanged
    - `l`: unchanged
    - `refa`: unchanged
    - `refb`: unchanged

### Returned ###

- `astrom`: Star-independent astrometry parameters:
    - `pmt`: unchanged
    - `eb`: unchanged
    - `eh`: unchanged
    - `em`: unchanged
    - `v`: unchanged
    - `bm1`: unchanged
    - `bpn`: unchanged
    - `along`: unchanged
    - `xpl`: unchanged
    - `ypl`: unchanged
    - `sphi`: unchanged
    - `cphi`: unchanged
    - `diurab`: unchanged
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: unchanged
    - `refb`: unchanged

### Notes ###

1. This function exists to enable sidereal-tracking applications to
   avoid wasteful recomputation of the bulk of the astrometry
   parameters:  only the Earth rotation is updated.

2. For targets expressed as equinox based positions, such as
   classical geocentric apparent (RA,Dec), the supplied theta can be
   Greenwich apparent sidereal time rather than Earth rotation
   angle.

3. The function [`aper13`](@ref) can be used instead of the present
   function, and starts from UT1 rather than ERA itself.

4. This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion,
   parallax, light deflection, and aberration.  From GCRS to CIRS
   comprises frame bias and precession-nutation.  From CIRS to
   observed takes account of Earth rotation, polar motion, diurnal
   aberration and parallax (unless subsumed into the ICRS <-> GCRS
   transformation), and atmospheric refraction.

"""
function aper(theta, astrom)
    ccall((:eraAper, liberfa), Cvoid,
          (Cdouble, Ref{ASTROM}),
          theta, astrom)
    return astrom
end

"""
    aper13(ut11, ut12, astrom)

In the star-independent astrometry parameters, update only the
Earth rotation angle.  The caller provides UT1, (n.b. not UTC).

### Given ###

- `ut11`: UT1 as a 2-part...
- `ut12`: ...Julian Date (Note 1)
- `astrom`: Star-independent astrometry parameters:
    - `pmt`: unchanged
    - `eb`: unchanged
    - `eh`: unchanged
    - `em`: unchanged
    - `v`: unchanged
    - `bm1`: unchanged
    - `bpn`: unchanged
    - `along`: Longitude + s' (radians)
    - `xpl`: unchanged
    - `ypl`: unchanged
    - `sphi`: unchanged
    - `cphi`: unchanged
    - `diurab`: unchanged
    - `l`: unchanged
    - `refa`: unchanged
    - `refb`: unchanged

### Returned ###

- `astrom`: Star-independent astrometry parameters:
    - `pmt`: unchanged
    - `eb`: unchanged
    - `eh`: unchanged
    - `em`: unchanged
    - `v`: unchanged
    - `bm1`: unchanged
    - `bpn`: unchanged
    - `along`: unchanged
    - `xpl`: unchanged
    - `ypl`: unchanged
    - `sphi`: unchanged
    - `cphi`: unchanged
    - `diurab`: unchanged
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: unchanged
    - `refb`: unchanged

### Notes ###

1. The UT1 date (n.b. not UTC) ut11+ut12 is a Julian Date,
   apportioned in any convenient way between the arguments ut11 and
   ut12.  For example, JD(UT1)=2450123.7 could be expressed in any
   of these ways, among others:

   |  `ut11`   | `ut12`  | Method               |
   |:----------|:--------|----------------------|
   | 2450123.7 |     0.0 | JD                   |
   | 2451545.0 | -1421.3 | J2000                |
   | 2400000.5 | 50123.2 | MJD                  |
   | 2450123.5 |     0.2 | date & time          |

   The JD method is the most natural and convenient to use in cases
   where the loss of several decimal digits of resolution is
   acceptable.  The J2000 and MJD methods are good compromises
   between resolution and convenience.  The date & time method is
   best matched to the algorithm used:  maximum precision is
   delivered when the ut11 argument is for 0hrs UT1 on the day in
   question and the ut12 argument lies in the range 0 to 1, or vice
   versa.

2. If the caller wishes to provide the Earth rotation angle itself,
   the function [`aper`](@ref) can be used instead.  One use of this
   technique is to substitute Greenwich apparent sidereal time and
   thereby to support equinox based transformations directly.

3. This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion,
   parallax, light deflection, and aberration.  From GCRS to CIRS
   comprises frame bias and precession-nutation.  From CIRS to
   observed takes account of Earth rotation, polar motion, diurnal
   aberration and parallax (unless subsumed into the ICRS <-> GCRS
   transformation), and atmospheric refraction.

### Called ###

- [`aper`](@ref): astrometry parameters: update ERA
- [`era00`](@ref): Earth rotation angle, IAU 2000

"""
function aper13(ut11, ut12, astrom)
    ccall((:eraAper13, liberfa), Cvoid,
          (Cdouble, Cdouble, Ref{ASTROM}),
          ut11, ut12, astrom)
    return astrom
end

"""
    apio(sp, theta, elong, phi, hm, xp, yp, refa, refb)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between CIRS and observed
coordinates.  The caller supplies the Earth orientation information
and the refraction constants as well as the site coordinates.

### Given ###

- `sp`: The TIO locator s' (radians, Note 1)
- `theta`: Earth rotation angle (radians)
- `elong`: Longitude (radians, east +ve, Note 2)
- `phi`: Geodetic latitude (radians, Note 2)
- `hm`: Height above ellipsoid (m, geodetic Note 2)
- `xp`, `yp`: Polar motion coordinates (radians, Note 3)
- `refa`: Refraction constant A (radians, Note 4)
- `refb`: Refraction constant B (radians, Note 4)

### Returned ###

- `astrom`: Star-independent astrometry parameters:
    - `pmt`: unchanged
    - `eb`: unchanged
    - `eh`: unchanged
    - `em`: unchanged
    - `v`: unchanged
    - `bm1`: unchanged
    - `bpn`: unchanged
    - `along`: Longitude + s' (radians)
    - `xp1`: Polar motion xp wrt local meridian (radians)
    - `yp1`: Polar motion yp wrt local meridian (radians)
    - `sphi`: Sine of geodetic latitude
    - `cphi`: Cosine of geodetic latitude
    - `diurab`: Magnitude of diurnal aberration vector
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: Refraction constant A (radians)
    - `refb`: Refraction constant B (radians)

### Notes ###

1. sp, the TIO locator s', is a tiny quantity needed only by the
   most precise applications.  It can either be set to zero or
   predicted using the ERFA function [`sp00`](@ref).

2. The geographical coordinates are with respect to the [`WGS84`](@ref)
   reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
   longitude required by the present function is east-positive
   (i.e. right-handed), in accordance with geographical convention.

3. The polar motion xp,yp can be obtained from IERS bulletins.  The
   values are the coordinates (in radians) of the Celestial
   Intermediate Pole with respect to the International Terrestrial
   Reference System (see IERS Conventions 2003), measured along the
   meridians 0 and 90 deg west respectively.  For many applications,
   xp and yp can be set to zero.

   Internally, the polar motion is stored in a form rotated onto the
   local meridian.

4. The refraction constants refa and refb are for use in a
   ``dZ = A*\\tan(Z)+B*\\tan^3(Z)`` model, where Z is the observed
   (i.e. refracted) zenith distance and dZ is the amount of
   refraction.

5. It is advisable to take great care with units, as even unlikely
   values of the input parameters are accepted and processed in
   accordance with the models used.

6. In cases where the caller does not wish to provide the Earth
   rotation information and refraction constants, the function
   [`apio13`](@ref) can be used instead of the present function.  This
   starts from UTC and weather readings etc. and computes suitable
   values using other ERFA functions.

7. This is one of several functions that inserts into the astrom
   structure star-independent parameters needed for the chain of
   astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

   The various functions support different classes of observer and
   portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

   Those with names ending in "13" use contemporary ERFA models to
   compute the various ephemerides.  The others accept ephemerides
   supplied by the caller.

   The transformation from ICRS to GCRS covers space motion,
   parallax, light deflection, and aberration.  From GCRS to CIRS
   comprises frame bias and precession-nutation.  From CIRS to
   observed takes account of Earth rotation, polar motion, diurnal
   aberration and parallax (unless subsumed into the ICRS <-> GCRS
   transformation), and atmospheric refraction.

8. The context structure astrom produced by this function is used by
   [`atioq`](@ref) and [`atoiq`](@ref).

### Called ###

- [`pvtob`](@ref): position/velocity of terrestrial station
- [`aper`](@ref): astrometry parameters: update ERA

"""
function apio(sp, theta, elong, phi, hm, xp, yp, refa, refb)
    astrom = ASTROM()
    ccall((:eraApio, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ref{ASTROM}),
          sp, theta, elong, phi, hm, xp, yp, refa, refb, astrom)
    return astrom
end

"""
    apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)

For a terrestrial observer, prepare star-independent astrometry
parameters for transformations between CIRS and observed
coordinates.  The caller supplies UTC, site coordinates, ambient air
conditions and observing wavelength.

### Given ###

- `utc1`: UTC as a 2-part...
- `utc2`: ...quasi Julian Date (Notes 1,2)
- `dut1`: UT1-UTC (seconds)
- `elong`: Longitude (radians, east +ve, Note 3)
- `phi`: Geodetic latitude (radians, Note 3)
- `hm`: Height above ellipsoid (m, geodetic Notes 4,6)
- `xp`, `yp`: Polar motion coordinates (radians, Note 5)
- `phpa`: Pressure at the observer (hPa = mB, Note 6)
- `tc`: Ambient temperature at the observer (deg C)
- `rh`: Relative humidity at the observer (range 0-1)
- `wl`: Wavelength (micrometers, Note 7)

### Returned ###

- `astrom`: Star-independent astrometry parameters:
    - `pmt`: unchanged
    - `eb`: unchanged
    - `eh`: unchanged
    - `em`: unchanged
    - `v`: unchanged
    - `bm1`: unchanged
    - `bpn`: unchanged
    - `along`: Longitude + s' (radians)
    - `xp1`: Polar motion xp wrt local meridian (radians)
    - `yp1`: Polar motion yp wrt local meridian (radians)
    - `sphi`: Sine of geodetic latitude
    - `cphi`: Cosine of geodetic latitude
    - `diurab`: Magnitude of diurnal aberration vector
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: Refraction constant A (radians)
    - `refb`: Refraction constant B (radians)

### Notes ###

1.  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function [`dtf2d`](@ref) to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

2.  The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the future
    to be trusted.  See [`dat`](@ref) for further details.

3.  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

4.  The geographical coordinates are with respect to the [`WGS84`](@ref)
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

5.  The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many applications,
    xp and yp can be set to zero.

    Internally, the polar motion is stored in a form rotated onto
    the local meridian.

6.  If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB), is
    available, an adequate estimate of hm can be obtained from the
    expression

    ```
    hm = -29.3 * tsl * log ( phpa / 1013.25 );
    ```

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    ```
    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    ```

    Note, however, that the refraction is nearly proportional to the
    pressure and that an accurate phpa value is important for
    precise work.

7.  The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

8.  It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

9.  In cases where the caller wishes to supply his own Earth
    rotation information and refraction constants, the function
    `apc*` can be used instead of the present function.

10. This is one of several functions that inserts into the astrom
    structure star-independent parameters needed for the chain of
    astrometric transformations ICRS <-> GCRS <-> CIRS <-> observed.

    The various functions support different classes of observer and
    portions of the transformation chain:

   | Functions                       | Observer       | Transformation          |
   | :-------------------------------| :------------- | :---------------------- |
   | [`apcg`](@ref) [`apcg13`](@ref) | geocentric     | ICRS <-> GCRS           |
   | [`apci`](@ref) [`apci13`](@ref) | terrestrial    | ICRS <-> CIRS           |
   | [`apco`](@ref) [`apco13`](@ref) | terrestrial    | ICRS <-> observed       |
   | [`apcs`](@ref) [`apcs13`](@ref) | space          | ICRS <-> GCRS           |
   | [`aper`](@ref) [`aper13`](@ref) | terrestrial    | update Earth rotation   |
   | [`apio`](@ref) [`apio13`](@ref) | terrestrial    | CIRS <-> observed       |

    Those with names ending in "13" use contemporary ERFA models to
    compute the various ephemerides.  The others accept ephemerides
    supplied by the caller.

    The transformation from ICRS to GCRS covers space motion,
    parallax, light deflection, and aberration.  From GCRS to CIRS
    comprises frame bias and precession-nutation.  From CIRS to
    observed takes account of Earth rotation, polar motion, diurnal
    aberration and parallax (unless subsumed into the ICRS <-> GCRS
    transformation), and atmospheric refraction.

11. The context structure astrom produced by this function is used
    by [`atioq`](@ref) and [`atoiq`](@ref).

### Called ###

- [`utctai`](@ref): UTC to TAI
- [`taitt`](@ref): TAI to TT
- [`utcut1`](@ref): UTC to UT1
- [`sp00`](@ref): the TIO locator s', IERS 2000
- [`era00`](@ref): Earth rotation angle, IAU 2000
- [`refco`](@ref): refraction constants for given ambient conditions
- [`apio`](@ref): astrometry parameters, CIRS-observed

"""
function apio13(utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)
    astrom = ASTROM()
    i = ccall((:eraApio13, liberfa), Cint,
              (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
               Cdouble, Cdouble, Cdouble, Cdouble, Ref{ASTROM}),
              utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl, astrom)
    if i == -1
        throw(ERFAException("unacceptable date"))
    elseif i == +1
        @warn "dubious year"
    end
    return astrom
end

"""
    atci13(rc, dc, pr, pd, px, rv, date1, date2)

Transform ICRS star data, epoch J2000.0, to CIRS.

### Given ###

- `rc`: ICRS right ascension at J2000.0 (radians, Note 1)
- `dc`: ICRS declination at J2000.0 (radians, Note 1)
- `pr`: RA proper motion (radians/year; Note 2)
- `pd`: Dec proper motion (radians/year)
- `px`: Parallax (arcsec)
- `rv`: Radial velocity (km/s, +ve if receding)
- `date1`: TDB as a 2-part...
- `date2`: ...Julian Date (Note 3)

### Returned ###

- `ri`, `di`: CIRS geocentric RA,Dec (radians)
- `eo`: Equation of the origins (ERA-GST, Note 5)

### Notes ###

1. Star data for an epoch other than J2000.0 (for example from the
   Hipparcos catalog, which has an epoch of J1991.25) will require a
   preliminary call to [`pmsafe`](@ref) before use.

2. The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

3. The TDB date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TDB)=2450123.7 could be expressed in any of these ways, among
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

4. The available accuracy is better than 1 milliarcsecond, limited
   mainly by the precession-nutation model that is used, namely
   IAU 2000A/2006.  Very close to solar system bodies, additional
   errors of up to several milliarcseconds can occur because of
   unmodeled light deflection;  however, the Sun's contribution is
   taken into account, to first order.  The accuracy limitations of
   the ERFA function [`epv00`](@ref) (used to compute Earth position and
   velocity) can contribute aberration errors of up to
   5 microarcseconds.  Light deflection at the Sun's limb is
   uncertain at the 0.4 mas level.

5. Should the transformation to (equinox based) apparent place be
   required rather than (CIO based) intermediate place, subtract the
   equation of the origins from the returned right ascension:
   RA = RI - EO. (The [`anp`](@ref) function can then be applied, as
   required, to keep the result in the conventional 0-2pi range.)

### Called ###

- [`apci13`](@ref): astrometry parameters, ICRS-CIRS, 2013
- [`atciq`](@ref): quick ICRS to CIRS

"""
function atci13(rc, dc, pr, pd, px, rv, date1, date2)
    ri = Ref{Cdouble}()
    di = Ref{Cdouble}()
    eo = Ref{Cdouble}()
    ccall((:eraAtci13, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
           Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
          rc, dc, pr, pd, px, rv, date1, date2, ri, di, eo)
    return ri[], di[], eo[]
end

"""
    atciq(rc, dc, pr, pd, px, rv, astrom)

Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
star-independent astrometry parameters.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions `apci[13]`, `apcg[13]`, `apco[13]` or `apcs[13]`.

If the parallax and proper motions are zero the [`atciqz`](@ref) function
can be used instead.

### Given ###

- `rc`, `dc`: ICRS RA,Dec at J2000.0 (radians)
- `pr`: RA proper motion (radians/year; Note 3)
- `pd`: Dec proper motion (radians/year)
- `px`: Parallax (arcsec)
- `rv`: Radial velocity (km/s, +ve if receding)
- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: Longitude + s' (radians)
    - `xp1`: Polar motion xp wrt local meridian (radians)
    - `yp1`: Polar motion yp wrt local meridian (radians)
    - `sphi`: Sine of geodetic latitude
    - `cphi`: Cosine of geodetic latitude
    - `diurab`: Magnitude of diurnal aberration vector
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: Refraction constant A (radians)
    - `refb`: Refraction constant B (radians)

### Returned ###

- `ri`, `di`: CIRS RA,Dec (radians)

### Notes ###

1. All the vectors are with respect to BCRS axes.

2. Star data for an epoch other than J2000.0 (for example from the
   Hipparcos catalog, which has an epoch of J1991.25) will require a
   preliminary call to [`pmsafe`](@ref) before use.

3. The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

### Called ###

- [`pmpx`](@ref): proper motion and parallax
- [`ldsun`](@ref): light deflection by the Sun
- [`ab`](@ref): stellar aberration
- [`rxp`](@ref): product of r-matrix and pv-vector
- [`c2s`](@ref): p-vector to spherical
- [`anp`](@ref): normalize angle into range 0 to 2pi

"""
function atciq(rc, dc, pr, pd, px, rv, astrom)
    ri = Ref{Cdouble}()
    di = Ref{Cdouble}()
    ccall((:eraAtciq, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ref{ASTROM},
           Ref{Cdouble}, Ref{Cdouble}),
          rc, dc, pr, pd, px, rv, astrom, ri, di)
    return ri[], di[]
end

"""
    atciqn(rc, dc, pr, pd, px, rv, astrom, b::Vector{LDBODY})

Quick ICRS, epoch J2000.0, to CIRS transformation, given precomputed
star-independent astrometry parameters plus a list of light-
deflecting bodies.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions `apci[13]`, `apcg[13]`, `apco[13]` or `apcs[13]`.

If the only light-deflecting body to be taken into account is the
Sun, the [`atciq`](@ref) function can be used instead.  If in addition the
parallax and proper motions are zero, the [`atciqz`](@ref) function can be
used.

### Given ###

- `rc`, `dc`: ICRS RA,Dec at J2000.0 (radians)
- `pr`: RA proper motion (radians/year; Note 3)
- `pd`: Dec proper motion (radians/year)
- `px`: Parallax (arcsec)
- `rv`: Radial velocity (km/s, +ve if receding)
- EraASTROM*   star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: Longitude + s' (radians)
    - `xp1`: Polar motion xp wrt local meridian (radians)
    - `yp1`: Polar motion yp wrt local meridian (radians)
    - `sphi`: Sine of geodetic latitude
    - `cphi`: Cosine of geodetic latitude
    - `diurab`: Magnitude of diurnal aberration vector
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: Refraction constant A (radians)
    - `refb`: Refraction constant B (radians)
- `n`: Number of bodies (Note 3)
- `b::Vector{LDBODY}`: Data for each of the n bodies (Notes 3,4):
    - `bm`: Mass of the body (solar masses, Note 5)
    - `dl`: Deflection limiter (Note 6)
    - `pv`: Barycentric PV of the body (au, au/day)

### Returned ###

- `ri`, `di`: CIRS RA,Dec (radians)

### Notes ###

1. Star data for an epoch other than J2000.0 (for example from the
   Hipparcos catalog, which has an epoch of J1991.25) will require a
   preliminary call to [`pmsafe`](@ref) before use.

2. The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

3. The struct b contains n entries, one for each body to be
   considered.  If n = 0, no gravitational light deflection will be
   applied, not even for the Sun.

4. The struct b should include an entry for the Sun as well as for
   any planet or other body to be taken into account.  The entries
   should be in the order in which the light passes the body.

5. In the entry in the b struct for body i, the mass parameter
   b[i].bm can, as required, be adjusted in order to allow for such
   effects as quadrupole field.

6. The deflection limiter parameter b[i].dl is phi^2/2, where phi is
   the angular separation (in radians) between star and body at
   which limiting is applied.  As phi shrinks below the chosen
   threshold, the deflection is artificially reduced, reaching zero
   for phi = 0.   Example values suitable for a terrestrial
   observer, together with masses, are as follows:

   | `body i` | `b[i].bm`  | `b[i].dl` |
   |:---------|:-----------|:----------|
   |  Sun     | 1.0        |  6e-6     |
   |  Jupiter | 0.00095435 |  3e-9     |
   |  Saturn  | 0.00028574 |  3e-10    |

7. For efficiency, validation of the contents of the b array is
   omitted.  The supplied masses must be greater than zero, the
   position and velocity vectors must be right, and the deflection
   limiter greater than zero.

### Called ###

- [`pmpx`](@ref): proper motion and parallax
- [`ldn`](@ref): light deflection by n bodies
- [`ab`](@ref): stellar aberration
- [`rxp`](@ref): product of r-matrix and pv-vector
- [`c2s`](@ref): p-vector to spherical
- [`anp`](@ref): normalize angle into range 0 to 2pi

"""
function atciqn(rc, dc, pr, pd, px, rv, astrom, b::Vector{LDBODY})
    ri = Ref{Cdouble}()
    di = Ref{Cdouble}()
    n = length(b)
    ccall((:eraAtciqn, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ref{ASTROM}, Cint,
           Ref{LDBODY}, Ref{Cdouble}, Ref{Cdouble}),
          rc, dc, pr, pd, px, rv, astrom, n, b, ri, di)
    return ri[], di[]
end

"""
    atciqz(rc, dc, astrom)

Quick ICRS to CIRS transformation, given precomputed star-
independent astrometry parameters, and assuming zero parallax and
proper motion.

Use of this function is appropriate when efficiency is important and
where many star positions are to be transformed for one date.  The
star-independent parameters can be obtained by calling one of the
functions `apci[13]`, `apcg[13]`, `apco[13]` or `apcs[13]`.

The corresponding function for the case of non-zero parallax and
proper motion is [`atciq`](@ref).

### Given ###

- `rc`, `dc`: ICRS astrometric RA,Dec (radians)
- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: Longitude + s' (radians)
    - `xp1`: Polar motion xp wrt local meridian (radians)
    - `yp1`: Polar motion yp wrt local meridian (radians)
    - `sphi`: Sine of geodetic latitude
    - `cphi`: Cosine of geodetic latitude
    - `diurab`: Magnitude of diurnal aberration vector
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: Refraction constant A (radians)
    - `refb`: Refraction constant B (radians)

### Returned ###

- `ri`, `di`: CIRS RA,Dec (radians)

### Note ###

   All the vectors are with respect to BCRS axes.

### References ###

- Urban, S. & Seidelmann, P. K. (eds), Explanatory Supplement to
    the Astronomical Almanac, 3rd ed., University Science Books
    (2013).

- Klioner, Sergei A., "A practical relativistic model for micro-
    arcsecond astrometry in space", Astr. J. 125, 1580-1597 (2003).

### Called ###

- [`s2c`](@ref): spherical coordinates to unit vector
- [`ldsun`](@ref): light deflection due to Sun
- [`ab`](@ref): stellar aberration
- [`rxp`](@ref): product of r-matrix and p-vector
- [`c2s`](@ref): p-vector to spherical
- [`anp`](@ref): normalize angle into range +/- pi

"""
function atciqz(rc, dc, astrom)
    ri = Ref{Cdouble}()
    di = Ref{Cdouble}()
    ccall((:eraAtciqz, liberfa), Cvoid,
          (Cdouble, Cdouble, Ref{ASTROM}, Ref{Cdouble}, Ref{Cdouble}),
          rc, dc, astrom, ri, di)
    return ri[], di[]
end

"""
    atco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)

ICRS RA,Dec to observed place.  The caller supplies UTC, site
coordinates, ambient air conditions and observing wavelength.

ERFA models are used for the Earth ephemeris, bias-precession-
nutation, Earth orientation and refraction.

### Given ###

- `rc`, `dc`: ICRS right ascension at J2000.0 (radians, Note 1)
- `pr`: RA proper motion (radians/year; Note 2)
- `pd`: Dec proper motion (radians/year)
- `px`: Parallax (arcsec)
- `rv`: Radial velocity (km/s, +ve if receding)
- `utc1`: UTC as a 2-part...
- `utc2`: ...quasi Julian Date (Notes 3-4)
- `dut1`: UT1-UTC (seconds, Note 5)
- `elong`: Longitude (radians, east +ve, Note 6)
- `phi`: Latitude (geodetic, radians, Note 6)
- `hm`: Height above ellipsoid (m, geodetic, Notes 6,8)
- `xp`, `yp`: Polar motion coordinates (radians, Note 7)
- `phpa`: Pressure at the observer (hPa = mB, Note 8)
- `tc`: Ambient temperature at the observer (deg C)
- `rh`: Relative humidity at the observer (range 0-1)
- `wl`: Wavelength (micrometers, Note 9)

### Returned ###

- `aob`: Observed azimuth (radians: N=0,E=90)
- `zob`: Observed zenith distance (radians)
- `hob`: Observed hour angle (radians)
- `dob`: Observed declination (radians)
- `rob`: Observed right ascension (CIO-based, radians)
- `eo`: Equation of the origins (ERA-GST)

### Notes ###

1.  Star data for an epoch other than J2000.0 (for example from the
    Hipparcos catalog, which has an epoch of J1991.25) will require
    a preliminary call to [`pmsafe`](@ref) before use.

2.  The proper motion in RA is dRA/dt rather than cos(Dec)*dRA/dt.

3.  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function [`dtf2d`](@ref) to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

4.  The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the
    future to be trusted.  See [`dat`](@ref) for further details.

5.  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

6.  The geographical coordinates are with respect to the [`WGS84`](@ref)
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

7.  The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many
    applications, xp and yp can be set to zero.

8.  If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB),
    is available, an adequate estimate of hm can be obtained from
    the expression

    ```
    hm = -29.3 * tsl * log ( phpa / 1013.25 );
    ```

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    ```
    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    ```

    Note, however, that the refraction is nearly proportional to
    the pressure and that an accurate phpa value is important for
    precise work.

9.  The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

10. The accuracy of the result is limited by the corrections for
    refraction, which use a simple ``A*tan(z) + B*tan^3(z)`` model.
    Providing the meteorological parameters are known accurately and
    there are no gross local effects, the predicted observed
    coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    (radio) for a zenith distance of less than 70 degrees, better
    than 30 arcsec (optical or radio) at 85 degrees and better
    than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

    Without refraction, the complementary functions [`atco13`](@ref) and
    [`atoc13`](@ref) are self-consistent to better than 1 microarcsecond
    all over the celestial sphere.  With refraction included,
    consistency falls off at high zenith distances, but is still
    better than 0.05 arcsec at 85 degrees.

11. "Observed" Az,ZD means the position that would be seen by a
    perfect geodetically aligned theodolite.  (Zenith distance is
    used rather than altitude in order to reflect the fact that no
    allowance is made for depression of the horizon.)  This is
    related to the observed HA,Dec via the standard rotation, using
    the geodetic latitude (corrected for polar motion), while the
    observed HA and RA are related simply through the Earth rotation
    angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    means the position that would be seen by a perfect equatorial
    with its polar axis aligned to the Earth's axis of rotation.

12. It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

### Called ###

- [`apco13`](@ref): astrometry parameters, ICRS-observed, 2013
- [`atciq`](@ref): quick ICRS to CIRS
- [`atioq`](@ref): quick CIRS to observed

"""
function atco13(rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)
    aob = Ref{Cdouble}()
    zob = Ref{Cdouble}()
    hob = Ref{Cdouble}()
    dob = Ref{Cdouble}()
    rob = Ref{Cdouble}()
    eo = Ref{Cdouble}()
    i = ccall((:eraAtco13, liberfa), Cint,
              (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
               Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
               Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble},
               Ref{Cdouble}),
              rc, dc, pr, pd, px, rv, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh,
              wl, aob, zob, hob, dob, rob, eo)
    if i == -1
        throw(ERFAException("unacceptable date"))
    elseif i == +1
        @warn "dubious year"
    end
    return aob[], zob[], hob[], dob[], rob[], eo[]
end

"""
    atic13(ri, di, date1, date2)

Transform star RA,Dec from geocentric CIRS to ICRS astrometric.

### Given ###

- `ri`, `di`: CIRS geocentric RA,Dec (radians)
- `date1`: TDB as a 2-part...
- `date2`: ...Julian Date (Note 1)

### Returned ###

- `rc`, `dc`: ICRS astrometric RA,Dec (radians)
- `eo`: Equation of the origins (ERA-GST, Note 4)

### Notes ###

1. The TDB date date1+date2 is a Julian Date, apportioned in any
   convenient way between the two arguments.  For example,
   JD(TDB)=2450123.7 could be expressed in any of these ways, among
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
   good compromises between resolution and convenience.  For most
   applications of this function the choice will not be at all
   critical.

   TT can be used instead of TDB without any significant impact on
   accuracy.

2. Iterative techniques are used for the aberration and light
   deflection corrections so that the functions [`atic13`](@ref) (or
   [`aticq`](@ref)) and [`atci13`](@ref) (or [`atciq`](@ref)) are accurate inverses;
   even at the edge of the Sun's disk the discrepancy is only about
   1 nanoarcsecond.

3. The available accuracy is better than 1 milliarcsecond, limited
   mainly by the precession-nutation model that is used, namely
   IAU 2000A/2006.  Very close to solar system bodies, additional
   errors of up to several milliarcseconds can occur because of
   unmodeled light deflection;  however, the Sun's contribution is
   taken into account, to first order.  The accuracy limitations of
   the ERFA function [`epv00`](@ref) (used to compute Earth position and
   velocity) can contribute aberration errors of up to
   5 microarcseconds.  Light deflection at the Sun's limb is
   uncertain at the 0.4 mas level.

4. Should the transformation to (equinox based) J2000.0 mean place
   be required rather than (CIO based) ICRS coordinates, subtract the
   equation of the origins from the returned right ascension:
   RA = RI - EO.  (The [`anp`](@ref) function can then be applied, as
   required, to keep the result in the conventional 0-2pi range.)

### Called ###

- [`apci13`](@ref): astrometry parameters, ICRS-CIRS, 2013
- [`aticq`](@ref): quick CIRS to ICRS astrometric

"""
function atic13(ri, di, date1, date2)
    rc = Ref{Cdouble}()
    dc = Ref{Cdouble}()
    eo = Ref{Cdouble}()
    ccall((:eraAtic13, liberfa), Cvoid,
          (Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
          ri, di, date1, date2, rc, dc, eo)
    return rc[], dc[], eo[]
end

"""
    aticq(ri, di, astrom)

Quick CIRS RA,Dec to ICRS astrometric place, given the star-
independent astrometry parameters.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.
The star-independent astrometry parameters can be obtained by
calling one of the functions `apci[13]`, `apcg[13]`, `apco[13]`
or `apcs[13]`.

### Given ###

- `ri`, `di`: CIRS RA,Dec (radians)
- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: Longitude + s' (radians)
    - `xp1`: Polar motion xp wrt local meridian (radians)
    - `yp1`: Polar motion yp wrt local meridian (radians)
    - `sphi`: Sine of geodetic latitude
    - `cphi`: Cosine of geodetic latitude
    - `diurab`: Magnitude of diurnal aberration vector
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: Refraction constant A (radians)
    - `refb`: Refraction constant B (radians)

### Returned ###

- `rc`, `dc`: ICRS astrometric RA,Dec (radians)

### Notes ###

1. Only the Sun is taken into account in the light deflection
   correction.

2. Iterative techniques are used for the aberration and light
   deflection corrections so that the functions [`atic13`](@ref) (or
   [`aticq`](@ref)) and [`atci13`](@ref) (or [`atciq`](@ref)) are accurate inverses;
   even at the edge of the Sun's disk the discrepancy is only about
   1 nanoarcsecond.

### Called ###

- [`s2c`](@ref): spherical coordinates to unit vector
- [`trxp`](@ref): product of transpose of r-matrix and p-vector
- [`zp`](@ref): zero p-vector
- [`ab`](@ref): stellar aberration
- [`ldsun`](@ref): light deflection by the Sun
- [`c2s`](@ref): p-vector to spherical
- [`anp`](@ref): normalize angle into range +/- pi

"""
function aticq(ri, di, astrom)
    rc = Ref{Cdouble}()
    dc = Ref{Cdouble}()
    ccall((:eraAticq, liberfa), Cvoid,
          (Cdouble, Cdouble, Ref{ASTROM}, Ref{Cdouble}, Ref{Cdouble}),
          ri, di, astrom, rc, dc)
    return rc[], dc[]
end

"""
    aticqn(ri, di, astrom, b::Array{LDBODY})

Quick CIRS to ICRS astrometric place transformation, given the star-
independent astrometry parameters plus a list of light-deflecting
bodies.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.
The star-independent astrometry parameters can be obtained by
calling one of the functions `apci[13]`, `apcg[13]`, `apco[13]`
or `apcs[13]`.

### Given ###

- `ri`, `di`: CIRS RA,Dec (radians)
- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: Longitude + s' (radians)
    - `xp1`: Polar motion xp wrt local meridian (radians)
    - `yp1`: Polar motion yp wrt local meridian (radians)
    - `sphi`: Sine of geodetic latitude
    - `cphi`: Cosine of geodetic latitude
    - `diurab`: Magnitude of diurnal aberration vector
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: Refraction constant A (radians)
    - `refb`: Refraction constant B (radians)
- `n`: Number of bodies (Note 3)
- `b::Vector{LDBODY}`: Data for each of the n bodies (Notes 3,4):
    - `bm`: Mass of the body (solar masses, Note 5)
    - `dl`: Deflection limiter (Note 6)
    - `pv`: Barycentric PV of the body (au, au/day)

### Returned ###

- `rc`, `dc`: ICRS astrometric RA,Dec (radians)

### Notes ###

1. Iterative techniques are used for the aberration and light
   deflection corrections so that the functions [`aticqn`](@ref) and
   [`atciqn`](@ref) are accurate inverses; even at the edge of the Sun's
   disk the discrepancy is only about 1 nanoarcsecond.

2. If the only light-deflecting body to be taken into account is the
   Sun, the [`aticq`](@ref) function can be used instead.

3. The struct b contains n entries, one for each body to be
   considered.  If n = 0, no gravitational light deflection will be
   applied, not even for the Sun.

4. The struct b should include an entry for the Sun as well as for
   any planet or other body to be taken into account.  The entries
   should be in the order in which the light passes the body.

5. In the entry in the b struct for body i, the mass parameter
   b[i].bm can, as required, be adjusted in order to allow for such
   effects as quadrupole field.

6. The deflection limiter parameter b[i].dl is phi^2/2, where phi is
   the angular separation (in radians) between star and body at
   which limiting is applied.  As phi shrinks below the chosen
   threshold, the deflection is artificially reduced, reaching zero
   for phi = 0.   Example values suitable for a terrestrial
   observer, together with masses, are as follows:

   | `body i` | `b[i].bm`  | `b[i].dl` |
   |:---------|:-----------|:----------|
   |  Sun     | 1.0        |  6e-6     |
   |  Jupiter | 0.00095435 |  3e-9     |
   |  Saturn  | 0.00028574 |  3e-10    |

7. For efficiency, validation of the contents of the b array is
   omitted.  The supplied masses must be greater than zero, the
   position and velocity vectors must be right, and the deflection
   limiter greater than zero.

### Called ###

- [`s2c`](@ref): spherical coordinates to unit vector
- [`trxp`](@ref): product of transpose of r-matrix and p-vector
- [`zp`](@ref): zero p-vector
- [`ab`](@ref): stellar aberration
- [`ldn`](@ref): light deflection by n bodies
- [`c2s`](@ref): p-vector to spherical
- [`anp`](@ref): normalize angle into range +/- pi

"""
function aticqn(ri, di, astrom, b::Array{LDBODY})
    rc = Ref{Cdouble}()
    dc = Ref{Cdouble}()
    n = length(b)
    ccall((:eraAticqn, liberfa), Cvoid,
          (Cdouble, Cdouble, Ref{ASTROM}, Cint, Ref{LDBODY}, Ref{Cdouble}, Ref{Cdouble}),
          ri, di, astrom, n, b, rc, dc)
    return rc[], dc[]
end

"""
    atio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)

CIRS RA,Dec to observed place.  The caller supplies UTC, site
coordinates, ambient air conditions and observing wavelength.

### Given ###

- `ri`: CIRS right ascension (CIO-based, radians)
- `di`: CIRS declination (radians)
- `utc1`: UTC as a 2-part...
- `utc2`: ...quasi Julian Date (Notes 1,2)
- `dut1`: UT1-UTC (seconds, Note 3)
- `elong`: Longitude (radians, east +ve, Note 4)
- `phi`: Geodetic latitude (radians, Note 4)
- `hm`: Height above ellipsoid (m, geodetic Notes 4,6)
- `xp`, `yp`: Polar motion coordinates (radians, Note 5)
- `phpa`: Pressure at the observer (hPa = mB, Note 6)
- `tc`: Ambient temperature at the observer (deg C)
- `rh`: Relative humidity at the observer (range 0-1)
- `wl`: Wavelength (micrometers, Note 7)

### Returned ###

- `aob`: Observed azimuth (radians: N=0,E=90)
- `zob`: Observed zenith distance (radians)
- `hob`: Observed hour angle (radians)
- `dob`: Observed declination (radians)
- `rob`: Observed right ascension (CIO-based, radians)

### Notes ###

1.  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function [`dtf2d`](@ref) to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

2.  The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the
    future to be trusted.  See [`dat`](@ref) for further details.

3.  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

4.  The geographical coordinates are with respect to the [`WGS84`](@ref)
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

5.  The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many
    applications, xp and yp can be set to zero.

6.  If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB), is
    available, an adequate estimate of hm can be obtained from the
    expression

    ```
    hm = -29.3 * tsl * log ( phpa / 1013.25 );
    ```

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    ```
    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    ```

    Note, however, that the refraction is nearly proportional to
    the pressure and that an accurate phpa value is important for
    precise work.

7.  The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

8.  "Observed" Az,ZD means the position that would be seen by a
    perfect geodetically aligned theodolite.  (Zenith distance is
    used rather than altitude in order to reflect the fact that no
    allowance is made for depression of the horizon.)  This is
    related to the observed HA,Dec via the standard rotation, using
    the geodetic latitude (corrected for polar motion), while the
    observed HA and RA are related simply through the Earth rotation
    angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    means the position that would be seen by a perfect equatorial
    with its polar axis aligned to the Earth's axis of rotation.

9.  The accuracy of the result is limited by the corrections for
    refraction, which use a simple ``A*tan(z) + B*tan^3(z)`` model.
    Providing the meteorological parameters are known accurately and
    there are no gross local effects, the predicted astrometric
    coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    (radio) for a zenith distance of less than 70 degrees, better
    than 30 arcsec (optical or radio) at 85 degrees and better
    than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

10. The complementary functions [`atio13`](@ref) and [`atoi13`](@ref) are self-
    consistent to better than 1 microarcsecond all over the
    celestial sphere.

11. It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

### Called ###

- [`apio13`](@ref): astrometry parameters, CIRS-observed, 2013
- [`atioq`](@ref): quick CIRS to observed

"""
function atio13(ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)
    aob = Ref{Cdouble}()
    zob = Ref{Cdouble}()
    hob = Ref{Cdouble}()
    dob = Ref{Cdouble}()
    rob = Ref{Cdouble}()
    i = ccall((:eraAtio13, liberfa), Cint,
              (Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
               Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble},
               Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}),
              ri, di, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl, aob,
              zob, hob, dob, rob)
    if i == -1
        throw(ERFAException("unacceptable date"))
    elseif i == +1
        @warn "dubious year"
    end
    return aob[], zob[], hob[], dob[], rob[]
end

"""
    atioq(ri, di, astrom)

Quick CIRS to observed place transformation.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.
The star-independent astrometry parameters can be obtained by
calling `apio[13]` or `apco[13]`.

### Given ###

- `ri`: CIRS right ascension
- `di`: CIRS declination
- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: Longitude + s' (radians)
    - `xp1`: Polar motion xp wrt local meridian (radians)
    - `yp1`: Polar motion yp wrt local meridian (radians)
    - `sphi`: Sine of geodetic latitude
    - `cphi`: Cosine of geodetic latitude
    - `diurab`: Magnitude of diurnal aberration vector
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: Refraction constant A (radians)
    - `refb`: Refraction constant B (radians)

### Returned ###

- `aob`: Observed azimuth (radians: N=0,E=90)
- `zob`: Observed zenith distance (radians)
- `hob`: Observed hour angle (radians)
- `dob`: Observed declination (radians)
- `rob`: Observed right ascension (CIO-based, radians)

### Notes ###

1. This function returns zenith distance rather than altitude in
   order to reflect the fact that no allowance is made for
   depression of the horizon.

2. The accuracy of the result is limited by the corrections for
   refraction, which use a simple ``A*tan(z) + B*tan^3(z)`` model.
   Providing the meteorological parameters are known accurately and
   there are no gross local effects, the predicted observed
   coordinates should be within 0.05 arcsec (optical) or 1 arcsec
   (radio) for a zenith distance of less than 70 degrees, better
   than 30 arcsec (optical or radio) at 85 degrees and better
   than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

   Without refraction, the complementary functions [`atioq`](@ref) and
   [`atoiq`](@ref) are self-consistent to better than 1 microarcsecond all
   over the celestial sphere.  With refraction included, consistency
   falls off at high zenith distances, but is still better than
   0.05 arcsec at 85 degrees.

3. It is advisable to take great care with units, as even unlikely
   values of the input parameters are accepted and processed in
   accordance with the models used.

4. The CIRS RA,Dec is obtained from a star catalog mean place by
   allowing for space motion, parallax, the Sun's gravitational lens
   effect, annual aberration and precession-nutation.  For star
   positions in the ICRS, these effects can be applied by means of
   the [`atci13`](@ref) (etc.) functions.  Starting from classical "mean
   place" systems, additional transformations will be needed first.

5. "Observed" Az,El means the position that would be seen by a
   perfect geodetically aligned theodolite.  This is obtained from
   the CIRS RA,Dec by allowing for Earth orientation and diurnal
   aberration, rotating from equator to horizon coordinates, and
   then adjusting for refraction.  The HA,Dec is obtained by
   rotating back into equatorial coordinates, and is the position
   that would be seen by a perfect equatorial with its polar axis
   aligned to the Earth's axis of rotation.  Finally, the RA is
   obtained by subtracting the HA from the local ERA.

6. The star-independent CIRS-to-observed-place parameters in ASTROM
   may be computed with `apio[13]` or `apco[13]`.  If nothing has
   changed significantly except the time, `aper[13]` may be used to
   perform the requisite adjustment to the astrom structure.

### Called ###

- [`s2c`](@ref): spherical coordinates to unit vector
- [`c2s`](@ref): p-vector to spherical
- [`anp`](@ref): normalize angle into range 0 to 2pi

"""
function atioq(ri, di, astrom)
    aob = Ref{Cdouble}()
    zob = Ref{Cdouble}()
    hob = Ref{Cdouble}()
    dob = Ref{Cdouble}()
    rob = Ref{Cdouble}()
    ccall((:eraAtioq, liberfa), Cvoid,
          (Cdouble, Cdouble, Ref{ASTROM}, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble},
           Ref{Cdouble}, Ref{Cdouble}),
          ri, di, astrom, aob, zob, hob, dob, rob)
    return aob[], zob[], hob[], dob[], rob[]
end

"""
    atoc13(typeofcoordinates, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)

Observed place at a groundbased site to to ICRS astrometric RA,Dec.
The caller supplies UTC, site coordinates, ambient air conditions
and observing wavelength.

### Given ###

- `type`: Type of coordinates - "R", "H" or "A" (Notes 1,2)
- `ob1`: Observed Az, HA or RA (radians; Az is N=0,E=90)
- `ob2`: Observed ZD or Dec (radians)
- `utc1`: UTC as a 2-part...
- `utc2`: ...quasi Julian Date (Notes 3,4)
- `dut1`: UT1-UTC (seconds, Note 5)
- `elong`: Longitude (radians, east +ve, Note 6)
- `phi`: Geodetic latitude (radians, Note 6)
- `hm`: Height above ellipsoid (m, geodetic Notes 6,8)
- `xp`, `yp`: Polar motion coordinates (radians, Note 7)
- `phpa`: Pressure at the observer (hPa = mB, Note 8)
- `tc`: Ambient temperature at the observer (deg C)
- `rh`: Relative humidity at the observer (range 0-1)
- `wl`: Wavelength (micrometers, Note 9)

### Returned ###

- `rc`, `dc`: ICRS astrometric RA,Dec (radians)

### Notes ###

1.  "Observed" Az,ZD means the position that would be seen by a
    perfect geodetically aligned theodolite.  (Zenith distance is
    used rather than altitude in order to reflect the fact that no
    allowance is made for depression of the horizon.)  This is
    related to the observed HA,Dec via the standard rotation, using
    the geodetic latitude (corrected for polar motion), while the
    observed HA and RA are related simply through the Earth rotation
    angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    means the position that would be seen by a perfect equatorial
    with its polar axis aligned to the Earth's axis of rotation.

2.  Only the first character of the type argument is significant.
    "R" or "r" indicates that ob1 and ob2 are the observed right
    ascension and declination;  "H" or "h" indicates that they are
    hour angle (west +ve) and declination;  anything else ("A" or
    "a" is recommended) indicates that ob1 and ob2 are azimuth
    (north zero, east 90 deg) and zenith distance.

3.  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function [`dtf2d`](@ref) to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

4.  The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the
    future to be trusted.  See [`dat`](@ref) for further details.

5.  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

6.  The geographical coordinates are with respect to the [`WGS84`](@ref)
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

7.  The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many
    applications, xp and yp can be set to zero.

8.  If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB), is
    available, an adequate estimate of hm can be obtained from the
    expression

    ```
    hm = -29.3 * tsl * log ( phpa / 1013.25 );
    ```

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    ```
    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    ```

    Note, however, that the refraction is nearly proportional to
    the pressure and that an accurate phpa value is important for
    precise work.

9.  The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

10. The accuracy of the result is limited by the corrections for
    refraction, which use a simple ``A*tan(z) + B*tan^3(z)`` model.
    Providing the meteorological parameters are known accurately and
    there are no gross local effects, the predicted astrometric
    coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    (radio) for a zenith distance of less than 70 degrees, better
    than 30 arcsec (optical or radio) at 85 degrees and better
    than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

    Without refraction, the complementary functions [`atco13`](@ref) and
    [`atoc13`](@ref) are self-consistent to better than 1 microarcsecond
    all over the celestial sphere.  With refraction included,
    consistency falls off at high zenith distances, but is still
    better than 0.05 arcsec at 85 degrees.

11. It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

### Called ###

- [`apco13`](@ref): astrometry parameters, ICRS-observed
- [`atoiq`](@ref): quick observed to CIRS
- [`aticq`](@ref): quick CIRS to ICRS

"""
function atoc13(typeofcoordinates, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)
    rc = Ref{Cdouble}()
    dc = Ref{Cdouble}()
    if !(typeofcoordinates in ("R", "r", "H", "h", "A", "a"))
        typeofcoordinates = "A"
    end
    i = ccall((:eraAtoc13, liberfa), Cint,
              (Cstring, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
               Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble},
               Ref{Cdouble}),
              typeofcoordinates, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa,
              tk, rh, wl, rc, dc)
    if i == -1
        throw(ERFAException("unacceptable date"))
    elseif i == +1
        @warn "dubious year"
    end
    rc[], dc[]
end

"""
    atoi13(typeofcoordinates, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)

Observed place to CIRS.  The caller supplies UTC, site coordinates,
ambient air conditions and observing wavelength.

### Given ###

- `type`: Type of coordinates - "R", "H" or "A" (Notes 1,2)
- `ob1`: Observed Az, HA or RA (radians; Az is N=0,E=90)
- `ob2`: Observed ZD or Dec (radians)
- `utc1`: UTC as a 2-part...
- `utc2`: ...quasi Julian Date (Notes 3,4)
- `dut1`: UT1-UTC (seconds, Note 5)
- `elong`: Longitude (radians, east +ve, Note 6)
- `phi`: Geodetic latitude (radians, Note 6)
- `hm`: Height above the ellipsoid (meters, Notes 6,8)
- `xp`, `yp`: Polar motion coordinates (radians, Note 7)
- `phpa`: Pressure at the observer (hPa = mB, Note 8)
- `tc`: Ambient temperature at the observer (deg C)
- `rh`: Relative humidity at the observer (range 0-1)
- `wl`: Wavelength (micrometers, Note 9)

### Returned ###

- `ri`: CIRS right ascension (CIO-based, radians)
- `di`: CIRS declination (radians)

### Notes ###

1.  "Observed" Az,ZD means the position that would be seen by a
    perfect geodetically aligned theodolite.  (Zenith distance is
    used rather than altitude in order to reflect the fact that no
    allowance is made for depression of the horizon.)  This is
    related to the observed HA,Dec via the standard rotation, using
    the geodetic latitude (corrected for polar motion), while the
    observed HA and RA are related simply through the Earth rotation
    angle and the site longitude.  "Observed" RA,Dec or HA,Dec thus
    means the position that would be seen by a perfect equatorial
    with its polar axis aligned to the Earth's axis of rotation.

2.  Only the first character of the type argument is significant.
    "R" or "r" indicates that ob1 and ob2 are the observed right
    ascension and declination;  "H" or "h" indicates that they are
    hour angle (west +ve) and declination;  anything else ("A" or
    "a" is recommended) indicates that ob1 and ob2 are azimuth
    (north zero, east 90 deg) and zenith distance.

3.  utc1+utc2 is quasi Julian Date (see Note 2), apportioned in any
    convenient way between the two arguments, for example where utc1
    is the Julian Day Number and utc2 is the fraction of a day.

    However, JD cannot unambiguously represent UTC during a leap
    second unless special measures are taken.  The convention in the
    present function is that the JD day represents UTC days whether
    the length is 86399, 86400 or 86401 SI seconds.

    Applications should use the function [`dtf2d`](@ref) to convert from
    calendar date and time of day into 2-part quasi Julian Date, as
    it implements the leap-second-ambiguity convention just
    described.

4.  The warning status "dubious year" flags UTCs that predate the
    introduction of the time scale or that are too far in the
    future to be trusted.  See [`dat`](@ref) for further details.

5.  UT1-UTC is tabulated in IERS bulletins.  It increases by exactly
    one second at the end of each positive UTC leap second,
    introduced in order to keep UT1-UTC within +/- 0.9s.  n.b. This
    practice is under review, and in the future UT1-UTC may grow
    essentially without limit.

6.  The geographical coordinates are with respect to the [`WGS84`](@ref)
    reference ellipsoid.  TAKE CARE WITH THE LONGITUDE SIGN:  the
    longitude required by the present function is east-positive
    (i.e. right-handed), in accordance with geographical convention.

7.  The polar motion xp,yp can be obtained from IERS bulletins.  The
    values are the coordinates (in radians) of the Celestial
    Intermediate Pole with respect to the International Terrestrial
    Reference System (see IERS Conventions 2003), measured along the
    meridians 0 and 90 deg west respectively.  For many
    applications, xp and yp can be set to zero.

8.  If hm, the height above the ellipsoid of the observing station
    in meters, is not known but phpa, the pressure in hPa (=mB), is
    available, an adequate estimate of hm can be obtained from the
    expression

    ```
    hm = -29.3 * tsl * log ( phpa / 1013.25 );
    ```

    where tsl is the approximate sea-level air temperature in K
    (See Astrophysical Quantities, C.W.Allen, 3rd edition, section
    52).  Similarly, if the pressure phpa is not known, it can be
    estimated from the height of the observing station, hm, as
    follows:

    ```
    phpa = 1013.25 * exp ( -hm / ( 29.3 * tsl ) );
    ```

    Note, however, that the refraction is nearly proportional to
    the pressure and that an accurate phpa value is important for
    precise work.

9.  The argument wl specifies the observing wavelength in
    micrometers.  The transition from optical to radio is assumed to
    occur at 100 micrometers (about 3000 GHz).

10. The accuracy of the result is limited by the corrections for
    refraction, which use a simple ``A*tan(z) + B*tan^3(z)`` model.
    Providing the meteorological parameters are known accurately and
    there are no gross local effects, the predicted astrometric
    coordinates should be within 0.05 arcsec (optical) or 1 arcsec
    (radio) for a zenith distance of less than 70 degrees, better
    than 30 arcsec (optical or radio) at 85 degrees and better
    than 20 arcmin (optical) or 30 arcmin (radio) at the horizon.

    Without refraction, the complementary functions [`atio13`](@ref) and
    [`atoi13`](@ref) are self-consistent to better than 1 microarcsecond
    all over the celestial sphere.  With refraction included,
    consistency falls off at high zenith distances, but is still
    better than 0.05 arcsec at 85 degrees.

12. It is advisable to take great care with units, as even unlikely
    values of the input parameters are accepted and processed in
    accordance with the models used.

### Called ###

- [`apio13`](@ref): astrometry parameters, CIRS-observed, 2013
- [`atoiq`](@ref): quick observed to CIRS

"""
function atoi13(typeofcoordinates, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa, tk, rh, wl)
    ri = Ref{Cdouble}()
    di = Ref{Cdouble}()
    i = ccall((:eraAtoi13, liberfa), Cint,
              (Cstring, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble,
               Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Ref{Cdouble},
               Ref{Cdouble}),
              typeofcoordinates, ob1, ob2, utc1, utc2, dut1, elong, phi, hm, xp, yp, phpa,
              tk, rh, wl, ri, di)
    if i == -1
        throw(ERFAException("unacceptable date"))
    elseif i == +1
        @warn "dubious year"
    end
    return ri[], di[]
end

"""
    atoiq(typeofcoordinates, ob1, ob2, astrom)

Quick observed place to CIRS, given the star-independent astrometry
parameters.

Use of this function is appropriate when efficiency is important and
where many star positions are all to be transformed for one date.
The star-independent astrometry parameters can be obtained by
calling `apio[13]` or `apco[13]`.

### Given ###

- `type`: Type of coordinates: "R", "H" or "A" (Note 1)
- `ob1`: Observed Az, HA or RA (radians; Az is N=0,E=90)
- `ob2`: Observed ZD or Dec (radians)
- `astrom`: Star-independent astrometry parameters:
    - `pmt`: PM time interval (SSB, Julian years)
    - `eb`: SSB to observer (vector, au)
    - `eh`: Sun to observer (unit vector)
    - `em`: Distance from Sun to observer (au)
    - `v`: Barycentric observer velocity (vector, c)
    - `bm1`: ``\\sqrt{1-|v|^2}`` Reciprocal of Lorenz factor
    - `bpn`: Bias-precession-nutation matrix
    - `along`: Longitude + s' (radians)
    - `xp1`: Polar motion xp wrt local meridian (radians)
    - `yp1`: Polar motion yp wrt local meridian (radians)
    - `sphi`: Sine of geodetic latitude
    - `cphi`: Cosine of geodetic latitude
    - `diurab`: Magnitude of diurnal aberration vector
    - `l`: "Local" Earth rotation angle (radians)
    - `refa`: Refraction constant A (radians)
    - `refb`: Refraction constant B (radians)

### Returned ###

- `ri`: CIRS right ascension (CIO-based, radians)
- `di`: CIRS declination (radians)

### Notes ###

1. "Observed" Az,El means the position that would be seen by a
   perfect geodetically aligned theodolite.  This is related to
   the observed HA,Dec via the standard rotation, using the geodetic
   latitude (corrected for polar motion), while the observed HA and
   RA are related simply through the Earth rotation angle and the
   site longitude.  "Observed" RA,Dec or HA,Dec thus means the
   position that would be seen by a perfect equatorial with its
   polar axis aligned to the Earth's axis of rotation.  By removing
   from the observed place the effects of atmospheric refraction and
   diurnal aberration, the CIRS RA,Dec is obtained.

2. Only the first character of the type argument is significant.
   "R" or "r" indicates that ob1 and ob2 are the observed right
   ascension and declination;  "H" or "h" indicates that they are
   hour angle (west +ve) and declination;  anything else ("A" or
   "a" is recommended) indicates that ob1 and ob2 are azimuth (north
   zero, east 90 deg) and zenith distance.  (Zenith distance is used
   rather than altitude in order to reflect the fact that no
   allowance is made for depression of the horizon.)

3. The accuracy of the result is limited by the corrections for
   refraction, which use a simple ``A*tan(z) + B*tan^3(z)`` model.
   Providing the meteorological parameters are known accurately and
   there are no gross local effects, the predicted observed
   coordinates should be within 0.05 arcsec (optical) or 1 arcsec
   (radio) for a zenith distance of less than 70 degrees, better
   than 30 arcsec (optical or radio) at 85 degrees and better than
   20 arcmin (optical) or 30 arcmin (radio) at the horizon.

   Without refraction, the complementary functions [`atioq`](@ref) and
   [`atoiq`](@ref) are self-consistent to better than 1 microarcsecond all
   over the celestial sphere.  With refraction included, consistency
   falls off at high zenith distances, but is still better than
   0.05 arcsec at 85 degrees.

4. It is advisable to take great care with units, as even unlikely
   values of the input parameters are accepted and processed in
   accordance with the models used.

### Called ###

- [`s2c`](@ref): spherical coordinates to unit vector
- [`c2s`](@ref): p-vector to spherical
- [`anp`](@ref): normalize angle into range 0 to 2pi

"""
function atoiq(typeofcoordinates, ob1, ob2, astrom)
    ri = Ref{Cdouble}()
    di = Ref{Cdouble}()
    ccall((:eraAtoiq, liberfa),
          Cvoid, (Cstring, Cdouble, Cdouble, Ref{ASTROM}, Ref{Cdouble}, Ref{Cdouble}),
          typeofcoordinates, ob1, ob2, astrom, ri, di)
    return ri[], di[]
end

"""
    anp(a)

Normalize angle into the range 0 <= a < 2pi.

!!! warning "Deprecated"
    Use `Base.mod2pi` instead.

### Given ###

- `a`: Angle (radians)

### Returned ###

- Angle in range 0-2pi

"""
anp

@deprecate anp mod2pi
_anp(a) = ccall((:eraAnp, liberfa), Cdouble, (Cdouble,), a)

"""
    anpm(a)

Normalize angle into the range -pi <= a < +pi.

### Given ###

- `a`: Angle (radians)

### Returned ###

- Angle in range +/-pi

"""
anpm(a) = ccall((:eraAnpm, liberfa), Cdouble, (Cdouble,), a)

"""
    a2af(ndp, a)

Decompose radians into degrees, arcminutes, arcseconds, fraction.

### Given ###

- `ndp`: Resolution (Note 1)
- `angle`: Angle in radians

### Returned ###

- `sign`: '+' or '-'
- `idmsf`: Degrees, arcminutes, arcseconds, fraction

### Called ###

- [`d2tf`](@ref): decompose days to hms

### Notes ###

1. The argument ndp is interpreted as follows:

   | ndp |     resolution       |
   |:----|:---------------------|
   |  :  | ...0000 00 00        |
   | -7  |    1000 00 00        |
   | -6  |     100 00 00        |
   | -5  |      10 00 00        |
   | -4  |       1 00 00        |
   | -3  |       0 10 00        |
   | -2  |       0 01 00        |
   | -1  |       0 00 10        |
   |  0  |       0 00 01        |
   |  1  |       0 00 00.1      |
   |  2  |       0 00 00.01     |
   |  3  |       0 00 00.001    |
   |  :  |       0 00 00.000... |

2. The largest positive useful value for ndp is determined by the
   size of angle, the format of doubles on the target platform, and
   the risk of overflowing idmsf[3].  On a typical platform, for
   angle up to 2pi, the available floating-point precision might
   correspond to ndp=12.  However, the practical limit is typically
   ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
   only 16 bits.

3. The absolute value of angle may exceed 2pi.  In cases where it
   does not, it is up to the caller to test for and handle the
   case where angle is very nearly 2pi and rounds up to 360 degrees,
   by testing for idmsf[0]=360 and setting idmsf[0-3] to zero.

"""
a2af

"""
    a2tf(ndp, a)

Decompose radians into hours, minutes, seconds, fraction.

### Given ###

- `ndp`: Resolution (Note 1)
- `angle`: Angle in radians

### Returned ###

- `sign`: '+' or '-'
- `ihmsf`: Hours, minutes, seconds, fraction

### Called ###

- [`d2tf`](@ref): decompose days to hms

### Notes ###

1. The argument ndp is interpreted as follows:

   | ndp |     resolution       |
   |:----|:---------------------|
   |  :  | ...0000 00 00        |
   | -7  |    1000 00 00        |
   | -6  |     100 00 00        |
   | -5  |      10 00 00        |
   | -4  |       1 00 00        |
   | -3  |       0 10 00        |
   | -2  |       0 01 00        |
   | -1  |       0 00 10        |
   |  0  |       0 00 01        |
   |  1  |       0 00 00.1      |
   |  2  |       0 00 00.01     |
   |  3  |       0 00 00.001    |
   |  :  |       0 00 00.000... |

2. The largest positive useful value for ndp is determined by the
   size of angle, the format of doubles on the target platform, and
   the risk of overflowing ihmsf[3].  On a typical platform, for
   angle up to 2pi, the available floating-point precision might
   correspond to ndp=12.  However, the practical limit is typically
   ndp=9, set by the capacity of a 32-bit int, or ndp=4 if int is
   only 16 bits.

3. The absolute value of angle may exceed 2pi.  In cases where it
   does not, it is up to the caller to test for and handle the
   case where angle is very nearly 2pi and rounds up to 24 hours,
   by testing for ihmsf[0]=24 and setting ihmsf[0-3] to zero.

"""
a2tf

for name in ("a2af",
             "a2tf")
    f = Symbol(name)
    fc = "era" * uppercasefirst(name)
    @eval begin
        function ($f)(ndp, a)
            s = Ref{Cchar}('+')
            i = zeros(Cint, 4)
            ccall(($fc, liberfa), Cvoid,
                  (Cint, Cdouble, Ref{Cchar}, Ptr{Cint}),
                  ndp, a, s, i)
            return Char(s[]), i[1], i[2], i[3], i[4]
        end
    end
end

"""
    af2a(s, ideg, iamin, asec)

Convert degrees, arcminutes, arcseconds to radians.

### Given ###

- `s`: Sign:  '-' = negative, otherwise positive
- `ideg`: Degrees
- `iamin`: Arcminutes
- `asec`: Arcseconds

### Returned ###

- `rad`: Angle in radians

### Notes ###

1.  The result is computed even if any of the range checks fail.

2.  Negative ideg, iamin and/or asec produce a warning status, but
    the absolute value is used in the conversion.

3.  If there are multiple errors, the status value reflects only the
    first, the smallest taking precedence.

"""
function af2a(s, ideg, iamin, asec)
    rad = Ref{Cdouble}()
    i = ccall((:eraAf2a, liberfa), Cint,
                (Cchar, Cint, Cint, Cdouble, Ref{Cdouble}),
                s, ideg, iamin, asec, rad)
    if i == 1
        throw(ERFAException("ideg outside range 0-359"))
    elseif i == 2
        throw(ERFAException("iamin outside range 0-59"))
    elseif i == 3
        throw(ERFAException("asec outside range 0-59.999..."))
    end
    return rad[]
end

