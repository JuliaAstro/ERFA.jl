"""
    icrs2g(dr, dd)

Transformation from ICRS to Galactic Coordinates.

### Given ###

- `dr`: ICRS right ascension (radians)
- `dd`: ICRS declination (radians)

### Returned ###

- `dl`: Galactic longitude (radians)
- `db`: Galactic latitude (radians)

### Notes ###

1. The IAU 1958 system of Galactic coordinates was defined with
   respect to the now obsolete reference system FK4 B1950.0.  When
   interpreting the system in a modern context, several factors have
   to be taken into account:

   - The inclusion in FK4 positions of the E-terms of aberration.

   - The distortion of the FK4 proper motion system by differential
     Galactic rotation.

   - The use of the B1950.0 equinox rather than the now-standard
     J2000.0.

   - The frame bias between ICRS and the J2000.0 mean place system.

   The Hipparcos Catalogue (Perryman & ESA 1997) provides a rotation
   matrix that transforms directly between ICRS and Galactic
   coordinates with the above factors taken into account.  The
   matrix is derived from three angles, namely the ICRS coordinates
   of the Galactic pole and the longitude of the ascending node of
   the galactic equator on the ICRS equator.  They are given in
   degrees to five decimal places and for canonical purposes are
   regarded as exact.  In the Hipparcos Catalogue the matrix
   elements are given to 10 decimal places (about 20 microarcsec).
   In the present ERFA function the matrix elements have been
   recomputed from the canonical three angles and are given to 30
   decimal places.

2. The inverse transformation is performed by the function [`g2icrs`](@ref).

### Called ###

- [`anp`](@ref): normalize angle into range 0 to 2pi
- [`anpm`](@ref): normalize angle into range +/- pi
- [`s2c`](@ref): spherical coordinates to unit vector
- [`rxp`](@ref): product of r-matrix and p-vector
- [`c2s`](@ref): p-vector to spherical

### Reference ###

- Perryman M.A.C. & ESA, 1997, ESA SP-1200, The Hipparcos and Tycho
    catalogues.  Astrometric and photometric star catalogues
    derived from the ESA Hipparcos Space Astrometry Mission.  ESA
    Publications Division, Noordwijk, Netherlands.

"""
function icrs2g(a, b)
    r1 = Ref{Cdouble}()
    r2 = Ref{Cdouble}()
    ccall((:eraIcrs2g, liberfa), Cvoid,
            (Cdouble, Cdouble, Ref{Cdouble}, Ref{Cdouble}),
            a, b, r1, r2)
    return r1[], r2[]
end

"""
    ir()

Initialize an r-matrix to the identity matrix.

!!! warning "Deprecated"
    Use `Array{Float64}(LinearAlgebra.I, 3, 3)` instead.

### Returned ###

- `r`: r-matrix
"""
ir

function _ir()
    r = Array{Cdouble}(undef, 3, 3)
    ccall((:eraIr, liberfa), Cvoid, (Ptr{Cdouble},), r)
    return r
end

@deprecate ir() Array{Float64}(I, 3, 3)

