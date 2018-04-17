"""
    jd2cal(dr, dd)

Julian Date to Gregorian year, month, day, and fraction of a day.

### Given ###

* `dj1`, `dj2`: Julian Date (Notes 1, 2)

### Returned (arguments) ###

   iy        int      year
   im        int      month
   id        int      day
   fd        double   fraction of day

### Returned (function value) ###

             int      status:
                         0 = OK
                        -1 = unacceptable date (Note 1)

### Notes ###

1) The earliest valid date is -68569.5 (-4900 March 1).  The
   largest value accepted is 1e9.

2) The Julian Date is apportioned in any convenient way between
   the arguments dj1 and dj2.  For example, JD=2450123.7 could
   be expressed in any of these ways, among others:

          dj1             dj2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

3) In early eras the conversion is from the "proleptic Gregorian
   calendar";  no account is taken of the date(s) of adoption of
   the Gregorian calendar, nor is the AD/BC numbering convention
   observed.

### Reference ###

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 12.92 (p604).

"""
function jd2cal(d1, d2)
    iy = Ref{Cint}(0)
    imo = Ref{Cint}(0)
    id = Ref{Cint}(0)
    fd = Ref(0.0)
    i = ccall((:eraJd2cal, liberfa), Cint,
              (Cdouble, Cdouble, Ref{Cint}, Ref{Cint}, Ref{Cint}, Ref{Cdouble}),
              d1, d2, iy, imo, id, fd)
    @assert i == 0
    iy[], imo[], id[], fd[]
end

"""
    jdcalf(dr, dd)

Julian Date to Gregorian Calendar, expressed in a form convenient
for formatting messages:  rounded to a specified precision.

### Given ###

* `ndp`: Number of decimal places of days in fraction
* `dj1`, `dj2`: Dj1+dj2 = Julian Date (Note 1)

### Returned ###

* `iymdf`: Year, month, day, fraction in Gregorian
                      calendar

### Returned (function value) ###

             int      status:
                        -1 = date out of range
                         0 = OK
                        +1 = NDP not 0-9 (interpreted as 0)

### Notes ###

1) The Julian Date is apportioned in any convenient way between
   the arguments dj1 and dj2.  For example, JD=2450123.7 could
   be expressed in any of these ways, among others:

           dj1            dj2

       2450123.7           0.0       (JD method)
       2451545.0       -1421.3       (J2000 method)
       2400000.5       50123.2       (MJD method)
       2450123.5           0.2       (date & time method)

2) In early eras the conversion is from the "Proleptic Gregorian
   Calendar";  no account is taken of the date(s) of adoption of
   the Gregorian Calendar, nor is the AD/BC numbering convention
   observed.

3) Refer to the function eraJd2cal.

4) NDP should be 4 or less if internal overflows are to be
   avoided on machines which use 16-bit integers.

### Called ###

* `eraJd2cal`: JD to Gregorian calendar

### Reference ###

   Explanatory Supplement to the Astronomical Almanac,
   P. Kenneth Seidelmann (ed), University Science Books (1992),
   Section 12.92 (p604).

"""
function jdcalf(ndp, d1, d2)
    iymdf = Int32[0, 0, 0, 0]
    i = ccall((:eraJdcalf, liberfa), Cint,
              (Cint, Cdouble, Cdouble, Ptr{Cint}),
              ndp, d1, d2, iymdf)
    @assert i == 0
    iymdf[1], iymdf[2], iymdf[3], iymdf[4]
end
