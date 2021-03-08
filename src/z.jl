"""
    zp(p)

Zero a p-vector.

!!! warning "Deprecated"
    Use `fill!(p, 0.0)` instead.

### Returned ###

- `p`: zero p-vector

"""
zp

function _zp(p)
    ccall((:eraZp, liberfa), Cvoid, (Ptr{Cdouble},), p)
    return p
end

@deprecate zp(p) fill!(p, 0.0)

"""
    zpv(pv)

Zero a pv-vector.

!!! warning "Deprecated"
    Use `fill!.(pv, 0.0)` instead.

### Returned ###

- `p`: zero pv-vector

"""
zpv

function _zpv(pv)
    _pv = array_to_cmatrix(pv; n=3)
    ccall((:eraZpv, liberfa), Cvoid, (Ptr{Cdouble},), _pv)
    return cmatrix_to_array(_pv)
end

@deprecate zpv(pv) fill!.(pv, 0.0)

"""
    zr(r)

Initialize an r-matrix to the null matrix.

!!! warning "Deprecated"
    Use `fill!(r, 0.0)` instead.

### Returned ###

- `r`: r-matrix

"""
zr

function _zr(r)
    ccall((:eraZr, liberfa), Cvoid, (Ptr{Cdouble},), r)
    return r
end

@deprecate zr(r) fill!(r, 0.0)

