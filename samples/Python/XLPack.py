# ****************************************
# *                                      *
# *  XLPack Numerical Library            *
# *  Version 7.0 (January 31, 2023)      *
# *  (C) 2014-2023  K Technologies       *
# *                                      *
# ****************************************/
import numpy as np
from ctypes import *
import struct
if struct.calcsize('P') == 8:
    dll = CDLL('XLPack.dll')
    dllsp = CDLL('XLPackSp.dll')
    dllpde = CDLL('XLPackPde.dll')
else:
    dll = CDLL('XLPack_32.dll')
    dllsp = CDLL('XLPackSp_32.dll')
    dllpde = CDLL('XLPackPde_32.dll')
errno = c_int(0)
def d1num(i):
    dll._d1num.restype = c_double
    dll._d1num.argtype = [c_int]
    return dll._d1num(c_int(i))
def factorial(x):
    dll.x_factorial.restype = c_double
    dll.x_factorial.argtype = [c_uint, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_factorial(c_uint(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def chebs(c, n, x):
    dll.x_chebs.restype = c_double
    dll.x_chebs.argtype = [POINTER(c_double), c_int, c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_chebs(c.ctypes.data_as(POINTER(c_double)), c_int(n), c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def sqrt1pm1(x):
    dll.x_sqrt1pm1.restype = c_double
    dll.x_sqrt1pm1.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_sqrt1pm1(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def powm1(x, y):
    dll.x_powm1.restype = c_double
    dll.x_powm1.argtype = [c_double, c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_powm1(c_double(x), c_double(y), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def cos_pi(x):
    dll.x_cos_pi.restype = c_double
    dll.x_cos_pi.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_cos_pi(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def sin_pi(x):
    dll.x_sin_pi.restype = c_double
    dll.x_sin_pi.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_sin_pi(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def li(x):
    dll.x_li.restype = c_double
    dll.x_li.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_li(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def ei(x):
    dll.x_ei.restype = c_double
    dll.x_ei.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_ei(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def e1(x):
    dll.x_e1.restype = c_double
    dll.x_e1.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_e1(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def digamma(x):
    dll.x_digamma.restype = c_double
    dll.x_digamma.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_digamma(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besj0(x):
    dll.x_besj0.restype = c_double
    dll.x_besj0.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_besj0(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besj1(x):
    dll.x_besj1.restype = c_double
    dll.x_besj1.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_besj1(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besjnu(nu, x):
    dll.x_besjnu.restype = c_double
    dll.x_besjnu.argtype = [c_double, c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_besjnu(c_double(nu), c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besy0(x):
    dll.x_besy0.restype = c_double
    dll.x_besy0.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_besy0(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besy1(x):
    dll.x_besy1.restype = c_double
    dll.x_besy1.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_besy1(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besynu(nu, x):
    dll.x_besynu.restype = c_double
    dll.x_besynu.argtype = [c_double, c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_besynu(c_double(nu), c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besi0(x, kode = 1):
    dll.x_besi0.restype = c_double
    dll.x_besi0.argtype = [c_double, POINTER(c_int)]
    dll.x_besi0e.restype = c_double
    dll.x_besi0e.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    if kode != 2:
        result = dll.x_besi0(c_double(x), byref(xerrno))
    else:
        result = dll.x_besi0e(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besi1(x, kode = 1):
    dll.x_besi1.restype = c_double
    dll.x_besi1.argtype = [c_double, POINTER(c_int)]
    dll.x_besi1e.restype = c_double
    dll.x_besi1e.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    if kode != 2:
        result = dll.x_besi1(c_double(x), byref(xerrno))
    else:
        result = dll.x_besi1e(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besinu(nu, x):
    dll.x_besinu.restype = c_double
    dll.x_besinu.argtype = [c_double, c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_besinu(c_double(nu), c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besk0(x, kode = 1):
    dll.x_besk0.restype = c_double
    dll.x_besk0.argtype = [c_double, POINTER(c_int)]
    dll.x_besk0e.restype = c_double
    dll.x_besk0e.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    if kode != 2:
        result = dll.x_besk0(c_double(x), byref(xerrno))
    else:
        result = dll.x_besk0e(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besk1(x, kode = 1):
    dll.x_besk1.restype = c_double
    dll.x_besk1.argtype = [c_double, POINTER(c_int)]
    dll.x_besk1e.restype = c_double
    dll.x_besk1e.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    if kode != 2:
        result = dll.x_besk1(c_double(x), byref(xerrno))
    else:
        result = dll.x_besk1e(c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def besknu(nu, x):
    dll.x_besknu.restype = c_double
    dll.x_besknu.argtype = [c_double, c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_besknu(c_double(nu), c_double(x), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def celli1(k):
    dll.x_celli1.restype = c_double
    dll.x_celli1.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_celli1(c_double(k), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def celli2(k):
    dll.x_celli2.restype = c_double
    dll.x_celli2.argtype = [c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_celli2(c_double(k), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def celli3(n, k):
    dll.x_celli3.restype = c_double
    dll.x_celli3.argtype = [c_double, c_double, POINTER(c_int)]
    xerrno = c_int()
    result = dll.x_celli3(c_double(n), c_double(k), byref(xerrno))
    if xerrno.value != 0:
        errno.value = xerrno.value
    return result
def dconst(i):
    dll._dconst.restype = c_double
    dll._dconst.argtype = [c_int]
    return dll._dconst(c_int(i))
def dlange(norm, m, n, a):
    dll._dlange.restype = c_double
    dll._dlange.argtype = [c_char, c_int, c_int, c_int, POINTER(c_double), POINTER(c_double)]
    anorm = 0.0
    info = 0
    if not norm in 'Mm1OoIiFfEe':
        info = -1
    elif m < 0:
        info = -2
    elif n < 0:
        info = -3
    elif a.ndim == 1:
        lda = m
        if a.shape[0] < lda*n:
            info = -4
    elif a.ndim == 2:
        lda = a.shape[1]
        if a.shape[0] < n:
            info = -4
    else:
        info = -4
    if info == 0:
        work = (c_double * m)()
        anorm = dll._dlange(c_char(norm.encode()), c_int(m), c_int(n), c_int(lda), a.ctypes.data_as(POINTER(c_double)), work)
    return (anorm, info)
def dlansy(norm, uplo, n, a):
    dll._dlansy.restype = c_double
    dll._dlansy.argtype = [c_char, c_char, c_int, c_int, POINTER(c_double), POINTER(c_double)]
    anorm = 0.0
    info = 0
    if not norm in 'Mm1OoIiFfEe':
        info = -1
    elif not uplo in 'UuLl':
        info = -2
    elif n < 0:
        info = -3
    elif a.ndim == 1:
        lda = n
        if a.shape[0] < lda*n:
            info = -4
    elif a.ndim == 2:
        lda = a.shape[1]
        if a.shape[0] < n:
            info = -4
    else:
        info = -4
    if info == 0:
        work = (c_double * n)()
        anorm = dll._dlansy(c_char(norm.encode()), c_char(uplo.encode()), c_int(n), c_int(lda), a.ctypes.data_as(POINTER(c_double)), work)
    return (anorm, info)
def dgesv(n, a, ipiv, b, nrhs = 1):
    dll._dgesv.argtype = [c_int, c_int, c_int, POINTER(c_double), POINTER(c_int), c_int, POINTER(c_double), POINTER(c_int)]
    if a.ndim == 1:
        lda = n
        if a.shape[0] < lda*n:
            return -2
    elif a.ndim == 2:
        lda = a.shape[1]
        if a.shape[0] < n:
            return -2
    else:
        return -2
    if b.ndim == 1:
        ldb = n
        if b.shape[0] < ldb*nrhs:
            return -4
    elif b.ndim == 2:
        ldb = b.shape[1]
        if b.shape[0] < nrhs:
            return -4
    else:
        return -4
    info = c_int()
    dll._dgesv(c_int(n), c_int(nrhs), c_int(lda), a.ctypes.data_as(POINTER(c_double)), ipiv.ctypes.data_as(POINTER(c_int)), c_int(ldb), b.ctypes.data_as(POINTER(c_double)), byref(info))
    if info.value == -2:
        info.value = -5
    elif info.value == -3:
        info.value = -2
    elif info.value == -6:
        info.value = -4
    return info.value
def dgecon(norm, n, a, anorm):
    dll._dgecon.argtype = [c_char, c_int, c_int, POINTER(c_double), c_double, POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_int)]
    if n < 0:
        return 0.0, -2
    if n == 0:
        return 1.0, 0
    if a.ndim == 1:
        lda = n
        if a.shape[0] < lda*n:
            return 0.0, -3
    elif a.ndim == 2:
        lda = a.shape[1]
        if a.shape[0] < n:
            return 0.0, -3
    else:
        return 0.0, -3
    work = (c_double * 4*n)()
    iwork = (c_int * n)()
    rcond = c_double(0.0)
    info = c_int()
    dll._dgecon(c_char(norm.encode()), c_int(n), c_int(lda), a.ctypes.data_as(POINTER(c_double)), c_double(anorm), byref(rcond), work, iwork, byref(info))
    if info.value == -5:
        info.value = -4
    return (rcond.value, info.value)
def dposv(uplo, n, a, b, nrhs = 1):
    dll._dposv.argtype = [c_char, c_int, c_int, c_int, POINTER(c_double), c_int, POINTER(c_double), POINTER(c_int)]
    if a.ndim == 1:
        lda = n
        if a.shape[0] < lda*n:
            return -3
    elif a.ndim == 2:
        lda = a.shape[1]
        if a.shape[0] < n:
            return -3
    else:
        return -3
    if b.ndim <= 1:
        ldb = min(n, b.shape[0])
    else:
        ldb = b.shape[1]
    info = c_int()
    dll._dposv(c_char(uplo.encode()), c_int(n), c_int(nrhs), c_int(lda), a.ctypes.data_as(POINTER(c_double)), c_int(ldb), b.ctypes.data_as(POINTER(c_double)), byref(info))
    if info.value == -6:
        info.value = -5
    return info.value
def dpocon(uplo, n, a, anorm):
    dll._dpocon.argtype = [c_char, c_int, c_int, POINTER(c_double), c_double, POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_int)]
    if n < 0:
        return 0.0, -2
    if n == 0:
        return 1.0, 0
    if a.ndim == 1:
        lda = n
        if a.shape[0] < lda*n:
            return 0.0, -3
    elif a.ndim == 2:
        lda = a.shape[1]
        if a.shape[0] < n:
            return 0.0, -3
    else:
        return 0.0, -3
    work = (c_double * 3*n)()
    iwork = (c_int * n)()
    rcond = c_double(0.0)
    info = c_int()
    dll._dpocon(c_char(uplo.encode()), c_int(n), c_int(lda), a.ctypes.data_as(POINTER(c_double)), c_double(anorm), byref(rcond), work, iwork, byref(info))
    if info.value == -5:
        info.value = -4
    return rcond.value, info.value
def dsyev(jobz, uplo, n, a, w):
    dll._dsyev.argtype = [c_char, c_char, c_int, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int, POINTER(c_int)]
    if a.ndim == 1:
        lda = n
        if a.shape[0] < lda*n:
            return -4
    elif a.ndim == 2:
        lda = a.shape[1]
        if a.shape[0] < n:
            return -4
    else:
        return -4
    lwork = -1
    work0 = c_double(0.0)
    info = c_int()
    dll._dsyev(c_char(jobz.encode()), c_char(uplo.encode()), c_int(n), c_int(lda), a.ctypes.data_as(POINTER(c_double)), w.ctypes.data_as(POINTER(c_double)), byref(work0), c_int(lwork), byref(info))
    if info.value == 0:
        lwork = int(work0.value)
        work = (c_double * lwork)()
        dll._dsyev(c_char(jobz.encode()), c_char(uplo.encode()), c_int(n), c_int(lda), a.ctypes.data_as(POINTER(c_double)), w.ctypes.data_as(POINTER(c_double)), work, c_int(lwork), byref(info))
    return info.value
def dgels(trans, m, n, a, b, nrhs = 1):
    dll._dgels.argtype = [c_char, c_int, c_int, c_int, c_int, POINTER(c_double), c_int, POINTER(c_double), POINTER(c_double), c_int, POINTER(c_int)]
    if a.ndim == 1:
        lda = m
        if a.shape[0] < lda*n:
            return -4
    elif a.ndim == 2:
        lda = a.shape[1]
        if a.shape[0] < n:
            return -4
    else:
        return -4
    if b.ndim == 1:
        ldb = max(m, n)
        if b.shape[0] < ldb*nrhs:
            return -5
    elif b.ndim == 2:
        ldb = b.shape[1]
        if b.shape[0] < nrhs:
            return -5
    else:
        return -5
    lwork = -1
    work0 = c_double(0.0)
    info = c_int()
    dll._dgels(c_char(trans.encode()), c_int(m), c_int(n), c_int(nrhs), c_int(lda), a.ctypes.data_as(POINTER(c_double)), c_int(ldb), b.ctypes.data_as(POINTER(c_double)), byref(work0), c_int(lwork), byref(info))
    if info.value == 0:
        lwork = int(work0.value)
        work = (c_double * lwork)()
        dll._dgels(c_char(trans.encode()), c_int(m), c_int(n), c_int(nrhs), c_int(lda), a.ctypes.data_as(POINTER(c_double)), c_int(ldb), b.ctypes.data_as(POINTER(c_double)), work, c_int(lwork), byref(info))
    if info.value == -7:
        info.value = -6
    return info.value
def dgecov(job, n, a, ci):
    dll._dgecov.argtype = [c_int, c_int, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_int)]
    if a.ndim == 1:
        lda = n
        if a.shape[0] < lda*n:
            return -3
    elif a.ndim == 2:
        lda = a.shape[1]
        if a.shape[0] < n:
            return -3
    else:
        return -3
    info = c_int()
    dll._dgecov(c_int(job), c_int(n), c_int(lda), a.ctypes.data_as(POINTER(c_double)), ci.ctypes.data_as(POINTER(c_double)), byref(info))
    return info.value
def pchse(n, x, f, d, incfd = 1):
    dll._pchse.argtype = [c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int, POINTER(c_int)]
    lwork = 2*n
    work = (c_double * lwork)()
    info = c_int()
    dll._pchse(c_int(n), x.ctypes.data_as(POINTER(c_double)), f.ctypes.data_as(POINTER(c_double)), d.ctypes.data_as(POINTER(c_double)), c_int(incfd), work, c_int(lwork), byref(info))
    return info.value
def pchfe(n, x, f, d, ne, xe, fe, incfd = 1, skip = 0):
    dll._pchfe.argtype = [c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int, c_int, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_int)]
    info = c_int()
    dll._pchfe(c_int(n), x.ctypes.data_as(POINTER(c_double)), f.ctypes.data_as(POINTER(c_double)), d.ctypes.data_as(POINTER(c_double)), c_int(incfd), c_int(skip), c_int(ne), xe.ctypes.data_as(POINTER(c_double)), fe.ctypes.data_as(POINTER(c_double)), byref(info))
    return info.value
def pchia(n, x, f, d, a, b, incfd = 1, skip = 0):
    dll._pchia.restype = c_double
    dll._pchia.argtype = [c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int, c_int, c_double, c_double, POINTER(c_int)]
    info = c_int()
    s = dll._pchia(c_int(n), x.ctypes.data_as(POINTER(c_double)), f.ctypes.data_as(POINTER(c_double)), d.ctypes.data_as(POINTER(c_double)), c_int(incfd), c_int(skip), c_double(a), c_double(b), byref(info))
    return (s, info.value)
def rpzero2(n, a, rr, ri, s, iflag = 0, maxiter = 100):
    dll._rpzero2.argtype = [c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int, c_int, POINTER(c_int), POINTER(c_double), POINTER(c_double), POINTER(c_int)]
    iter = c_int(0)
    lwork = 8*n + 6
    work = (c_double * lwork)()
    info = c_int()
    dll._rpzero2(c_int(n), a.ctypes.data_as(POINTER(c_double)), rr.ctypes.data_as(POINTER(c_double)), ri.ctypes.data_as(POINTER(c_double)), c_int(iflag), c_int(maxiter), byref(iter), s.ctypes.data_as(POINTER(c_double)), work, byref(info))
    return (iter.value, info.value)
def dfzero(f, b, c, r, re = 1.0e-10, ae = 1.0e-10):
    dll._dfzero_r.argtype = [POINTER(c_double), POINTER(c_double), c_double, c_double, c_double, POINTER(c_int), POINTER(c_double), c_double, POINTER(c_int)]
    b1 = c_double(b)
    c1 = c_double(c)
    info = c_int()
    xx = c_double(0.0)
    yy = c_double(0.0)
    irev = c_int(0)
    while True:
        dll._dfzero_r(byref(b1), byref(c1), c_double(r), c_double(re), c_double(ae), byref(info), byref(xx), yy, byref(irev))
        if irev.value == 0:
            break
        yy = c_double(f(xx.value))
    return (b1.value, info.value)
def hybrd1(f, n, x, fvec, xtol = 1.0e-10):
    dll._hybrd1_r.argtype = [c_int, POINTER(c_double), POINTER(c_double), c_double, POINTER(c_double), c_int, POINTER(c_int), POINTER(c_double), POINTER(c_double), POINTER(c_int)]
    if n < 1:
        return -2
    if xtol < 0:
        return -5
    lwork = int(n*(3*n + 13)/2)
    work = (c_double * lwork)()
    info = c_int(0)
    xx = np.empty(n)
    yy = np.empty(n)
    irev = c_int(0)
    while True:
        dll._hybrd1_r(c_int(n), x.ctypes.data_as(POINTER(c_double)), fvec.ctypes.data_as(POINTER(c_double)), c_double(xtol), work, c_int(lwork), byref(info), xx.ctypes.data_as(POINTER(c_double)), yy.ctypes.data_as(POINTER(c_double)), byref(irev))
        if irev.value == 1 or irev.value == 2:
            iflag = 1
        elif irev.value == 3 or irev.value == 4:
            iflag = 2
        elif irev.value == 0:
            break
        f(n, xx, yy, iflag)
    return info.value
def dfmin(a, b, f, tol = 1.0e-10):
    dll._dfmin_r.argtype = [c_double, c_double, c_double, POINTER(c_double), c_double, POINTER(c_int)]
    xx = c_double(0)
    yy = c_double(0)
    irev = c_int(0)
    while True:
        dll._dfmin_r(c_double(a), c_double(b), c_double(tol), byref(xx), yy, byref(irev))
        if irev.value == 0:
            break
        yy = c_double(f(xx.value))
    return (xx.value, 0)
def optif0(n, x, f, xpls):
    dll._optif0_r.argtype = [c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int, POINTER(c_int), POINTER(c_double), c_double, POINTER(c_int)]
    if n < 1:
        return 0.0, -1
    fpls = c_double(0)
    lwork = n*(n + 10)
    work = (c_double * lwork)()
    info = c_int(0)
    xx = np.empty(n)
    yy = c_double(0)
    irev = c_int(0)
    while True:
        dll._optif0_r(c_int(n), x.ctypes.data_as(POINTER(c_double)), xpls.ctypes.data_as(POINTER(c_double)), byref(fpls), work, c_int(lwork), byref(info), xx.ctypes.data_as(POINTER(c_double)), yy, byref(irev))
        if irev.value >= 1 and irev.value <= 20:
            yy = c_double(f(n, xx))
        elif irev.value == 0:
            break
    return fpls.value, info.value
def qk15(f, a, b):
    dll._qk15_r.argtype = [c_double, c_double, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), c_double, POINTER(c_int)]
    result = c_double(0)
    abserr = c_double(0)
    resabs = c_double(0)
    resasc = c_double(0)
    xx = c_double(0)
    yy = c_double(0)
    irev = c_int(0)
    while True:
        dll._qk15_r(c_double(a), c_double(b), byref(result), byref(abserr), byref(resabs), byref(resasc), byref(xx), yy, byref(irev))
        if irev.value >= 1 and irev.value <= 5:
            yy = c_double(f(xx.value))
        elif irev.value == 0:
            break
    return result.value, abserr.value
def qag(f, a, b, epsabs = 1.0e-10, epsrel = 1.0e-10, key = 1, limit = 100):
    dll._qag_r.argtype = [c_double, c_double, c_double, c_double, c_int, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_double), c_int, POINTER(c_int), POINTER(c_int), POINTER(c_double), c_double, POINTER(c_int)]
    result = c_double(0)
    abserr = c_double(0)
    neval = c_int(0)
    last = c_int(0)
    lwork = 4*limit
    work = (c_double * lwork)()
    liwork = limit
    iwork = (c_int * liwork)()
    info = c_int(0)
    xx = c_double(0)
    yy = c_double(0)
    irev = c_int(0)
    while True:
        dll._qag_r(c_double(a), c_double(b), c_double(epsabs), c_double(epsrel), c_int(key), c_int(limit), byref(result), byref(abserr), byref(neval), byref(last), work, c_int(lwork), iwork, byref(info), byref(xx), yy, byref(irev))
        if irev.value >= 1 and irev.value <= 15:
            yy = c_double(f(xx.value))
        elif irev.value == 0:
            break
    if info.value == -6:
        info.value = -7
    return result.value, abserr.value, info.value
def qagi(f, bound, inf, epsabs = 1.0e-10, epsrel = 1.0e-10, limit = 100):
    dll._qagi_r.argtype = [c_double, c_int, c_double, c_double, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_double), c_int, POINTER(c_int), POINTER(c_int), POINTER(c_double), c_double, POINTER(c_int)]
    result = c_double(0)
    abserr = c_double(0)
    neval = c_int(0)
    last = c_int(0)
    lwork = 4*limit
    work = (c_double * lwork)()
    liwork = limit
    iwork = (c_int * liwork)()
    info = c_int(0)
    xx = c_double(0)
    yy = c_double(0)
    irev = c_int(0)
    while True:
        dll._qagi_r(c_double(bound), c_int(inf), c_double(epsabs), c_double(epsrel), c_int(limit), byref(result), byref(abserr), byref(neval), byref(last), work, c_int(lwork), iwork, byref(info), byref(xx), yy, byref(irev))
        if irev.value >= 1 and irev.value <= 18:
            yy = c_double(f(xx.value))
        elif irev.value == 0:
            break
    if info.value == -5:
        info.value = -6
    return result.value, abserr.value, info.value
def derkfa(info, n, f, t, y, tout, tend, wsave, iwsave, mode = 2, rtol = 1.0e-10, atol = 1.0e-10):
    dll._derkfa_r.argtype = [c_int, POINTER(c_double), POINTER(c_double), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, POINTER(c_double), c_int, POINTER(c_int), c_int, POINTER(c_int), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_int)]
    if info < 0 or info > 10:
        return t, -1
    if n < 1:
        return t, -2
    if y.ndim != 1 or y.shape[0] < n:
        return t, -5
    if wsave.ndim != 1:
        return t, -8
    lwsave = wsave.shape[0]
    if iwsave.ndim != 1:
        return t, -9
    liwsave = iwsave.shape[0]
    t1 = c_double(t)
    rtol1 = c_double(rtol)
    atol1 = c_double(atol)
    itol = c_int(0)
    info1 = c_int(info)
    tt = c_double(0)
    yy = np.empty(n)
    yyp = np.empty(n)
    irev = c_int(0)
    while True:
        dll._derkfa_r(c_int(n), byref(t1), y.ctypes.data_as(POINTER(c_double)), c_double(tout), c_double(tend), byref(rtol1), byref(atol1), itol, c_int(mode), wsave.ctypes.data_as(POINTER(c_double)), c_int(-lwsave), iwsave.ctypes.data_as(POINTER(c_int)), c_int(-liwsave), byref(info1), byref(tt), yy.ctypes.data_as(POINTER(c_double)), yyp.ctypes.data_as(POINTER(c_double)), byref(irev))
        if irev.value >= 1 and irev.value <= 10:
            f(n, tt.value, yy, yyp)
        elif irev.value == 0:
            break
    if info1.value == -14:
        info1.value = -1
    elif info1.value == -9:
        info1.value = -10
    elif info1.value == -11:
        info1.value = -8
    elif info1.value == -13:
        info1.value = -9
    elif info1.value == -6:
        info1.value = -11
    elif info1.value == -7:
        info1.value = -12
    return t1.value, info1.value
def dopn43(info, n, f2, t, y, yp, tout, tend, wsave, iwsave, mode = 2, rtol = 1.0e-10, atol = 1.0e-10):
    dll._dopn43_r.argtype = [c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), c_double, c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, POINTER(c_double), c_int, POINTER(c_int), c_int, POINTER(c_int), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_int)]
    if info < 0 or info > 10:
        return t, -1
    if n < 1:
        return t, -2
    if y.ndim != 1 or y.shape[0] < n:
        return t, -5
    if wsave.ndim != 1:
        return t, -8
    lwsave = wsave.shape[0]
    if iwsave.ndim != 1:
        return t, -9
    liwsave = iwsave.shape[0]
    t1 = c_double(t)
    rtol1 = c_double(rtol)
    atol1 = c_double(atol)
    itol = c_int(0)
    info1 = c_int(info)
    tt = c_double(0)
    yy = np.empty(n)
    yypp = np.empty(n)
    irev = c_int(0)
    while True:
        dll._dopn43_r(c_int(n), byref(t1), y.ctypes.data_as(POINTER(c_double)), yp.ctypes.data_as(POINTER(c_double)), c_double(tout), c_double(tend), byref(rtol1), byref(atol1), itol, c_int(mode), wsave.ctypes.data_as(POINTER(c_double)), c_int(-lwsave), iwsave.ctypes.data_as(POINTER(c_int)), c_int(-liwsave), byref(info1), byref(tt), yy.ctypes.data_as(POINTER(c_double)), yypp.ctypes.data_as(POINTER(c_double)), byref(irev))
        if irev.value >= 1 and irev.value <= 2:
            f2(n, tt.value, yy, yypp)
        elif irev.value == 0:
            break
    if info1.value == -15:
        info1.value = -1
    elif info1.value == -10:
        info1.value = -11
    elif info1.value == -12:
        info1.value = -9
    elif info1.value == -14:
        info1.value = -10
    elif info1.value == -7:
        info1.value = -12
    elif info1.value == -8:
        info1.value = -13
    return t1.value, info1.value
def rfft1f(n, r, wsave, inc = 1):
    dll._rfft1f.argtype = [c_int, c_int, POINTER(c_double), c_int, POINTER(c_double), c_int, POINTER(c_double), c_int, POINTER(c_int)]
    if r.ndim != 1:
        return -2
    if wsave.ndim != 1:
        return -3
    lr = r.shape[0]
    lwsave = wsave.shape[0]
    lwork = n
    work = (c_double * lwork)()
    info = c_int()
    dll._rfft1f(c_int(n), c_int(inc), r.ctypes.data_as(POINTER(c_double)), c_int(lr), wsave.ctypes.data_as(POINTER(c_double)), c_int(lwsave), work, c_int(lwork), byref(info))
    if info.value == -2:
        info.value = -4
    elif info.value == -4:
        info.value = -2
    elif info.value == -6:
        info.value = -3
    return info.value
def rfft1b(n, r, wsave, inc = 1):
    dll._rfft1b.argtype = [c_int, c_int, POINTER(c_double), c_int, POINTER(c_double), c_int, POINTER(c_double), c_int, POINTER(c_int)]
    if r.ndim != 1:
        return -2
    if wsave.ndim != 1:
        return -3
    lr = r.shape[0]
    lwsave = wsave.shape[0]
    lwork = n
    work = (c_double * lwork)()
    info = c_int()
    dll._rfft1b(c_int(n), c_int(inc), r.ctypes.data_as(POINTER(c_double)), c_int(lr), wsave.ctypes.data_as(POINTER(c_double)), c_int(lwsave), work, c_int(lwork), byref(info))
    if info.value == -2:
        info.value = -4
    elif info.value == -4:
        info.value = -2
    elif info.value == -6:
        info.value = -3
    return info.value
def rfft1i(n, wsave):
    dll._rfft1i.argtype = [c_int, POINTER(c_double), c_int, POINTER(c_int)]
    if wsave.ndim != 1:
        return -3
    lwsave = wsave.shape[0]
    info = c_int()
    dll._rfft1i(c_int(n), wsave.ctypes.data_as(POINTER(c_double)), c_int(lwsave), byref(info))
    if info.value == -3:
        info.value = -2
    return info.value
def lmdif1(f, m, n, x, fvec, tol = 1.0e-10):
    dll._lmdif1_r.argtype = [c_int, c_int, POINTER(c_double), POINTER(c_double), c_double, POINTER(c_double), c_int, POINTER(c_int), POINTER(c_int), POINTER(c_double), POINTER(c_double), POINTER(c_int)]
    if m < n:
        return -2
    if n < 1:
        return -3
    if tol < 0:
        return -6
    lwork = n*(m + 5) + m
    work = (c_double * lwork)()
    liwork = n
    iwork = (c_int * liwork)()
    info = c_int(0)
    xx = np.empty(n)
    yy = np.empty(m)
    irev = c_int(0)
    while True:
        dll._lmdif1_r(c_int(m), c_int(n), x.ctypes.data_as(POINTER(c_double)), fvec.ctypes.data_as(POINTER(c_double)), c_double(tol), work, c_int(lwork), iwork, byref(info), xx.ctypes.data_as(POINTER(c_double)), yy.ctypes.data_as(POINTER(c_double)), byref(irev))
        if irev.value == 1 or irev.value == 2:
            iflag = 1
        elif irev.value == 3:
            iflag = 2
        elif irev.value == 0:
            break
        f(m, n, xx, yy, iflag)
    return info.value
def init_genrand(s):
    dll._init_genrand.argtype = [c_uint]
    dll._init_genrand(c_ulong(s))
def genrand_int32():
    dll._genrand_int32.restype = c_int
    dll._genrand_int32.argtype = []
    return dll._genrand_int32()
def genrand_int31():
    dll._genrand_int31.restype = c_int
    dll._genrand_int31.argtype = []
    return dll._genrand_int31()
def genrand_res53():
    dll._genrand_res53.restype = c_double
    dll._genrand_res53.argtype = []
    return dll._genrand_res53()
def dlamch(c):
    dll._dlamch.restype = c_double
    dll._dlamch.argtype = [c_char]
    return dll._dlamch(c_wchar(c))
def get_errno():
    return errno.value
def set_errno(n):
    errno.value = n
def csr_dusmv(trans, m, n, alpha, val, rowptr, colind, base, x, incx, beta, y, incy):
    dllsp._csr_dusmv.argtype = [c_char, c_int, c_int, c_double, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_double), c_int, c_double, POINTER(c_double), c_int, POINTER(c_int)]
    info = c_int(0)
    dllsp._csr_dusmv(c_char(trans.encode()), c_int(m), c_int(n), c_double(alpha), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), x.ctypes.data_as(POINTER(c_double)), c_int(incx), c_double(beta), y.ctypes.data_as(POINTER(c_double)), c_int(incy), byref(info))
    return info.value
def csr_dussv(uplo, trans, diag, n, val, rowptr, colind, base, x, incx):
    dllsp._csr_dussv.argtype = [c_char, c_char, c_char, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_double), c_int, POINTER(c_int)]
    info = c_int(0)
    dllsp._csr_dussv(c_char(uplo.encode()), c_char(trans.encode()), c_char(diag.encode()), c_int(n), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), x.ctypes.data_as(POINTER(c_double)), c_int(incx), byref(info))
    return info.value
def csr_dusmm(trans, order, m, n, l, alpha, val, rowptr, colind, base, b, beta, c):
    dllsp._csr_dusmm.argtype = [c_char, c_char, c_int, c_int, c_int, c_double, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, c_int, POINTER(c_double), c_double, c_int, POINTER(c_double), POINTER(c_int)]
    if l < 0:
        return -5
    if l == 0:
        return 0
    if b.ndim == 2:
        ldb = b.shape[1]
        if b.shape[0] < l:
            return -11
    else:
        return -11
    if c.ndim == 2:
        ldc = c.shape[1]
        if c.shape[0] < l:
            return -12
    else:
        return -12
    info = c_int(0)
    dllsp._csr_dusmm(c_char(trans.encode()), c_char(order.encode()), c_int(m), c_int(n), c_int(l), c_double(alpha), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), c_int(ldb), b.ctypes.data_as(POINTER(c_double)), c_double(beta), c_int(ldc), c.ctypes.data_as(POINTER(c_double)), byref(info))
    if info.value == -14:
        info.value = -13
    return info.value
def csr_dussm(uplo, trans, diag, order, n, nrhs, val, rowptr, colind, base, x):
    dllsp._csr_dussm.argtype = [c_char, c_char, c_char, c_char, c_int, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, c_int, POINTER(c_double), POINTER(c_int)]
    if nrhs < 0:
        return -6
    if nrhs == 0:
        return 0
    if x.ndim == 2:
        ldx = x.shape[1]
        if x.shape[0] < nrhs:
            return -11
    else:
        return -11
    info = c_int(0)
    dllsp._csr_dussm(c_char(uplo.encode()), c_char(trans.encode()), c_char(diag.encode()), c_char(order.encode()), c_int(n), c_int(nrhs), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), c_int(ldx), x.ctypes.data_as(POINTER(c_double)), byref(info))
    return info.value
def ssr_dusmv(uplo, n, alpha, val, rowptr, colind, base, x, incx, beta, y, incy):
    dllsp._ssr_dusmv.argtype = [c_char, c_int, c_double, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_double), c_int, c_double, POINTER(c_double), c_int, POINTER(c_int)]
    info = c_int(0)
    dllsp._ssr_dusmv(c_char(uplo.encode()), c_int(n), c_double(alpha), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), x.ctypes.data_as(POINTER(c_double)), c_int(incx), c_double(beta), y.ctypes.data_as(POINTER(c_double)), c_int(incy), byref(info))
    return info.value
def csr_dense(m, n, val, rowptr, colind, base, a):
    dllsp._csr_dense.argtype = [c_int, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, c_int, POINTER(c_double), POINTER(c_int)]
    if n < 0:
        return -2
    if n == 0:
        return 0
    if a.ndim == 2:
        lda = a.shape[1]
        if a.shape[0] < n:
            return -7
    else:
        return -7
    info = c_int(0)
    dllsp._csr_dense(c_int(m), c_int(n), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), c_int(lda), a.ctypes.data_as(POINTER(c_double)), byref(info))
    return info.value
def dense_csr(m, n, a, maxnnz, val, rowptr, colind, base):
    dllsp._dense_csr.argtype = [c_int, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, c_int, POINTER(c_double), POINTER(c_int)]
    if n < 0:
        return -2
    if n == 0:
        return 0
    if a.ndim == 2:
        lda = a.shape[1]
        if a.shape[0] < n:
            return -3
    else:
        return -3
    info = c_int(0)
    dllsp._dense_csr(c_int(m), c_int(n), c_int(lda), a.ctypes.data_as(POINTER(c_double)), c_int(maxnnz), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), byref(info))
    return info.value
def csr_coo(m, n, val, rowptr, colind, base, val2, rowind2, colind2, base2):
    dllsp._csr_coo.argtype = [c_int, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_int)]
    info = c_int(0)
    dllsp._csr_coo(c_int(m), c_int(n), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), val2.ctypes.data_as(POINTER(c_double)), rowind2.ctypes.data_as(POINTER(c_int)), colind2.ctypes.data_as(POINTER(c_int)), c_int(base2), byref(info))
    return info.value
def coo_csr(m, n, nnz, val, rowind, colind, base, val2, rowptr2, colind2, base2):
    dllsp._coo_csr.argtype = [c_int, c_int, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_int)]
    info = c_int(0)
    dllsp._coo_csr(c_int(m), c_int(n), c_int(nnz), val.ctypes.data_as(POINTER(c_double)), rowind.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), val2.ctypes.data_as(POINTER(c_double)), rowptr2.ctypes.data_as(POINTER(c_int)), colind2.ctypes.data_as(POINTER(c_int)), c_int(base2), byref(info))
    return info.value
def csr_trans(m, n, val, ptr, ind, base, val2, ptr2, ind2, base2):
    dllsp._csc_csr.argtype = [c_int, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_int)]
    info = c_int(0)
    dllsp._csc_csr(c_int(n), c_int(m), val.ctypes.data_as(POINTER(c_double)), ptr.ctypes.data_as(POINTER(c_int)), ind.ctypes.data_as(POINTER(c_int)), c_int(base), val2.ctypes.data_as(POINTER(c_double)), ptr2.ctypes.data_as(POINTER(c_int)), ind2.ctypes.data_as(POINTER(c_int)), c_int(base2), byref(info))
    return info.value
def csr_ssr(uplo, n, val, rowptr, colind, base, val2, rowptr2, colind2, base2):
    dllsp._csr_ssr.argtype = [c_char, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_int)]
    if val2.ndim == 1:
        nnz2 = val2.shape[0]
    else:
        return -7
    if colind2.ndim == 1:
        if nnz2 > colind2.shape[0]:
            nnz2 = colind2.shape[0]
    else:
        return -9
    info = c_int(0)
    dllsp._csr_ssr(c_char(uplo.encode()), c_int(n), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), c_int(nnz2), val2.ctypes.data_as(POINTER(c_double)), rowptr2.ctypes.data_as(POINTER(c_int)), colind2.ctypes.data_as(POINTER(c_int)), c_int(base2), byref(info))
    if info.value <= -8:
        info.value = info.value + 1
    return info.value
def ssr_csr(uplo, n, val, rowptr, colind, base, val2, rowptr2, colind2, base2):
    dllsp._ssr_csr.argtype = [c_char, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_int)]
    if val2.ndim == 1:
        nnz2 = val2.shape[0]
    else:
        return -7
    if colind2.ndim == 1:
        if nnz2 > colind2.shape[0]:
            nnz2 = colind2.shape[0]
    else:
        return -9
    info = c_int(0)
    dllsp._ssr_csr(c_char(uplo.encode()), c_int(n), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), c_int(nnz2), val2.ctypes.data_as(POINTER(c_double)), rowptr2.ctypes.data_as(POINTER(c_int)), colind2.ctypes.data_as(POINTER(c_int)), c_int(base2), byref(info))
    if info.value <= -8:
        info.value = info.value + 1
    return info.value
def bicg1(n, val, rowptr, colind, b, x, tol = 1.0e-10, maxiter = 500):
    dllsp._bicg1.argtype = [c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_double), POINTER(c_double), c_double, c_int, POINTER(c_int), POINTER(c_double), c_int, POINTER(c_double), POINTER(c_int)]
    if n < 1:
        return -1, 0, 0.0
    lwork = 8*n
    work = (c_double * lwork)()
    iter = c_int(0)
    res = c_double(0)
    info = c_int(0)
    dllsp._bicg1(c_int(n), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), b.ctypes.data_as(POINTER(c_double)), x.ctypes.data_as(POINTER(c_double)), c_double(tol), c_int(maxiter), byref(iter), byref(res), c_int(lwork), work, byref(info))
    if info.value == -11:
        info.value = -9
    return info.value, iter.value, res.value
def cg1(uplo, n, val, rowptr, colind, b, x, tol = 1.0e-10, maxiter = 500):
    dllsp._cg1.argtype = [c_char, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_double), POINTER(c_double), c_double, c_int, POINTER(c_int), POINTER(c_double), c_int, POINTER(c_double), POINTER(c_int)]
    if n < 1:
        return -2, 0, 0.0
    lwork = 5*n
    work = (c_double * lwork)()
    iter = c_int(0)
    res = c_double(0)
    info = c_int(0)
    dllsp._cg1(c_char(uplo.encode()), c_int(n), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), b.ctypes.data_as(POINTER(c_double)), x.ctypes.data_as(POINTER(c_double)), c_double(tol), c_int(maxiter), byref(iter), byref(res), c_int(lwork), work, byref(info))
    if info.value == -12:
        info.value = -10
    return info.value, iter.value, res.value
def fem2p(n, ne, x, y, knc, p, q, f, nb1, ib, bv, nb2, ks2, alpha, beta, val, rowptr, colind, rhs, base = 0):
    dllpde._fem2p.argtype = [c_int, c_int, POINTER(c_double), POINTER(c_double), c_int, POINTER(c_int), POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int, POINTER(c_int), POINTER(c_double), c_int, c_int, POINTER(c_int), c_int, POINTER(c_double), c_int, POINTER(c_double), c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int)]
    if knc.ndim == 2:
        ldknc = knc.shape[1]
    else:
        return -5
    if ks2.ndim == 2:
        ldks2 = ks2.shape[1]
    else:
        return -13
    if alpha.ndim == 2:
        ldalpha = alpha.shape[1]
    else:
        return -14
    if beta.ndim == 2:
        ldbeta = beta.shape[1]
    else:
        return -15
    if val.ndim == 1:
        mxnnz = val.shape[0]
    else:
        return -16
    if colind.ndim == 1:
        if mxnnz > colind.shape[0]:
            mxnnz = colind.shape[0]
    else:
        return -18
    if rowptr.ndim != 1:
        return -17
    iwork = (c_int * 2 * n)()
    info = c_int(0)
    dllpde._fem2p(c_int(n), c_int(ne), x.ctypes.data_as(POINTER(c_double)), y.ctypes.data_as(POINTER(c_double)), c_int(ldknc), knc.ctypes.data_as(POINTER(c_int)), p.ctypes.data_as(POINTER(c_double)), q.ctypes.data_as(POINTER(c_double)), f.ctypes.data_as(POINTER(c_double)), c_int(nb1), ib.ctypes.data_as(POINTER(c_int)), bv.ctypes.data_as(POINTER(c_double)), c_int(nb2), c_int(ldks2), ks2.ctypes.data_as(POINTER(c_int)), c_int(ldalpha), alpha.ctypes.data_as(POINTER(c_double)), c_int(ldbeta), beta.ctypes.data_as(POINTER(c_double)), c_int(mxnnz), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), c_int(base), rhs.ctypes.data_as(POINTER(c_double)), iwork, byref(info))
    if info.value == -10:
        info.value = -9
    elif info.value == -13:
        info.value = -12
    elif info.value == -14:
        info.value = -13
    elif info.value == -16:
        info.value = -14
    elif info.value == -18:
        info.value = -15
    elif info.value == -20:
        info.value = -16
    elif info.value == -24:
        info.value = -20
    return info.value
def mesh23(nx, ny, x, y, knc, ks, lb, sclx = 1.0, scly = 1.0):
    dllpde._mesh23.argtype = [c_int, c_int, POINTER(c_int), c_double, c_double, POINTER(c_double), POINTER(c_double), POINTER(c_int), c_int, POINTER(c_int), c_int, POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)]
    if knc.ndim == 2:
        ldknc = knc.shape[1]
    else:
        return -5
    if ks.ndim == 2:
        ldks = ks.shape[1]
    else:
        return -6
    n = c_int(0)
    ne = c_int(0)
    nb = c_int(0)
    info = c_int(0)
    dllpde._mesh23(c_int(nx), c_int(ny), byref(n), c_double(sclx), c_double(scly), x.ctypes.data_as(POINTER(c_double)), y.ctypes.data_as(POINTER(c_double)), byref(ne), c_int(ldknc), knc.ctypes.data_as(POINTER(c_int)), c_int(ldks), ks.ctypes.data_as(POINTER(c_int)), lb.ctypes.data_as(POINTER(c_int)), byref(nb), byref(info))
    if info.value == -9:
        info.value = -5
    elif info.value == -11:
        info.value = -6
    return (info.value, n.value, ne.value, nb.value)
def mm_read1(fname, val, rowptr, colind, base = 0, sort = 0):
    dllsp._mm_read.argtype = [POINTER(c_char), POINTER(c_char), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_double), c_int, POINTER(c_int), c_int, POINTER(c_int), c_int, c_int, c_int, c_int, c_int, POINTER(c_int)]
    matcode = np.empty(4, dtype=bytes)
    nrow = c_int(0)
    ncol = c_int(0)
    nnz = c_int(0)
    if val.ndim == 1:
        lval = val.shape[0]
    else:
        return -2
    if rowptr.ndim == 1:
        lrowptr = rowptr.shape[0]
    else:
        return -3
    if colind.ndim == 1:
        lcolind = colind.shape[0]
    else:
        return -4
    skip = c_int(1)
    format = c_int(0)
    info = c_int(0)
    dllsp._mm_read(c_char_p(fname.encode()), matcode.ctypes.data_as(POINTER(c_char)), byref(nrow), byref(ncol), byref(nnz), val.ctypes.data_as(POINTER(c_double)), c_int(lval), rowptr.ctypes.data_as(POINTER(c_int)), c_int(lrowptr), colind.ctypes.data_as(POINTER(c_int)), c_int(lcolind), skip, c_int(base), format, c_int(sort), byref(info))
    if info.value == 0:
        if matcode[1] == b'C' and matcode[2] == b'R' and (matcode[3] == b'G' or matcode[3] == b'S'):
            skip = c_int(0)
            dllsp._mm_read(c_char_p(fname.encode()), matcode.ctypes.data_as(POINTER(c_char)), byref(nrow), byref(ncol), byref(nnz), val.ctypes.data_as(POINTER(c_double)), c_int(lval), rowptr.ctypes.data_as(POINTER(c_int)), c_int(lrowptr), colind.ctypes.data_as(POINTER(c_int)), c_int(lcolind), skip, c_int(base), format, c_int(sort), byref(info))
        else:
            info.value = 5
    if info.value == -7:
        info.value = -2
    elif info.value == -9:
        info.value = -3
    elif info.value == -11:
        info.value = -4
    elif info.value == -13:
        info.value = -5
    elif info.value == -15:
        info.value = -6
    return (info.value, nrow.value, ncol.value)
def mm_read_info(fname):
    dllsp._mm_read.argtype = [POINTER(c_char), POINTER(c_char), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_double), c_int, POINTER(c_int), c_int, POINTER(c_int), c_int, c_int, c_int, c_int, c_int, POINTER(c_int)]
    matcode = np.empty(4, dtype=bytes)
    nrow = c_int(0)
    ncol = c_int(0)
    nnz = c_int(0)
    val = c_double(0.0)
    lval = 0
    rowptr = c_int(0)
    lrowptr = 0
    colind = c_int(0)
    lcolind = 0
    skip = c_int(1)
    base = c_int(0)
    format = c_int(0)
    sort = c_int(0)
    info = c_int(0)
    dllsp._mm_read(c_char_p(fname.encode()), matcode.ctypes.data_as(POINTER(c_char)), byref(nrow), byref(ncol), byref(nnz), byref(val), lval, byref(rowptr), lrowptr, byref(colind), lcolind, skip, base, format, sort, byref(info))
    mform = 0
    dtype = 0
    mshape = 0
    if info.value == 0:
        if matcode[1] == b'A':
            mform = 1
        elif matcode[1] != b'C':
            info = 4
        if matcode[2] == b'C':
            dtype = 1
        elif matcode[2] == b'P':
            dtype = 2
        elif matcode[2] == b'I':
            dtype = 3
        elif matcode[2] != b'R':
            info = 4
        if matcode[3] == b'S':
            mshape = 1
        elif matcode[3] == b'Z':
            mshape = 2
        elif matcode[3] == b'H':
            mshape = 3
        elif matcode[3] != b'G':
            info = 4
    return (info.value, nrow.value, ncol.value, nnz.value, dtype, mshape, mform)
def mm_write1(fname, nrow, ncol, val, rowptr, colind, fchk = 0, mshape = 0):
    dllsp._mm_write.argtype = [POINTER(c_char), POINTER(c_char), c_int, c_int, c_int, POINTER(c_double), c_int, POINTER(c_int), c_int, POINTER(c_int), c_int, c_int, c_int, c_int, POINTER(c_int)]
    matcode = np.array(['M', 'C', 'R', 'G'], dtype=bytes)
    nnz = rowptr[nrow] - rowptr[0]
    if val.ndim == 1:
        lval = val.shape[0]
    else:
        return -4
    if rowptr.ndim == 1:
        lrowptr = rowptr.shape[0]
    else:
        return -5
    if colind.ndim == 1:
        lcolind = colind.shape[0]
    else:
        return -6
    base = c_int(0)
    format = c_int(0)
    if mshape == 1:
        matcode[3] = 'S'
    elif mshape == 2:
        matcode[3] = 'Z'
    elif mshape != 0:
        return -8
    info = c_int(0)
    dllsp._mm_write(c_char_p(fname.encode()), matcode.ctypes.data_as(POINTER(c_char)), c_int(nrow), c_int(ncol), c_int(nnz), val.ctypes.data_as(POINTER(c_double)), c_int(lval), rowptr.ctypes.data_as(POINTER(c_int)), c_int(lrowptr), colind.ctypes.data_as(POINTER(c_int)), c_int(lcolind), base, c_int(fchk), format, byref(info))
    if info.value == -3:
        info.value = -2
    elif info.value == -4:
        info.value = -3
    elif info.value == -5 or info.value == -12:
        info.value = -5
    elif info.value == -7:
        info.value = -4
    elif info.value == -9:
        info.value = -5
    elif info.value == -11:
        info.value = -6
    elif info.value == -13:
        info.value = -7
    return info.value
def writevtkug(fname, n, x, y, z, ne, kc, u):
    dllpde._writevtkug.argtype = [POINTER(c_char), c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int, c_int, POINTER(c_int), POINTER(c_double), POINTER(c_int)]
    if kc.ndim == 2:
        ldkc = kc.shape[1]
    else:
        return -7
    info = c_int(0)
    dllpde._writevtkug(c_char_p(fname.encode()), c_int(n), x.ctypes.data_as(POINTER(c_double)), y.ctypes.data_as(POINTER(c_double)), z.ctypes.data_as(POINTER(c_double)), c_int(ne), c_int(ldkc), kc.ctypes.data_as(POINTER(c_int)), u.ctypes.data_as(POINTER(c_double)), byref(info))
    return info.value
def readgmsh22(fname, x, y, z, kc, lb):
    dllpde._readgmsh22.argtype = [POINTER(c_char), POINTER(c_int), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_int), c_int, POINTER(c_int), c_int, POINTER(c_int), POINTER(c_int)]
    n = c_int(0)
    ne = c_int(0)
    if kc.ndim == 2:
        ldkc = kc.shape[1]
    else:
        return -5
    if lb.ndim == 2:
        ldlb = lb.shape[1]
    else:
        return -6
    info = c_int(0)
    dllpde._readgmsh22(c_char_p(fname.encode()), byref(n), x.ctypes.data_as(POINTER(c_double)), y.ctypes.data_as(POINTER(c_double)), z.ctypes.data_as(POINTER(c_double)), byref(ne), c_int(ldkc), kc.ctypes.data_as(POINTER(c_int)), c_int(ldlb), lb.ctypes.data_as(POINTER(c_int)), byref(info))
    if info.value == -9:
        info.value = -8
    return (info.value, n.value, ne.value)
def writegmsh22(fname, n, x, y, z, ne, kc, lb):
    dllpde._writegmsh22.argtype = [POINTER(c_char), c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), c_int, c_int, POINTER(c_int), c_int, POINTER(c_int), POINTER(c_int)]
    if kc.ndim == 2:
        ldkc = kc.shape[1]
    else:
        return -7
    if lb.ndim == 2:
        ldlb = lb.shape[1]
    else:
        return -8
    info = c_int(0)
    dllpde._writegmsh22(c_char_p(fname.encode()), c_int(n), x.ctypes.data_as(POINTER(c_double)), y.ctypes.data_as(POINTER(c_double)), z.ctypes.data_as(POINTER(c_double)), c_int(ne), c_int(ldkc), kc.ctypes.data_as(POINTER(c_int)), c_int(ldlb), lb.ctypes.data_as(POINTER(c_int)), byref(info))
    if info.value == -9:
        info.value = -8
    return info.value
def csr_check(m, n, val, rowptr, colind, result):
    dllsp._csr_check.argtype = [c_int, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_int), POINTER(c_int)]
    dllsp._csx_check_sym.argtype = [c_int, c_int, POINTER(c_double), POINTER(c_int), POINTER(c_int), POINTER(c_int)]
    info = c_int(0)
    dllsp._csr_check(c_int(m), c_int(n), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), result.ctypes.data_as(POINTER(c_int)), byref(info))
    if info.value == 0:
        result[10] = 0
        if m == n:
            info1 = c_int(0)
            dllsp._csx_check_sym(c_int(n), val.ctypes.data_as(POINTER(c_double)), rowptr.ctypes.data_as(POINTER(c_int)), colind.ctypes.data_as(POINTER(c_int)), byref(info1))
            if info1.value == 1:
                result[10] = 1
    return info.value
