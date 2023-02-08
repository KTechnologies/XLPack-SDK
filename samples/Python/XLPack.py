# ****************************************
# *                                      *
# *  XLPack Python Numerical Library     *
# *  Version 6.1 (December 1, 2022)      *
# *  (C) 2014-2022  K Technologies       *
# *                                      *
# ****************************************/
import os
import struct
import numpy as np
from ctypes import *
if os.name == 'posix':
    dll = CDLL('/Library/Application Support/Microsoft/Office365/User Content.localized/Add-Ins.localized/XLPack.dylib')
elif struct.calcsize('P') == 8:
    dll = CDLL('XLPack.dll')
else:
    dll = CDLL('XLPack_32.dll')
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
def derkf(info, n, f, t, y, tout, wsave, iwsave, rtol = 1.0e-10, atol = 1.0e-10, mode = 0):
    dll._derkf_r.argtype = [c_int, POINTER(c_double), POINTER(c_double), c_double, POINTER(c_double), POINTER(c_double), c_int, c_int, POINTER(c_double), c_int, POINTER(c_int), c_int, POINTER(c_int), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_int)]
    info1 = c_int(info)
    t1 = c_double(t)
    if n < 1:
        return t, -2
    if y.ndim != 1 or y.shape[0] < n:
        return t, -5
    if wsave.ndim != 1:
        return t, -7
    lwsave = wsave.shape[0]
    if iwsave.ndim != 1:
        return t, -8
    liwsave = iwsave.shape[0]
    rtol1 = c_double(rtol)
    atol1 = c_double(atol)
    itol = c_int(0)
    tt = c_double(0)
    yy = np.empty(n)
    yyp = np.empty(n)
    irev = c_int(0)
    while True:
        dll._derkf_r(c_int(n), byref(t1), y.ctypes.data_as(POINTER(c_double)), c_double(tout), byref(rtol1), byref(atol1), itol, c_int(mode), wsave.ctypes.data_as(POINTER(c_double)), c_int(lwsave), iwsave.ctypes.data_as(POINTER(c_int)), c_int(liwsave), byref(info1), byref(tt), yy.ctypes.data_as(POINTER(c_double)), yyp.ctypes.data_as(POINTER(c_double)), byref(irev))
        if irev.value >= 1 and irev.value <= 11:
            f(n, tt.value, yy, yyp)
        elif irev.value == 0:
            break
    if info1.value == -13:
        info1.value = -1
    elif info1.value == -10:
        info1.value = -7
    elif info1.value == -12:
        info1.value = -8
    elif info1.value == -5:
        info1.value = -9
    elif info1.value == -6:
        info1.value = -10
    elif info1.value < 0:
        info1.value = info1.value - 2
    return t1.value, info1.value
def derkf_int(n, t, y, wsave):
    dll._derkf_int.argtype = [c_int, c_double, POINTER(c_double), POINTER(c_double)]
    if n < 1 or y.ndim != 1 or wsave.ndim != 1:
        return
    dll._derkf_int(c_int(n), c_double(t), y.ctypes.data_as(POINTER(c_double)), wsave.ctypes.data_as(POINTER(c_double)))
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
