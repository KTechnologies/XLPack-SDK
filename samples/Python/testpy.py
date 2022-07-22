# *********************************************
# *                                           *
# *  Experimental Python interface to XLPack  *
# *  Test program                             *
# *  Version 6.0 (June 14, 2022)              *
# *  (C) 2014-2022  K Technologies            *
# *                                           *
# *********************************************

import numpy as np
from math import *
import platform
version = platform.python_version_tuple()
if platform.system() == 'Windows' and int(version[0]) >= 3 and int(version[1]) >= 8:
    # The following statements are used to run with python 3.8 or later
    import os
    # Set XLPack.dll install directory
    os.add_dll_directory(os.path.expanduser('~/AppData/Local/Programs/XLPack'))
import struct
if struct.calcsize('P') == 8:
    from XLPack import *
else:
    from XLPack_32 import *

def bj0(x):
    # besj0 by Chebyshev series approximation (0 <= x <= 4)
    bj0cs = np.array([
        0.100254161968939137,
        -0.665223007764405132,
        0.248983703498281314,
        -0.0332527231700357697,
        0.0023114179304694015,
        -0.0000991127741995080,
        0.0000028916708643998,
        -0.0000000612108586630,
        0.0000000009838650793,
        -0.0000000000124235515,
         0.0000000000001265433,
        -0.0000000000000010619,
        0.0000000000000000074])
    return chebs(bj0cs, bj0cs.size, 0.125*x*x - 1)

def TestSf():
    print('** d1num **')
    for i in range(5):
        print(i, d1num(i))
    print('** factorial **')
    print(factorial(10))
    print('** chebs **')
    print(bj0(0.5), besj0(0.5))
    print('** powm1 **')
    print(powm1(2.0, 2.0e-13), 2.0**2.0e-13 - 1)
    print('** sqrt1pm1 **')
    print(sqrt1pm1(2.0e-13), sqrt(1.0 + 2.0e-13) - 1)
    print('** cos_pi **')
    print(cos_pi(9000.5), cos(9000.5*np.pi))
    print('** sin_pi **')
    print(sin_pi(9000.0), sin(9000.0*np.pi))
    print('** li **')
    print(li(2.0))
    print('** ei **')
    print(ei(1.0))
    print('** e1 **')
    print(e1(1.0))
    print('** digamma **')
    print(digamma(3.5))
    print('** besj0 **')
    print(besj0(1.0))
    print('** besj1 **')
    print(besj1(1.0))
    print('** besjnu **')
    print(besjnu(1.0, 1.0))
    print('** besy0 **')
    print(besy0(1.0))
    print('** besy1 **')
    print(besy1(1.0))
    print('** besynu **')
    print(besynu(1.0, 1.0))
    print('** besi0 **')
    print(besi0(1.0))
    print(besi0(1.0, 2))
    print('** besi1 **')
    print(besi1(1.0))
    print(besi1(1.0, 2))
    print('** besinu **')
    print(besinu(1.0, 1.0))
    print('** besk0 **')
    print(besk0(1.0))
    print(besk0(1.0, 2))
    print('** besk1 **')
    print(besk1(1.0))
    print(besk1(1.0, 2))
    print('** besknu **')
    print(besknu(1.0, 1.0))
    print('** celli1 **')
    print(celli1(0.5))
    print('** celli2 **')
    print(celli2(0.5))
    print('** celli3 **')
    print(celli3(0.7, 0.5))
    # Check errno
    print('** celli3 ** (check errno)')
    set_errno(0)
    print(celli3(0.5, 0.5))
    print(get_errno())
    set_errno(0)
    print(celli3(0.5, 1.0))
    print(get_errno())

def TestDconst():
    print('** Dconst **')
    for i in range(36):
        print(dconst(i))

def TestDgesv():
    print ('** dgesv **')
    n = 3
    a = np.array([
        [0.2, -0.32, -0.8],
        [-0.11, 0.81, -0.92],
        [-0.93, 0.37, -0.29]])
    b = np.array([-0.3727, 0.4319, -1.4247])
    ipiv = np.empty(n, dtype=int)
    anorm, info = dlange('1', n, n, a)
    info = dgesv(n, a, ipiv, b)
    print(b, info)
    rcond, info = dgecon('1', n, a, anorm)
    print(rcond, info)

def TestDposv():
    print ('** dposv **')
    n = 3
    a = np.array([
        [2.2, 0.0, 0.0],
        [-0.11, 2.93, 0.0],
        [-0.32, 0.81, 2.37]])
    b = np.array([-1.566, -2.8425, -1.1765])
    anorm, info = dlansy('1', 'U', n, a)
    info = dposv('U', n, a, b)
    print(b, info)
    rcond, info = dpocon('U', n, a, anorm)
    print(rcond, info)

def TestDsyev():
    print ('** dsyev **')
    n = 3
    a = np.array([
        [2.2, 0.0, 0.0],
        [-0.11, 2.93, 0.0],
        [-0.32, 0.81, 2.37]])
    w = np.empty(3)
    info = dsyev('V', 'U', n, a, w)
    print(w, info)
    print(a)

def TestDgels():
    print ('** dgels **')
    m = 4
    n = 2
    x = np.array([0.2, 118.2, 337.4, 884.6])
    y = np.array([0.1, 118.1, 338.8, 888.0])
    a = np.empty((n, m))
    for i in range(m):
        a[0][i] = 1
        a[1][i] = x[i]
    info = dgels('N', m, n, a, y)
    print(y[0], y[1], info)
    if info == 0:
        ci = np.empty(n)
        info = dgecov(0, n, a, ci)
        s = 0.0
        for i in range(n, m):
            s = s + y[i]**2
        s = s/(m - n)
        print(sqrt(s*ci[0]), sqrt(s*ci[1]), info)

def TestPchse():
    print ('** pchse **')
    n = 4
    x = np.array([0.1, 0.11, 0.12, 0.13])
    y = np.array([-2.3026, -2.2073, -2.1203, -2.0402])
    d = np.empty(n)
    info = pchse(n, x, y, d)
    print(info)
    if info == 0:
        ne = 2
        xe = np.array([0.115, 0.125])
        ye = np.empty(ne)
        info = pchfe(n, x, y, d, ne, xe, ye)
        print(ye, info)

def TestPchia():
    print ('** pchia **')
    n = 7
    x = np.array([-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    y = np.array([0.5, 1, 0.5, 0.2, 0.1, 0.05882, 0.03846])
    d = np.empty(n)
    info = pchse(n, x, y, d)
    print(info)
    if info == 0:
        a = 0.0
        b = 4.0
        s, info = pchia(n, x, y, d, a, b)
        print(s, info)

def TestRpzero2():
    print ('** rpzero2 **')
    n = 5;
    a = np.array([1.0, 0.0, 2.0, 2.0, -15.0, 10.0])
    rr = np.empty(n)
    ri = np.empty(n)
    s = np.empty(n)
    iter, info = rpzero2(n, a, rr, ri, s)
    for i in range(n):
        print(rr[i], ri[i], s[i])
    print(iter, info)

def FDfzero(x):
    return (x*x - 2.0)*x - 5.0

def TestDfzero():
    print ('** dfzero **')
    b = 1.0
    c = 3.0
    r = b
    print(b, c)
    x, info = dfzero(FDfzero, b, c, r)
    print(x, info)

def FHybrd1(n, x, fvec, iflag):
    if iflag == 1 or iflag == 2:
        fvec[0] = x[0]**2 - x[1] - 1
        fvec[1] = (x[0] - 2)**2 + (x[1] - 0.5)**2 - 1

def TestHybrd1():
    print ('** hybrd1 **')
    n = 2
    x = np.array([0.0, 0.0])
    fvec = np.empty(n)
    info = hybrd1(FHybrd1, n, x, fvec)
    print(x)
    print(fvec)
    print(info)

def FDfmin(x):
    return x**3 - 2*x - 5

def TestDfmin():
    print ('** dfmin **')
    a = 0.0
    b = 2.0
    x, info = dfmin(a, b, FDfmin)
    print(x, info)

def FOptif0(n, x):
    return 100*(x[1] - x[0]**2)**2 + (1 - x[0])**2

def TestOptif0():
    print ('** optif0 **')
    n = 2
    x = np.array([-1.2, 1.0])
    xpls = np.empty(n)
    fpls, info = optif0(n, x, FOptif0, xpls)
    print(xpls)
    print(fpls)
    print(info)

def FQk15(x):
    return 1/(1 + x**2)

def TestQk15():
    print ('** qk15 **')
    a = 0
    b = 4
    s, abserr = qk15(FQk15, a, b)
    print(s, abserr)

def FQag(x):
    return 1/(1 + x**2)

def TestQag():
    print ('** qag **')
    a = 0
    b = 4
    s, abserr, info = qag(FQag, a, b)
    print(s, abserr, info)

def FQagi(x):
    return 1/(1 + x**2)

def TestQagi():
    print ('** qagi **')
    bound = 0.0
    inf = 2
    s, abserr, info = qagi(FQagi, bound, inf)
    print('qagi [-inf, +inf]')
    print(s, abserr, info)
    inf = 1
    s, abserr, info = qagi(FQagi, bound, inf)
    print('qagi [0, +inf]')
    print(s, abserr, info)
    inf = -1
    s, abserr, info = qagi(FQagi, bound, inf)
    print('qagi [-inf, 0]')
    print(s, abserr, info)

def FDerkf(n, t, y, yp):
    yp[0] = -2*y[0] + y[1] - cos(t)
    yp[1] = 2*y[0] - 3*y[1] + 3*cos(t) - sin(t)

def TestDerkf():
    print ('** Derkf **')
    n = 2
    wsave = np.empty(7*n + 20)
    iwsave = np.empty(20, dtype=int)
    t = 0
    y = np.array([1.0, 2.0])
    tfinal = 10.0
    tprint = 1.0
    info = 0
    while True:
        tout = t + tprint
        t, info = derkf(info, n, FDerkf, t, y, tout, wsave, iwsave)
        print(t, y)
        if t >= tfinal or info != 1:
            break
    if info != 1:
        print('Error: info =', info)

def TestDerkf_2():
    print ('** Derkf (2) (dense output) **')
    n = 2
    wsave = np.empty(9*n + 20)
    iwsave = np.empty(20, dtype=int)
    y_temp = np.empty(n)
    t = 0
    y = np.array([1.0, 2.0])
    tfinal = 10.0
    tprint = 1.0
    tout = tprint
    info = 0
    while True:
        t, info = derkf(info, n, FDerkf, t, y, tfinal, wsave, iwsave, mode = 2)
        if info != 1 and info != 2:
            print('Error: info =', info)
            break
        while t >= tout:
            derkf_int(n, tout, y_temp, wsave)
            print(tout, y_temp)
            tout = tout + tprint
        if t >= tfinal:
            break

def TestRfft1():
    # Initialization
    print('** Rfft1 **')
    seed = 13
    init_genrand(seed)
    n = 5
    lwsave = n + int(log(n)/log(2)) + 4
    wsave = np.empty(lwsave)
    info = rfft1i(n, wsave)
    if info != 0:
        print('Error during initialization: info =', info)
        return 1
    # Generate test data
    inc = 1
    lr = inc*(n - 1) + 1
    r = np.empty(lr)
    rcopy = np.empty(lr)
    for i in range(n):
        k = inc*i
        r[k] = genrand_res53()
        rcopy[k] = r[k]
    # Forward transform
    info = rfft1f(n, r, wsave)
    if info != 0:
        print('Error in Rfft1f: info =', info)
        return 2
    # Backward transform
    info = rfft1b(n, r, wsave)
    if info != 0:
        print('Error in Rfft1b: info =', info)
        return 3
    # Check results
    diff = 0.0
    for i in range(n):
        k = inc*i
        if abs(r[k] - rcopy[k]) > diff:
            diff = abs(r[k] - rcopy[k])
        print(rcopy[k], r[k], abs(r[k] - rcopy[k]))
    print('diff(max) =', diff)

# Data for TestLmdif1
xdata = np.array([77.6, 114.9, 141.1, 190.8, 239.9, 289.0, 332.8, 378.4,
    434.8, 477.3, 536.8, 593.1, 689.1, 760.0])
ydata = np.array([10.07, 14.73, 17.94, 23.93, 29.61, 35.18, 40.02, 44.82,
    50.76, 55.05, 61.01, 66.4, 75.47, 81.78])

def FLmdif1(m, n, x, fvec, iflag):
    if iflag == 1 or iflag == 2:
        for i in range(m):
            fvec[i] = ydata[i] - x[0]*(1 - exp(-xdata[i]*x[1]))

def TestLmdif1():
    print ('** lmdif1 **')
    m = 14
    n = 2
    x = np.array([500.0, 0.0001])
    fvec = np.empty(m)
    info = lmdif1(FLmdif1, m, n, x, fvec)
    print(x)
    print(fvec)
    print(info)

def TestGenrandInt31():
    print('** genrand_int31 **')
    seed = 13
    print('seed =', seed)
    init_genrand(seed)
    for i in range(10):
        print(genrand_int31())

def TestGenrandInt32():
    print('** genrand_int32 **')
    seed = 13
    print('seed =', seed)
    init_genrand(seed)
    for i in range(10):
        print(genrand_int32())

def TestGenrandRes53():
    print('** genrand_res53 **')
    seed = 13
    print('seed =', seed)
    init_genrand(seed)
    for i in range(10):
        print(genrand_res53())

def TestDlamch():
    print ('** dlamch **')
    for c in 'esbpnrmulo':
        print (c, dlamch(c))

# -----------------

TestSf()
TestDconst()
TestDgesv()
TestDposv()
TestDsyev()
TestDgels()
TestPchse()
TestPchia()
TestRpzero2()
TestDfzero()
TestHybrd1()
TestDfmin()
TestOptif0()
TestQk15()
TestQag()
TestQagi()
TestDerkf()
TestDerkf_2()
TestRfft1()
TestLmdif1()
TestGenrandInt31()
TestGenrandInt32()
TestGenrandRes53()
TestDlamch()
