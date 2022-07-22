# *********************************************
# *                                           *
# *  XLPack Lite Julia Library Test Program   *
# *  Version 6.0 (February 12, 2022)          *
# *  (C) 2014-2022  K Technologies            *
# *                                           *
# *********************************************/

include("XLPack.jl")

using .XLPack

function bj0(x)
    # besj0 by Chebyshev series approximation (0 <= x <= 4)
    bj0cs = [
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
        0.0000000000000000074 ]
    chebs(bj0cs, size(bj0cs, 1), 0.125*x*x - 1)
end

function TestSF()
    println("** d1num **")
    for i = 1 : 4
        println(i, ": ", d1num(i))
    end
    println("** chebs **")
    println(bj0(0.5), "  ", besj0(0.5))
    println("** powm1 **")
    println(powm1(2.0, 2.0e-13), "  ", 2.0^2.0e-13 - 1)
    println("** sqrt1pm1 **")
    println(sqrt1pm1(2.0e-13), "  ", sqrt(1.0 + 2.0e-13) - 1)
    println("** cos_pi **")
    println(cos_pi(9000.5), "  ", cos(9000.5*pi))
    println("** sin_pi **")
    println(sin_pi(9000.0), "  ", sin(9000.0*pi))
    println("** li **")
    println(li(2.0))
    println("** ei **")
    println(ei(1.0))
    println("** e1 **")
    println(e1(1.0))
    println("** digamma **")
    println(digamma(3.5))
    println("** besj0 **")
    println(besj0(1.0))
    println("** besj1 **")
    println(besj1(1.0))
    println("** besjnu **")
    println(besjnu(1.0, 1.0))
    println("** besy0 **")
    println(besy0(1.0))
    println("** besy1 **")
    println(besy1(1.0))
    println("** besynu **")
    println(besynu(1.0, 1.0))
    println("** besi0 **")
    println(besi0(1.0))
    println(besi0(1.0, true))
    println("** besi1 **")
    println(besi1(1.0))
    println(besi1(1.0, true))
    println("** besinu **")
    println(besinu(1.0, 1.0))
    println("** besk0 **")
    println(besk0(1.0))
    println(besk0(1.0, true))
    println("** besk1 **")
    println(besk1(1.0))
    println(besk1(1.0, true))
    println("** besknu **")
    println(besknu(1.0, 1.0))
    println("** celli1 **")
    println(celli1(0.5))
    println("** celli2 **")
    println(celli2(0.5))
    println("** celli3 **")
    println(celli3(0.7, 0.5))
end

function TestDconst()
    println("** dconst **")
    for i = 1 : 35
        println(i, ": ", dconst(i))
    end
end

function TestDgesv()
    println("** dgesv **")
    n = 3
    a = [ 0.2 -0.11 -0.93;
         -0.32 0.81 0.37;
         -0.8 -0.92 -0.29 ]
    b = [ -0.3727, 0.4319, -1.4247 ]
    ipiv = Vector{Cint}(undef, n)
    anorm, info = dlange('1', n, n, a)
    info = dgesv(n, a, ipiv, b)
    println("x = ", b, ", info = ", info)
    if info == 0
        rcond, info = dgecon('1', n, a, anorm)
        println("rcond = ", rcond, ", info = ", info)
    end
end

function TestDposv()
    println("** dposv **")
    n = 3
    a = [ 2.2 -0.11 -0.32;
          0.0 2.93 0.81;
          0.0 0.0 2.37 ]
    b = [ -1.566, -2.8425, -1.1765 ]
    anorm, info = dlansy('1', 'U', n, a)
    info = dposv('U', n, a, b)
    println("x = ", b, ", info = ", info)
    if info == 0
        rcond, info = dpocon('U', n, a, anorm)
        println("rcond = ", rcond, ", info = ", info)
    end
end

function TestDsyev()
    println("** dsyev **")
    n = 3
    a = [ 2.2 -0.11 -0.32;
          0.0 2.93 0.81;
          0.0 0.0 2.37 ]
    w = Vector{Cdouble}(undef, n)
    info = dsyev('V', 'U', n, a, w)
    println("w = ", w, ", info = ", info)
    println("a:")
    println(a)
end

function TestDgels()
    println("** dgels (2) **")
    m = 4
    n = 3
    a = [ -1.06  0.48 -0.04;
          -1.19  0.73 -0.24;
           1.97 -0.89  0.56;
           0.68 -0.53  0.08 ]
    b = [ 0.3884, 0.1120, -0.3644, -0.0002 ]
    info = dgels('N', m, n, a, b)
    println(b[1], " ", b[2], " ", b[3], " ", info)
    if info == 0
        ci = Vector{Cdouble}(undef, n)
        info = dgecov(0, n, a, ci)
        println(ci, " ", info)
    end
end

function TestPchse()
    println("** pchse **")
    n = 4
    x = [ 0.1, 0.11, 0.12, 0.13 ]
    y = [ -2.3026, -2.2073, -2.1203, -2.0402 ]
    d = Vector{Cdouble}(undef, n)
    info = pchse(n, x, y, d)
    println(info)
    if info == 0
        ne = 2
        xe = [ 0.115, 0.125 ]
        ye = Vector{Cdouble}(undef, ne)
        info = pchfe(n, x, y, d, ne, xe, ye)
        println(ye)
        println(info)
    end
end

function TestPchia()
    println("** pchia **")
    n = 7
    x = [ -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 ]
    y = [ 0.5, 1, 0.5, 0.2, 0.1, 0.05882, 0.03846 ]
    d = Vector{Cdouble}(undef, n)
    info = pchse(n, x, y, d)
    println(info)
    if info == 0
        a = 0.0
        b = 4.0
        s, info = pchia(n, x, y, d, a, b)
        println(s, " ", info)
    end
end

function TestRpzero2()
    println("** rpzero2 **")
    n = 5
    a = [ 1.0, 0.0, 2.0, 2.0, -15.0, 10.0 ]
    rr = Vector{Cdouble}(undef, n)
    ri = Vector{Cdouble}(undef, n)
    s = Vector{Cdouble}(undef, n)
    iter, info = rpzero2(n, a, rr, ri, s)
    for i = 1:n
        println(rr[i], " ", ri[i], " ", s[i])
    end
    println("iter = ", iter, ", info = ", info)
end

function TestDfzero()
    println("** dfzero **")
    f(x) = (x*x - 2.0)*x - 5.0
    b = 1.0
    c = 3.0
    r = b
    println("b = ", b, ", c = ", c)
    x, info = dfzero(f, b, c, r)
    println("x = ", x, ", info = ", info)
end

function fhybrd1(n, x, fvec, iflag)
    fvec[1] = x[1]^2 - x[2] - 1
    fvec[2] = (x[1] - 2)^2 + (x[2] - 0.5)^2 - 1
end

function TestHybrd1()
    println("** hybrd1 **")
    n = 2
    x = [ 0.0, 0.0 ]
    fvec = Vector{Cdouble}(undef, n)
    info = hybrd1(fhybrd1, n, x, fvec)
    println(x)
    println(fvec)
    println(info)
end

function TestDfmin()
    println("** dfmin **")
    f(x) = x^3 - 2*x - 5
    a = 0.0
    b = 2.0
    x = dfmin(a, b, f)
    println(x)
end

function TestOptif0()
    println("** optif0 **")
    f(n, x) = 100*(x[2] - x[1]^2)^2 + (1 - x[1])^2
    n = 2
    x = [ -1.2, 1.0 ]
    xpls = Vector{Cdouble}(undef, n)
    fpls, info = optif0(n, x, f, xpls)
    println(xpls)
    println(fpls)
    println(info)
end

function TestQk15()
    println("** qk15 **")
    f(x) = 1/(1 + x^2)
    a = 0
    b = 4
    s, abserr = qk15(f, a, b)
    println(s, " ", abserr)
end

function TestQag()
    println("** qag **")
    f(x) = 1/(1 + x^2)
    a = 0
    b = 4
    s, abserr, info = qag(f, a, b)
    println(s, " ", abserr, " ", info)
end

function TestQagi()
    println("** qagi **")
    f(x) = 1/(1 + x^2)
    bound = 0.0
    inf = 2
    s, abserr, info = qagi(f, bound, inf)
    println("qagi [-inf, +inf]")
    println(s, " ", abserr, " ", info)
    inf = 1
    s, abserr, info = qagi(f, bound, inf)
    println("qagi [0, +inf]")
    println(s, " ", abserr, " ", info)
    inf = -1
    s, abserr, info = qagi(f, bound, inf)
    println("qagi [-inf, 0]")
    println(s, " ", abserr, " ", info)
end

function TestDerkf()
    println("** derkf **")
    f(n, t, y, yp) = begin
        yp[1] = -2*y[1] + y[2] - cos(t)
        yp[2] = 2*y[1] - 3*y[2] + 3*cos(t) - sin(t)
    end
    n = 2
    wsave = Vector{Cdouble}(undef, 7*n + 20)
    iwsave = Vector{Cint}(undef, 20)
    t = 0
    y = [ 1.0, 2.0 ]
    tfinal = 10.0
    tprint = 1.0
    info = 0
    while true
        tout = t + tprint
        t, info = derkf(info, n, f, t, y, tout, wsave, iwsave)
        println(t, " ", y)
        if t >= tfinal || info != 1
            break
        end
    end
end

function TestDerkf_2()
    println("** derkf (2) (dense output) **")
    f(n, t, y, yp) = begin
        yp[1] = -2*y[1] + y[2] - cos(t)
        yp[2] = 2*y[1] - 3*y[2] + 3*cos(t) - sin(t)
    end
    n = 2
    wsave = Vector{Cdouble}(undef, 9*n + 20)
    iwsave = Vector{Cint}(undef, 20)
    y_temp = Vector{Cdouble}(undef, n)
    t = 0
    y = [ 1.0, 2.0 ]
    tfinal = 10.0
    tprint = 1.0
    tout = tprint
    info = 0
    while true
        t, info = derkf(info, n, f, t, y, tfinal, wsave, iwsave, mode = 2)
        if info != 1 && info != 2
            print("Error: info =", info)
            break
        end
        while t >= tout
            derkf_int(n, tout, y_temp, wsave)
            println(tout, " ", y_temp)
            tout = tout + tprint
        end
        if t >= tfinal
            break
        end
    end
end

function TestRfft1()
    # Initialization
    println("** Rfft1 **")
    seed = 13
    init_genrand(seed)
    n = 5
    lwsave::Cint = n + round(Cint, log(n)/log(2)) + 4
    wsave = Vector{Cdouble}(undef, lwsave)
    info = rfft1i(n, wsave)
    if info != 0
        println("Error during initialization: info = ", info)
        return 1
    end
    # Generate test data
    inc = 1
    lr = inc*(n - 1) + 1
    r = Vector{Cdouble}(undef, lr)
    rcopy = Vector{Cdouble}(undef, lr)
    for i = 1:n
        k = 1 + inc*(i - 1)
        r[k] = genrand_res53()
        rcopy[k] = r[k]
    end
    # Forward transform
    info = rfft1f(n, r, wsave)
    if info != 0
        println("Error in Rfft1f: info = ", info)
        return 2
    end
    # Backward transform
    info = rfft1b(n, r, wsave)
    if info != 0
        println("Error in Rfft1b: info = ", info)
        return 3
    end
    # Check results
    diff = 0.0
    for i = 1:n
        k = 1 + inc*(i - 1)
        if abs(r[k] - rcopy[k]) > diff
            diff = abs(r[k] - rcopy[k])
        end
        println(rcopy[k], " ", r[k], " ", abs(r[k] - rcopy[k]))
    end
    println("diff(max) = ", diff)
end

function flmdif1(m, n, x, fvec, iflag)
    xdata = [ 77.6, 239.9, 434.8, 760.0 ]
    ydata = [ 10.07, 29.61, 50.76, 81.78 ]
    for i = 1 : m
        fvec[i] = ydata[i] - x[1]*(1 - exp(-xdata[i]*x[2]))
    end
end

function TestLmdif1()
    println("** lmdif1 (2) **")
    m = 4
    n = 2
    x = [ 500.0, 0.0001 ]
    fvec = Vector{Cdouble}(undef, m)
    info = lmdif1(flmdif1, m, n, x, fvec)
    println(x)
    println(fvec)
    println(info)
end

function TestGenrandInt31()
    println("** genrand_int31 **")
    seed = 13
    println("seed = ", seed)
    init_genrand(seed)
    for i = 1 : 10
        println(genrand_int31())
    end
end

function TestGenrandInt32()
    println("** genrand_int32 **")
    seed = 13
    println("seed = ", seed)
    init_genrand(seed)
    for i = 1 : 10
        println(genrand_int32())
    end
end

function TestGenrandRes53()
    println("** genrand_res53 **")
    seed = 13
    println("seed = ", seed)
    init_genrand(seed)
    for i = 1 : 10
        println(genrand_res53())
    end
end

function TestDlamch()
    println("** dlamch **")
    for c in [ 'e', 's', 'b', 'p', 'n', 'r', 'm', 'u', 'l', 'o']
        println(c, ": ", dlamch(c))
    end
end

#---

TestSF()
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
