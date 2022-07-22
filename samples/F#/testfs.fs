(*****************************************
 *                                       *
 *  Experimental F# interface to XLPack  *
 *  Test program                         *
 *  Version 5.4 (November 20, 2020)      *
 *  (C) 2014-2020  K Technologies        *
 *                                       *
 *****************************************)

open XLPack

let TestSf() =

    let PrintFunc fmt arg result =
        match result with
            | (value, errno) -> printfn fmt arg value errno

    let PrintFunc2 fmt arg1 arg2 result =
        match result with
            | (value, errno) -> printfn fmt arg1 arg2 value errno

    let mutable result = 0.0, 0
    let mutable ix, x, x2, nu = 0, 0.0, 0.0, 0.0

    printfn "** D1num"
    for i = 1 to 4 do
        printfn "%d: %g" i (XLPack.D1num(i))

    printfn "** Factorial"
    ix <- 10
    result <- XLPack.Factorial(ix)
    PrintFunc "Factorial(%d) = %g, errno = %d" ix result

    printfn "** Li"
    x <- 2.0
    result <- XLPack.Li(x)
    PrintFunc "Li(%g) = %g, errno = %d" x result

    printfn "** Ei"
    x <- 1.0
    result <- XLPack.Ei(x)
    PrintFunc "Ei(%g) = %g, errno = %d" x result

    printfn "** E1"
    x <- 1.0
    result <- XLPack.E1(x)
    PrintFunc "E1(%g) = %g, errno = %d" x result

    printfn "** Digamma"
    x <- 3.5
    result <- XLPack.Digamma(x)
    PrintFunc "Digamma(%g) = %g, errno =%d" x result

    printfn "** Besj0"
    x <- 1.0
    result <- XLPack.Besj0(x)
    PrintFunc "Besj0(%g) = %g, errno = %d" x result

    printfn "** Besj1"
    x <- 1.0
    result <- XLPack.Besj1(x)
    PrintFunc "Besj1(%g) = %g, errno = %d" x result

    printfn "** Besjnu"
    nu <- 1.0; x <- 1.0
    result <- XLPack.Besjnu(nu, x)
    PrintFunc2 "Besjnu(%g, %g) = %g, errno = %d" nu x result

    printfn "** Besy0"
    x <- 1.0
    result <- XLPack.Besy0(x)
    PrintFunc "Besy0(%g) = %g, errno = %d" x result

    printfn "** Besy1"
    x <- 1.0
    result <- XLPack.Besy1(x)
    PrintFunc "Besy1(%g) = %g, errno = %d" x result

    printfn "** Besynu"
    nu <- 1.0; x <- 1.0
    result <- XLPack.Besynu(nu, x)
    PrintFunc2 "Besynu(%g, %g) = %g, errno = %d" nu x result

    printfn "** Besi0"
    x <- 1.0
    result <- XLPack.Besi0(x)
    PrintFunc "Besi0(%g) = %g, errno = %d" x result
    result <- XLPack.Besi0(x, 2)
    PrintFunc "Besi0(%g, 2) = %g, errno = %d" x result

    printfn "** Besi1"
    x <- 1.0
    result <- XLPack.Besi1(x)
    PrintFunc "Besi1(%g) = %g, errno = %d" x result
    result <- XLPack.Besi1(x, 2)
    PrintFunc "Besi1(%g, 2) = %g, errno = %d" x result

    printfn "** Besinu"
    nu <- 1.0; x <- 1.0
    result <- XLPack.Besinu(nu, x)
    PrintFunc2 "Besinu(%g, %g) = %g, errno = %d" nu x result

    printfn "** Besk0"
    x <- 1.0
    result <- XLPack.Besk0(x)
    PrintFunc "Besk0(%g) = %g, errno = %d" x result
    result <- XLPack.Besk0(x, 2)
    PrintFunc "Besk0(%g, 2) = %g, errno = %d" x result

    printfn "** Besk1"
    x <- 1.0
    result <- XLPack.Besk1(x)
    PrintFunc "Besk1(%g) = %g, errno = %d" x result
    result <- XLPack.Besk1(x, 2)
    PrintFunc "Besk1(%g, 2) = %g, errno = %d" x result

    printfn "** Besknu"
    nu <- 1.0; x <- 1.0
    result <- XLPack.Besknu(nu, x)
    PrintFunc2 "Besknu(%g, %g) = %g, errno = %d" nu x result

    printfn "** Celli1"
    x <- 0.5
    result <- XLPack.Celli1(x)
    PrintFunc "Celli1(%g) = %g, errno = %d" x result

    printfn "** Celli2"
    x <- 0.5
    result <- XLPack.Celli2(x)
    PrintFunc "Celli2(%g) = %g, errno = %d" x result

    printfn "** Celli3"
    x <- 0.7; x2 <- 0.5
    result <- XLPack.Celli3(x, x2)
    PrintFunc2 "Celli3(%g, %g) = %g, errno = %d" x x2 result

    printfn "** Dconst"
    for i = 0 to 35 do
        printfn "%d: %g" i (XLPack.Dconst(i))

let TestDgesv() =
    let n = 3
    let a = array2D [
        [ 0.2; -0.32; -0.8 ];
        [ -0.11; 0.81; -0.92 ];
        [ -0.93; 0.37; -0.29 ] ]
    let b = [| -0.3727; 0.4319; -1.4247 |]
    let ipiv = Array.create n 0
    let anorm = XLPack.Dlange('1', n, n, a)
    let info = XLPack.Dgesv(n, a, ipiv, b)
    printfn "** Dgesv"
    printfn "x = %A, info = %d" b info
    if info = 0 then
        let rcond, info2 = XLPack.Dgecon('1', n, a, anorm)
        printfn "rcond = %f, info = %d" rcond info2

let TestDposv() =
    let n = 3
    let a = array2D [
        [ 2.2; 0.0; 0.0 ];
        [ -0.11; 2.93; 0.0 ];
        [ -0.32; 0.81; 2.37 ] ]
    let b = [| -1.566; -2.8425; -1.1765 |]
    let anorm = XLPack.Dlansy('1', 'U', n, a)
    let info = XLPack.Dposv('U', n, a, b)
    printfn "** Dposv"
    printfn "x = %A, info = %d" b info
    if info = 0 then
        let rcond, info2 = XLPack.Dpocon('U', n, a, anorm)
        printfn "rcond = %f, info = %d" rcond info2

let TestDsyev() =
    let n = 3
    let a = array2D [
        [ 2.2; 0.0; 0.0 ];
        [ -0.11; 2.93; 0.0 ];
        [ -0.32; 0.81; 2.37 ] ]
    let w = Array.create n 0.0
    let info = XLPack.Dsyev('V', 'U', n, a, w)
    printfn "** Dsyev"
    printfn "Eigenvalues = %A" w
    printfn "Eigenvectors = %A" a
    printfn "info = %d" info

let TestDgels() =
    let m = 4
    let n = 2
    let x = [| 0.2; 118.2; 337.4; 884.6 |]
    let y = [| 0.1; 118.1; 338.8; 888.0 |]
    let a = Array2D.create n m 1.0
    a.[1, *] <- x.[*]
    let info = XLPack.Dgels('N', m, n, a, y)
    printfn "** Dgels"
    printfn "a0 = %f, a1 = %f" y.[0] y.[1]
    printfn "info = %d" info
    if info = 0 then
        let job = 0
        let ci = Array.create n 0.0
        let info2 = XLPack.Dgecov(job, n, a, ci)
        let sumsq = Array.fold (fun s y -> s + y*y) 0.0 y.[n..m-1]
        let s = sumsq/(double (m - n))
        let stddev = Array.map (fun x -> sqrt(s*x)) ci
        printfn "std. dev. = %A" stddev
        printfn "info = %d" info2

let TestPchse() =
    let n = 4
    let ne = 2
    let x = [| 0.1; 0.11; 0.12; 0.13 |]
    let y = [| 2.3026; 2.2073; 2.1203; 2.0402 |]
    let d = Array.create n 0.0
    let xe = [| 0.115; 0.125 |]
    let ye = Array.create ne 0.0
    let info = XLPack.Pchse(n, x, y, d)
    printfn "** Pchse"
    printfn "info = %d" info
    if info = 0 then
        let info2 = XLPack.Pchfe(n, x, y, d, ne, xe, ye)
        printfn "ln(%f) = %f" xe.[0] ye.[0]
        printfn "ln(%f) = %f" xe.[1] ye.[1]
        printfn "info = %d" info2

let TestPchia() =
    let n = 7
    let x = [| -1.0; 0.0; 1.0; 2.0; 3.0; 4.0; 5.0 |]
    let y = [| 0.5; 1.0; 0.5; 0.2; 0.1; 0.05882; 0.03846 |]
    let d = Array.create n 0.0
    let info = XLPack.Pchse(n, x, y, d)
    printfn "** Pchia"
    printfn "info = %d" info
    if info = 0 then
        let a, b = (0.0, 4.0)
        let s, info2 = XLPack.Pchia(n, x, y, d, a, b)
        printfn "s = %f, info = %d" s info2

let TestRpzero2() =
    let n = 5
    let a = [| 1.0; 0.0; 2.0; 2.0; -15.0; 10.0 |]
    let zr = Array.create n 0.0
    let zi = Array.create n 0.0
    let s = Array.create n 0.0
    let iter, info = XLPack.Rpzero2(n, a, zr, zi, s)
    printfn "** Rpzero2"
    for i = 0 to n - 1 do
        printfn "(%g, %g), %g" zr.[i] zi.[i] s.[i]
    printfn "iter = %d, info = %d" iter info

let TestDfzero() =
    let f(x: double) = (x*x - 2.0)*x - 5.0
    let mutable b, c = 1.0, 3.0
    let r = 1.0
    let info = XLPack.Dfzero(f, &b, &c, r)
    printfn "** Dfzero"
    printfn "x = %f, info = %d" b info

let TestHybrd1() =
    let f(n: int, x: double[], y: double[], iflag: int) =
        y.[0] <- 4.0*x.[0]*x.[0] + x.[1]*x.[1] - 16.0
        y.[1] <- x.[0]*x.[0] + x.[1]*x.[1] - 9.0
        0
    let n = 2
    let x = [| 1.0; 2.0 |]
    let fvec = Array.create 2 0.0
    let info = XLPack.Hybrd1(f, n, x, fvec)
    printfn "** Hybrd1"
    printfn "x = %A" x
    printfn "info = %d" info

let TestDfmin() =
    let f(x: double) = (x*x - 2.0)*x - 5.0
    let a, b = (0.0, 1.0)
    let x = XLPack.Dfmin(a, b, f)
    printfn "** Dfmin"
    printfn "x = %f" x

let TestOptif0() =
    let f(n: int, x: double[]) =
        100.0*(x.[1] - x.[0]*x.[0])*(x.[1] - x.[0]*x.[0]) + (1.0 - x.[0])*(1.0 - x.[0])
    let n = 2
    let x = [| -1.2; 1.0 |]
    let xpls = Array.create n 0.0
    let fpls, info = XLPack.Optif0(n, x, f, xpls)
    printfn "** Optif0"
    printfn "xpls = %A" xpls
    printfn "fpls = %g" fpls
    printfn "info = %d" info

let TestQk15() =
    let f(x: double) = 1.0/(1.0 + x*x)
    let a, b = (0.0, 4.0)
    let result, abserr = XLPack.Qk15(f, a, b)
    printfn "** Qk15"
    printfn "result = %f, abserr = %g" result abserr

let TestQag() =
    let f(x: double) = 1.0/(1.0 + x*x)
    let a, b = (0.0, 4.0)
    let epsabs, epsrel = 1.0e-10, 1.0e-10
    for key = 1 to 6 do
        let result, abserr, info = XLPack.Qag(f, a, b, epsabs, epsrel, key)
        printfn "** Qag (Key = %d)" key
        printfn "result = %f, abserr = %g, info = %d" result abserr info

let TestQagi() =
    let f(x: double) = 1.0/(1.0 + x*x)
    let bound = 0.0
    for inf = -1 to 1 do
        let result, abserr, info = XLPack.Qagi(f, bound, inf)
        if inf = 1 then printfn "** Qagi [0, +inf]"
        elif inf = -1 then printfn "** Qagi [-inf, 0]"
        else printfn "** Qagi [-inf, +inf]"
        printfn "result = %f, abserr = %g, info = %d" result abserr info

let TestDerkf() =
    let alfa = 3.141592653589793/4.0
    let alfasq = alfa*alfa
    let f (n: int, t: double, y: double[], yp: double[]) =
        let ysq = y.[0]*y.[0] + y.[1]*y.[1]
        let r = alfasq*(ysq**(-1.5))
        yp.[0] <- y.[2]
        yp.[1] <- y.[3]
        yp.[2] <- -y.[0]*r
        yp.[3] <- -y.[1]*r
    let n = 4
    let ecc = 0.25
    let lwork = 7*n + 20
    let work = Array.create lwork 0.0
    let liwork = 20
    let iwork = Array.create liwork 0
    let tfinal, tprint = 12.0, 1.0
    let y = [| 1.0 - ecc; 0.0; 0.0; alfa*(sqrt((1.0 + ecc)/(1.0 - ecc))) |]
    printfn "** Derkf"
    let mutable t = 0.0
    let mutable info = 0
    let mutable iter = true
    while iter do
        let tout = t + tprint
        let t1, info1 = XLPack.Derkf(info, n, f, t, y, tout, work, iwork)
        t <- t1
        info <- info1
        printfn "t = %f, y = %A" t y
        if t >= tfinal || info <> 1 then
            iter <- false
    printfn "info = %d" info

let TestDerkf_2() =
    let alfa = 3.141592653589793/4.0
    let alfasq = alfa*alfa
    let f (n: int, t: double, y: double[], yp: double[]) =
        let ysq = y.[0]*y.[0] + y.[1]*y.[1]
        let r = alfasq*(ysq**(-1.5))
        yp.[0] <- y.[2]
        yp.[1] <- y.[3]
        yp.[2] <- -y.[0]*r
        yp.[3] <- -y.[1]*r
    let n = 4
    let ecc = 0.25
    let lwork = 9*n + 20
    let work = Array.create lwork 0.0
    let liwork = 20
    let iwork = Array.create liwork 0
    let tfinal, tprint = 12.0, 1.0
    let y = [| 1.0 - ecc; 0.0; 0.0; alfa*(sqrt((1.0 + ecc)/(1.0 - ecc))) |]
    printfn "** Derkf (2) (dense output)"
    let mutable t = 0.0
    let mutable tout = t + tprint
    let mutable info = 0
    let mutable iter = true
    while iter do
        let t1, info1 = XLPack.Derkf(info, n, f, t, y, tfinal, work, iwork, mode = 2)
        t <- t1
        info <- info1
        if info = 1 || info = 2 then
            while t >= tout do
                let y1 = Array.create n 0.0
                XLPack.DerkfInt(n, tout, y1, work)
                printfn "t = %f, y = %A" tout y1
                tout <- tout + tprint
        else
            iter <- false
        if t >= tfinal then
            iter <- false
    printfn "info = %d" info

let TestRfft1() =
    let n = 10
    let r = Array.create n 0.0
    let rcopy = Array.create n 0.0
    let mutable info = 0
    // Initialization
    printfn "** Rfft1"
    let seed = 13ul
    XLPack.Init_Genrand(seed)
    let lwsave = n + (int ((log (double n))/(log 2.0))) + 4
    let wsave = Array.create lwsave 0.0
    info <- XLPack.Rfft1i(n, wsave)
    if info <> 0 then
        printfn "Error during initialization"
    // Generate test data
    if info = 0 then
        for i = 0 to n - 1 do
            r.[i] <- XLPack.Genrand_Res53()
            rcopy.[i] <- r.[i]
    // Forward transform
    if info = 0 then
        info <- XLPack.Rfft1f(n, r, wsave)
    if info <> 0 then
        printfn "Error in Rfft1f"
    // Backward transform
    if info = 0 then
        info <- XLPack.Rfft1b(n, r, wsave)
    if info <> 0 then
        printfn "Error in Rfft1b"
    // Check results
    if info = 0 then
        let mutable diff = 0.0
        for i = 0 to n - 1 do
            if abs(r.[i] - rcopy.[i]) > diff then
                diff <- abs(r.[i] - rcopy.[i])
            printfn "%g  %g  %g" rcopy.[i] r.[i] (abs(r.[i] - rcopy.[i]))
        printfn "diff(max) = %g" diff

let TestLmdif1() =
    let f(m: int, n: int, x: double[], y: double[], iflag: int) =
        let xdata = [| 77.6; 114.9; 141.1; 190.8; 239.9; 289.0; 332.8; 378.4; 434.8;
            477.3; 536.8; 593.1; 689.1; 760.0 |]
        let ydata = [| 10.07; 14.73; 17.94; 23.93; 29.61; 35.18; 40.02; 44.82; 50.76;
            55.05; 61.01;  66.4; 75.47; 81.78 |]
        for i = 0 to m - 1 do
            y.[i] <- ydata.[i] - x.[0]*(1.0 - exp(-xdata.[i]*x.[1]))
        0
    let m, n = 14, 2
    let x = [| 500.0; 0.0001 |]
    let lfvec = m * n
    let fvec = Array.create lfvec 0.0
    let info = XLPack.Lmdif1(f, m, n, x, fvec)
    printfn "** Lmdif1"
    printfn "x = %A" x
    printfn "info = %d" info

let TestRand() =
    let seed = 11ul
    printfn "** Random numbers : seed = %d" seed
    XLPack.Init_Genrand(seed)
    for i = 1 to 10 do
        printfn "%d  %d  %g" (XLPack.Genrand_Int32()) (XLPack.Genrand_Int31()) (XLPack.Genrand_Res53())

let TestDlamch() =
    printfn "** Dlamch"
    for c in [ 'e'; 's'; 'b'; 'p'; 'n'; 'r'; 'm'; 'u'; 'l'; 'o' ] do
        printfn "%c: %g" c (XLPack.Dlamch(c))

// Start Testing

TestSf()
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
TestRand()
TestDlamch()
