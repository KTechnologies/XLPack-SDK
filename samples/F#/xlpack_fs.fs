(*****************************************
 *                                       *
 *  Experimental F# interface to XLPack  *
 *  Version 6.1 (December 1, 2022)       *
 *  (C) 2014-2022  K Technologies        *
 *                                       *
 *****************************************)

namespace XLPack

open System.Runtime.InteropServices

//--- External function definitions

module External =

#if macOS
    [<Literal>]
    let Dll = "/Library/Application Support/Microsoft/Office365/User Content.localized/Add-Ins.localized/XLPack.dylib"
#else
#if Win64
    [<Literal>]
    let Dll = "XLPack.dll"
#else
    [<Literal>]
    let Dll = "XLPack_32.dll"
#endif
#endif

    type DBL2 = double[,]

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_d1num")>]
    extern double d1num(int i)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="x_factorial")>]
    extern double factorial(int ix, int& errno)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="x_li")>]
    extern double li(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="x_ei")>]
    extern double ei(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="x_e1")>]
    extern double e1(double x, int& errno)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="x_digamma")>]
    extern double digamma(double x, int& errno)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="x_besj0")>]
    extern double besj0(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="x_besj1")>]
    extern double besj1(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="x_besjnu")>]
    extern double besjnu(double nu, double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="x_besy0")>]
    extern double besy0(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="x_besy1")>]
    extern double besy1(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="x_besynu")>]
    extern double besynu(double nu, double x, int& errno)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besi0")>]
    extern double besi0(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besi0e")>]
    extern double besi0e(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besi1")>]
    extern double besi1(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besi1e")>]
    extern double besi1e(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besinu")>]
    extern double besinu(double nu, double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besk0")>]
    extern double besk0(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besk0e")>]
    extern double besk0e(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besk1")>]
    extern double besk1(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besk1e")>]
    extern double besk1e(double x, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besknu")>]
    extern double besknu(double nu, double x, int& errno)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_celli1")>]
    extern double celli1(double k, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_celli2")>]
    extern double celli2(double k, int& errno)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_celli3")>]
    extern double celli3(double n, double k, int& errno)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dconst")>]
    extern double dconst(int i)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dlange")>]
    extern double dlange(char norm, int m, int n, int lda, DBL2 a, double[] work)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dlansy")>]
    extern double dlansy(char norm, char uplo, int n, int lda, DBL2 a, double[] work)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dgesv")>]
    extern void dgesv(int n, int nrhs, int lda, DBL2 a, int[] ipiv, int ldb, double[] b, int& info)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dgecon")>]
    extern void dgecon(char norm, int n, int lda, DBL2 a, double anorm, double& rcond, double[] work, int[] iwork, int& info)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dposv")>]
    extern void dposv(char uplo, int n, int nrhs, int lda, DBL2 a, int ldb, double[] b, int& info)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dpocon")>]
    extern void dpocon(char uplo, int n, int lda, DBL2 a, double anorm, double& rcond, double[] work, int[] iwork, int& info)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dsyev")>]
    extern void dsyev(char jobz, char uplo, int n, int lda, DBL2 a, double[] w, double[] work, int lwork, int& info)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dgels")>]
    extern void dgels(char trans, int m, int n, int nrhs, int lda, DBL2 a, int ldb, double[] b, double[] work, int lwork, int& info)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dgecov")>]
    extern void dgecov(int job, int n, int lda, DBL2 a, double[] ci, int& info)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_pchse")>]
    extern void pchse(int n, double[] x, double[] f, double[] d, int incfd, double[] work, int lwork, int& info)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_pchfe")>]
    extern void pchfe(int n, double[] x, double[] f, double[] d, int incfd, int skip, int ne, double[] xe, double[] fe, int& info)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_pchia")>]
    extern double pchia(int n, double[] x, double[] f, double[] d, int incfd, int skip, double a, double b, int& info)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_rpzero2")>]
    extern void rpzero2(int n, double[] a, double[] rr, double[] ri, int iflag, int maxiter, int& iter, double[] s, double[] work, int& info)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dfzero_r")>]
    extern void dfzero_r(double& b, double& c, double r, double re, double ae, int& info, double& xx, double yy, int& irev)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_hybrd1_r")>]
    extern void hybrd1_r(int n, double[] x, double[] fvec, double xtol, double[] work, int lwork, int& info, double[] xx, double[] yy, int& irev)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dfmin_r")>]
    extern void dfmin_r(double a, double b, double tol, double& xx, double yy, int& irev)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_optif0_r")>]
    extern void optif0_r(int n, double[] x, double[] xpls, double& fpls, double[] work, int lwork, int& info, double[] xx, double yy, int& irev)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_qk15_r")>]
    extern void qk15_r(double a, double b, double& result, double& abserr, double& resabs, double& resasc, double& xx, double yy, int& irev)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_qag_r")>]
    extern void qag_r(double a, double b, double epsabs, double epsrel, int key, int limit, double& result, double& abserr, int& neval, int& last, double[] work, int lwork, int[] iwork, int& info, double& xx, double yy, int& irev)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_qagi_r")>]
    extern void qagi_r(double bound, int inf, double epsabs, double epsrel, int limit, double& result, double& abserr, int& neval, int& last, double[] work, int lwork, int[] iwork, int& info, double& xx, double yy, int& irev)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_derkf_r")>]
    extern void derkf_r(int n, double& t, double[] y, double tout, double& rtol, double& atol, int itol, int mode, double[] work, int lwork, int[] iwork, int liwork, int& info, double& tt, double[] yy, double[] yyp, int& irev)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_derkf_int")>]
    extern void derkf_int(int n, double t, double[] y, double[] work)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_rfft1f")>]
    extern void rfft1f(int n, int inc, double[] r, int lr, double[] wsave, int lwsave, double[] work, int lwork, int& info)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_rfft1b")>]
    extern void rfft1b(int n, int inc, double[] r, int lr, double[] wsave, int lwsave, double[] work, int lwork, int& info)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_rfft1i")>]
    extern void rfft1i(int n, double[] wsave, int lwsave, int& info)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_lmdif1_r")>]
    extern void lmdif1_r(int m, int n, double[] x, double[] fvec, double tol, double[] work, int lwork, int[] iwork, int& info, double[] xx, double[] yy, int& irev)

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_init_genrand")>]
    extern void init_genrand(uint32 s)
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_genrand_int32")>]
    extern uint32 genrand_int32()
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_genrand_int31")>]
    extern int genrand_int31()
    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_genrand_res53")>]
    extern double genrand_res53()

    [<DllImport(Dll, CallingConvention = CallingConvention.Cdecl, EntryPoint="_dlamch")>]
    extern double dlamch(char cmach)

//--- F# function definitions

type XLPack = class

    static member D1num(i: int) = External.d1num(i)

    static member Factorial(ix: int) =
        let mutable errno = 0
        let result = External.factorial(ix, &errno)
        (result, errno)

    static member Li(x: double) =
        let mutable errno = 0
        let result = External.li(x, &errno)
        (result, errno)

    static member Ei(x: double) =
        let mutable errno = 0
        let result = External.ei(x, &errno)
        (result, errno)

    static member E1(x: double) =
        let mutable errno = 0
        let result = External.e1(x, &errno)
        (result, errno)

    static member Digamma(x: double) =
        let mutable errno = 0
        let result = External.digamma(x, &errno)
        (result, errno)

    static member Besj0(x: double) =
        let mutable errno = 0
        let result = External.besj0(x, &errno)
        (result, errno)

    static member Besj1(x: double) =
        let mutable errno = 0
        let result = External.besj1(x, &errno)
        (result, errno)

    static member Besjnu(nu: double, x: double) =
        let mutable errno = 0
        let result = External.besjnu(nu, x, &errno)
        (result, errno)

    static member Besy0(x: double) =
        let mutable errno = 0
        let result = External.besy0(x, &errno)
        (result, errno)

    static member Besy1(x: double) =
        let mutable errno = 0
        let result = External.besy1(x, &errno)
        (result, errno)

    static member Besynu(nu: double, x: double) =
        let mutable errno = 0
        let result = External.besynu(nu, x, &errno)
        (result, errno)

    static member Besi0(x: double, ?kode: int) =
        let kode1 = defaultArg kode 1
        let mutable result, errno = 0.0, 0
        if kode1 = 2 then
            result <- External.besi0e(x, &errno)
        else
            result <- External.besi0(x, &errno)
        (result, errno)

    static member Besi1(x: double, ?kode: int) =
        let kode1 = defaultArg kode 1
        let mutable result, errno = 0.0, 0
        if kode1 = 2 then
            result <- External.besi1e(x, &errno)
        else
            result <- External.besi1(x, &errno)
        (result, errno)

    static member Besinu(nu: double, x: double) =
        let mutable errno = 0
        let result = External.besinu(nu, x, &errno)
        (result, errno)

    static member Besk0(x: double, ?kode: int) =
        let kode1 = defaultArg kode 1
        let mutable result, errno = 0.0, 0
        if kode1 = 2 then
            result <- External.besk0e(x, &errno)
        else
            result <- External.besk0(x, &errno)
        (result, errno)

    static member Besk1(x: double, ?kode: int) =
        let kode1 = defaultArg kode 1
        let mutable result, errno = 0.0, 0
        if kode1 = 2 then
            result <- External.besk1e(x, &errno)
        else
            result <- External.besk1(x, &errno)
        (result, errno)

    static member Besknu(nu: double, x: double) =
        let mutable errno = 0
        let result = External.besknu(nu, x, &errno)
        (result, errno)

    static member Celli1(k: double) =
        let mutable errno = 0
        let result = External.celli1(k, &errno)
        (result, errno)

    static member Celli2(k: double) =
        let mutable errno = 0
        let result = External.celli2(k, &errno)
        (result, errno)

    static member Celli3(n: double, k: double) =
        let mutable errno = 0
        let result = External.celli3(n, k, &errno)
        (result, errno)

    static member Dconst(i: int) = External.dconst(i)

    static member Dlange(norm: char, m: int, n: int, a: double[,]) =
        let lda = Array2D.length2 a
        let lwork = max m 1
        let work = Array.create lwork 0.0
        External.dlange(norm, m, n, lda, a, work)

    static member Dlansy(norm: char, uplo: char, n: int, a: double[,]) =
        let lda = Array2D.length2 a
        let lwork = max n 1
        let work = Array.create lwork 0.0
        External.dlansy(norm, uplo, n, lda, a, work)

    static member Dgesv(n: int, a: double[,], ipiv: int[], b: double[], ?nrhs: int) =
        let nrhs1 = defaultArg nrhs 1
        let lda = Array2D.length2 a
        let ldb = n
        let mutable info = 0
        External.dgesv(n, nrhs1, lda, a, ipiv, ldb, b, &info)
        info

    static member Dgecon(norm: char, n: int, a: double[,], anorm: double) =
        let lda = Array2D.length2 a
        let lwork = 4*(max n 1)
        let work = Array.create lwork 0.0
        let iwork = Array.create n 0
        let mutable rcond, info = 0.0, 0
        External.dgecon(norm, n, lda, a, anorm, &rcond, work, iwork, &info)
        (rcond, info)

    static member Dposv(uplo: char, n: int, a: double[,], b: double[], ?nrhs: int) =
        let nrhs1 = defaultArg nrhs 1
        let lda = Array2D.length2 a
        let ldb = n
        let mutable info = 0
        External.dposv(uplo, n, nrhs1, lda, a, ldb, b, &info)
        info

    static member Dpocon(uplo: char, n: int, a: double[,], anorm: double) =
        let lda = Array2D.length2 a
        let lwork = 3*(max n 1)
        let work = Array.create lwork 0.0
        let iwork = Array.create n 0
        let mutable rcond, info = 0.0, 0
        External.dpocon(uplo, n, lda, a, anorm, &rcond, work, iwork, &info)
        (rcond, info)

    static member Dsyev(jobz: char, uplo: char, n: int, a: double[,], w: double[]) =
        let lda = Array2D.length2 a
        let mutable lwork = -1
        let mutable work = Array.create 1 0.0
        let mutable info = 0
        External.dsyev(jobz, uplo, n, lda, a, w, work, lwork, &info)
        if info = 0 then
            lwork <- int work.[0]
            work <- Array.create lwork 0.0
            External.dsyev(jobz, uplo, n, lda, a, w, work, lwork, &info)
        info

    static member Dgels(trans: char, m: int, n: int, a: double[,], b: double[], ?nrhs: int) =
        let nrhs1 = defaultArg nrhs 1
        let lda = Array2D.length2 a
        let ldb = m
        let mutable lwork = -1
        let mutable work = Array.create 1 0.0
        let mutable info = 0
        External.dgels(trans, m, n, nrhs1, lda, a, ldb, b, work, lwork, &info)
        if info = 0 then
            lwork <- int work.[0]
            work <- Array.create lwork 0.0
            External.dgels(trans, m, n, nrhs1, lda, a, ldb, b, work, lwork, &info)
        info

    static member Dgecov(job: int, n: int, a: double[,], ci: double[]) =
        let lda = Array2D.length2 a
        let mutable info = 0
        External.dgecov(job, n, lda, a, ci, &info)
        info

    static member Pchse(n: int, x: double[], f: double[], d: double[], ?incfd: int) =
        let incfd1 = defaultArg incfd 1
        let lwork = 2*(max n 1)
        let work = Array.create lwork 0.0
        let mutable info = 0
        External.pchse(n, x, f, d, incfd1, work, lwork, &info)
        info

    static member Pchfe(n: int, x: double[], f: double[], d: double[], ne: int, xe: double[], fe: double[], ?incfd: int, ?skip: int) =
        let incfd1 = defaultArg incfd 1
        let skip1 = defaultArg skip 0
        let mutable info = 0
        External.pchfe(n, x, f, d, incfd1, skip1, ne, xe, fe, &info)
        info

    static member Pchia(n: int, x: double[], f: double[], d: double[], a: double, b: double, ?incfd: int, ?skip: int) =
        let incfd1 = defaultArg incfd 1
        let skip1 = defaultArg skip 0
        let mutable info = 0
        let mutable s = External.pchia(n, x, f, d, incfd1, skip1, a, b, &info)
        (s, info)

    static member Rpzero2(n: int, a: double[], rr: double[], ri: double[], s: double[], ?iflag: int, ?maxiter: int) =
        let iflag1 = defaultArg iflag 0
        let maxiter1 = defaultArg maxiter 25*n
        let lwork = 8*(max n 1) + 6
        let work = Array.create lwork 0.0
        let mutable iter = 0
        let mutable info = 0
        External.rpzero2(n, a, rr, ri, iflag1, maxiter1, &iter, s, work, &info)
        (iter, info)

    static member Dfzero(f: double -> double, b: byref<double>, c: byref<double>, r: double, ?re: double, ?ae: double) =
        let re1 = defaultArg re 1.0e-10
        let ae1 = defaultArg ae 1.0e-10
        let mutable info = 0
        let mutable irev, xx, yy = 0, 0.0, 0.0
        let mutable iter = true
        while iter do
            External.dfzero_r(&b, &c, r, re1, ae1, &info, &xx, yy, &irev)
            if irev <> 0 then
                yy <- f(xx)
            else
                iter <- false
        info

    static member Hybrd1(f: int * double[] * double[] * int -> int, n: int, x: double[], fvec: double[], ?xtol: double) =
        let xtol1 = defaultArg xtol 1.0e-10
        if n < 1 then
            -2
        elif xtol1 < 0.0 then
            -5
        else
            let lwork = n*(3*n + 13)/2
            let work = Array.create lwork 0.0
            let xx = Array.create n 0.0
            let yy = Array.create n 0.0
            let mutable info, irev = 0, 0
            let mutable iter = true
            while iter do
                External.hybrd1_r(n, x, fvec, xtol1, work, lwork, &info, xx, yy, &irev)
                if irev = 1 || irev = 2 then
                    if f(n, xx, yy, 1) < 0 then
                        info <- 5
                        irev <- 0
                elif irev = 3 || irev = 4 then
                    if f(n, xx, yy, 2) < 0 then
                        info <- 5
                        irev <- 0
                if irev = 0 then
                    iter <- false
            info

    static member Dfmin(a: double, b: double, f: double -> double, ?tol: double) =
        let tol1 = defaultArg tol 1.0e-10
        let mutable irev, xx, yy = 0, 0.0, 0.0
        let mutable iter = true
        while iter do
            External.dfmin_r(a, b, tol1, &xx, yy, &irev)
            if irev <> 0 then
                yy <- f(xx)
            else
                iter <- false
        xx

    static member Optif0(n: int, x: double[], f: int * double[] -> double, xpls: double[]) =
        if n < 1 then
            (0.0, -1)
        else
            let lwork = n*(n + 10)
            let work = Array.create lwork 0.0
            let mutable fpls, info, irev = 0.0, 0, 0
            let xx = Array.create n 0.0
            let mutable yy = 0.0
            let mutable iter = true
            while iter do
                External.optif0_r(n, x, xpls, &fpls, work, lwork, &info, xx, yy, &irev)
                if irev >= 1 && irev <= 20 then
                    yy <- f(n, xx)
                if irev = 0 then
                    iter <- false
            (fpls, info)

    static member Qk15(f: double -> double, a: double, b: double) =
        let mutable result, abserr, resabs, resasc = 0.0, 0.0, 0.0, 0.0
        let mutable irev, xx, yy = 0, 0.0, 0.0
        let mutable iter = true
        while iter do
            External.qk15_r(a, b, &result, &abserr, &resabs, &resasc, &xx, yy, &irev)
            if irev >= 1 && irev <= 5 then
                yy <- f(xx)
            else
                iter <- false
        (result, abserr)

    static member Qag(f: double -> double, a: double, b: double, ?epsabs: double, ?epsrel: double, ?key: int, ?limit: int) =
        let epsabs1 = defaultArg epsabs 1.0e-10
        let epsrel1 = defaultArg epsrel 1.0e-10
        let key1 = defaultArg key 1
        let limit1 = defaultArg limit 100
        if limit1 < 1 then
            (0.0, 0.0, -7)
        else
            let lwork = 4*limit1
            let work = Array.create lwork 0.0
            let liwork = limit1
            let iwork = Array.create liwork 0
            let mutable result, abserr, neval, last = 0.0, 0.0, 0, 0
            let mutable info, irev, xx, yy = 0, 0, 0.0, 0.0
            let mutable iter = true
            while iter do
                External.qag_r(a, b, epsabs1, epsrel1, key1, limit1, &result, &abserr, &neval, &last, work, lwork, iwork, &info, &xx, yy, &irev)
                if irev >= 1 && irev <= 15 then
                    yy <- f(xx)
                else
                    iter <- false
            (result, abserr, info)

    static member Qagi(f: double -> double, bound: double, inf: int, ?epsabs: double, ?epsrel: double, ?limit: int) =
        let epsabs1 = defaultArg epsabs 1.0e-10
        let epsrel1 = defaultArg epsrel 1.0e-10
        let limit1 = defaultArg limit 100
        if limit1 < 1 then
            (0.0, 0.0, -6)
        else
            let lwork = 4*limit1
            let work = Array.create lwork 0.0
            let liwork = limit1
            let iwork = Array.create liwork 0
            let mutable result, abserr, neval, last = 0.0, 0.0, 0, 0
            let mutable info, irev, xx, yy = 0, 0, 0.0, 0.0
            let mutable iter = true
            while iter do
                External.qagi_r(bound, inf, epsabs1, epsrel1, limit1, &result, &abserr, &neval, &last, work, lwork, iwork, &info, &xx, yy, &irev)
                if irev >= 1 && irev <= 18 then
                    yy <- f(xx)
                else
                    iter <- false
            (result, abserr, info)

    static member Derkf(info: int, n: int, f: int * double * double[] * double[] -> unit, t: double, y: double[], tout: double, work: double[], iwork: int[], ?rtol: double, ?atol: double, ?mode: int) =
        let mutable rtol1 = defaultArg rtol 1.0e-10
        let mutable atol1 = defaultArg atol 0.0
        let mode1 = defaultArg mode 0
        if n < 1 then
            (t, -1)
        else
            let itol = 0
            let lwork = Array.length work
            let liwork = Array.length iwork
            let mutable tt = 0.0
            let yy = Array.create n 0.0
            let yyp = Array.create n 0.0
            let mutable info1, t1, irev = info, t, 0
            let mutable iter = true
            while iter do
                External.derkf_r(n, &t1, y, tout, &rtol1, &atol1, itol, mode1, work, lwork, iwork, liwork, &info1, &tt, yy, yyp, &irev)
                if irev >= 1 && irev <= 11 then
                    f(n, tt, yy, yyp)
                else
                    iter <- false
            if info1 = -13 then
                info1 <- -1
            if info1 = -10 then
                info1 <- -7
            elif info1 = -12 then
                info1 <- -8
            elif info1 = -5 then
                info1 <- -9
            elif info1 = -6 then
                info1 <- -10
            elif info1 < 0 then
                info1 <- info1 - 2
            (t1, info1)

    static member DerkfInt(n: int, t: double, y: double[], work: double[]) =
        if not (n < 1 || Array.length y < n || Array.length work < 9*n + 20) then
            External.derkf_int(n, t, y, work)

    static member Rfft1f(n: int, r: double[], wsave: double[], ?inc: int) =
        let inc1 = defaultArg inc 1
        if n < 1 then
            -1
        else
            let lr = Array.length r
            let lwsave = Array.length wsave
            let lwork = n
            let work = Array.create lwork 0.0
            let mutable info = 0
            External.rfft1f(n, inc1, r, lr, wsave, lwsave, work, lwork, &info)
            info

    static member Rfft1b(n: int, r: double[], wsave: double[], ?inc: int) =
        let inc1 = defaultArg inc 1
        if n < 1 then
            -1
        else
            let lr = Array.length r
            let lwsave = Array.length wsave
            let lwork = n
            let work = Array.create lwork 0.0
            let mutable info = 0
            External.rfft1b(n, inc1, r, lr, wsave, lwsave, work, lwork, &info)
            info

    static member Rfft1i(n: int, wsave: double[]) =
        let lwsave = Array.length wsave
        let mutable info = 0
        External.rfft1i(n, wsave, lwsave, &info)
        info

    static member Lmdif1(f: int * int * double[] * double[] * int -> int, m: int, n: int, x: double[], fvec: double[], ?tol: double) =
        let tol1 = defaultArg tol 1.0e-10
        if m < n then
            -2
        elif n < 1 then
            -3
        elif tol1 < 0.0 then
            -6
        else
            let lwork = n*(m + 6) + 2*m
            let work = Array.create lwork 0.0
            let liwork = n
            let iwork = Array.create liwork 0
            let xx = Array.create n 0.0
            let yy = Array.create m 0.0
            let mutable info, irev = 0, 0
            let mutable iter = true
            while iter do
                External.lmdif1_r(m, n, x, fvec, tol1, work, lwork, iwork, &info, xx, yy, &irev)
                if irev = 1 || irev = 2 then
                    if f(m, n, xx, yy, 1) < 0 then
                        info <- 5
                        irev <- 0
                elif irev = 3 then
                    if f(m, n, xx, yy, 2) < 0 then
                        info <- 5
                        irev <- 0
                if irev = 0 then
                    iter <- false
            info

    static member Init_Genrand(s: uint32) = External.init_genrand(s)
    static member Genrand_Int32() = External.genrand_int32()
    static member Genrand_Int31() = External.genrand_int31()
    static member Genrand_Res53() = External.genrand_res53()

    static member Dlamch(cmach: char) = External.dlamch(cmach)

    end
