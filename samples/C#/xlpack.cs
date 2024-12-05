/*****************************************
 *                                       *
 *  Experimental C# interface to XLPack  *
 *  Version 7.0 (January 31, 2023)       *
 *  (C) 2014-2023  K Technologies        *
 *                                       *
 *****************************************/

using static System.Math;
using System.Runtime.InteropServices;

using System.Text;

public static class XLPack
{

#if Win32
	const string DLL = "XLPack_32.dll";
	const string DLLSP = "XLPackSp_32.dll";
	const string DLLPDE = "XLPackPde_32.dll";
#else
	const string DLL = "XLPack.dll";
	const string DLLSP = "XLPackSp.dll";
	const string DLLPDE = "XLPackPde.dll";
#endif

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_d1num")]
public extern static double d1num(int i);

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_factorial")]
public extern static double factorial(uint x, out int errno);

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_li")]
public extern static double li(double x, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_ei")]
public extern static double ei(double x, out int errno);

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_digamma")]
public extern static double digamma(double x, out int errno);

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besj0")]
public extern static double besj0(double x, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besj1")]
public extern static double besj1(double x, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besjnu")]
public extern static double besjnu(double nu, double x, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besy0")]
public extern static double besy0(double x, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besy1")]
public extern static double besy1(double x, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besynu")]
public extern static double besynu(double nu, double x, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besi0")]
public extern static double _besi0(double x, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besi0e")]
public extern static double _besi0e(double x, out int errno);

public static double besi0(double x, out int errno, int kode = 1)
{
	if (kode == 2)
		return _besi0e(x, out errno);
	else
		return _besi0(x, out errno);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besi1")]
public extern static double _besi1(double x, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besi1e")]
public extern static double _besi1e(double x, out int errno);

public static double besi1(double x, out int errno, int kode = 1)
{
	if (kode == 2)
		return _besi1e(x, out errno);
	else
		return _besi1(x, out errno);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besinu")]
public extern static double besinu(double nu, double x, out int errno);

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besk0")]
public extern static double _besk0(double x, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besk0e")]
public extern static double _besk0e(double x, out int errno);

public static double besk0(double x, out int errno, int kode = 1)
{
	if (kode == 2)
		return _besk0e(x, out errno);
	else
		return _besk0(x, out errno);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besk1")]
public extern static double _besk1(double x, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besk1e")]
public extern static double _besk1e(double x, out int errno);

public static double besk1(double x, out int errno, int kode = 1)
{
	if (kode == 2)
		return _besk1e(x, out errno);
	else
		return _besk1(x, out errno);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_besknu")]
public extern static double besknu(double nu, double x, out int errno);


[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_celli1")]
public extern static double celli1(double k, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_celli2")]
public extern static double celli2(double k, out int errno);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "x_celli3")]
public extern static double celli3(double n, double k, out int errno);

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dconst")]
public extern static double dconst(int i);

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dlange")]
private extern static double _dlange(char norm, int m, int n, int lda, double[,] a, double[] work);

public static double dlange(char norm, int m, int n, double[,] a)
{
	int lda = a.GetLength(1);
	double[] work = new double[Max(1, m)];
	return _dlange(norm, m, n, lda, a, work);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dlansy")]
private extern static double _dlansy(char norm, char uplo, int n, int lda, double[,] a, double[] work);

public static double dlansy(char norm, char uplo, int n, double[,] a)
{
	int lda = a.GetLength(1);
	double[] work = new double[Max(1, n)];
	return _dlansy(norm, uplo, n, lda, a, work);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dgesv")]
public extern static void _dgesv(int n, int nrhs, int lda, double[,] a, int[] ipiv, int ldb, double[] b, out int info);

public static void dgesv(int n, double[,] a, int[] ipiv, double[] b, out int info, int nrhs = 1)
{
	int lda = a.GetLength(1);
	int ldb = n;
	_dgesv(n, nrhs, lda, a, ipiv, ldb, b, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dgecon")]
private extern static void _dgecon(char norm, int n, int lda, double[,] a, double anorm, out double rcond, double[] work, int[] iwork, out int info);

public static void dgecon(char norm, int n, double[,] a, double anorm, out double rcond, out int info)
{
	int lda = a.GetLength(1);
	double[] work = new double[Max(1, 4*n)];
	int[] iwork = new int[Max(1, n)];
	_dgecon(norm, n, lda, a, anorm, out rcond, work, iwork, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dposv")]
public extern static void _dposv(char uplo, int n, int nrhs, int lda, double[,] a, int ldb, double[] b, out int info);

public static void dposv(char uplo, int n, double[,] a, double[] b, out int info, int nrhs = 1)
{
	int lda = a.GetLength(1);
	int ldb = n;
	_dposv(uplo, n, nrhs, lda, a, ldb, b, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dpocon")]
private extern static void _dpocon(char uplo, int n, int lda, double[,] a, double anorm, out double rcond, double[] work, int[] iwork, out int info);

public static void dpocon(char uplo, int n, double[,] a, double anorm, out double rcond, out int info)
{
	int lda = a.GetLength(1);
	double[] work = new double[Max(1, 3*n)];
	int[] iwork = new int[Max(1, n)];
	_dpocon(uplo, n, lda, a, anorm, out rcond, work, iwork, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dsyev")]
private extern static void _dsyev(char jobz, char uplo, int n, int lda, double[,] a, double[] w, double[] work, int lwork, out int info);

public static void dsyev(char jobz, char uplo, int n, double[,] a, double[] w, out int info)
{
	int lda = a.GetLength(1);
	int lwork = -1;
	double[] work = new double[1];
	_dsyev(jobz, uplo, n, lda, a, w, work, lwork, out info);
	if (info == 0) {
		lwork = (int)work[0];
		work = new double[lwork];
		_dsyev(jobz, uplo, n, lda, a, w, work, lwork, out info);
	}
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dgels")]
private extern static void _dgels(char trans, int m, int n, int nrhs, int lda, double[,] a, int ldb, double[] b, double[] work, int lwork, out int info);

public static void dgels(char trans, int m, int n, double[,] a, double[] b, out int info, int nrhs = 1)
{
	int lda = a.GetLength(1);
	int ldb = m;
	int lwork = -1;
	double[] work = new double[1];
	_dgels(trans, m, n, nrhs, lda, a, ldb, b, work, lwork, out info);
	if (info == 0) {
		lwork = (int)work[0];
		work = new double[lwork];
		_dgels(trans, m, n, nrhs, lda, a, ldb, b, work, lwork, out info);
	}
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dgecov")]
public extern static void _dgecov(int job, int n, int lda, double[,] a, double[] ci, out int info);

public static void dgecov(int job, int n, double[,] a, double[] ci, out int info)
{
	int lda = a.GetLength(1);
	_dgecov(job, n, lda, a, ci, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_pchse")]
private extern static void _pchse(int n, double[] x, double[] f, double[] d, int incfd, double[] work, int lwork, out int info);

public static void pchse(int n, double[] x, double[] f, double[] d, out int info, int incfd = 1)
{
	int lwork = Max(1, 2*n);
	double[] work = new double[lwork];
	_pchse(n, x, f, d, incfd, work, lwork, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_pchfe")]
public extern static void _pchfe(int n, double[] x, double[] f, double[] d, int incfd, int skip, int ne, double[] xe, double[] fe, out int info);

public static void pchfe(int n, double[] x, double[] f, double[] d, int ne, double[] xe, double[] fe, out int info, int incfd = 1, int skip = 0)
{
	_pchfe(n, x, f, d, incfd, skip, ne, xe, fe, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_pchia")]
public extern static double _pchia(int n, double[] x, double[] f, double[] d, int incfd, int skip, double a, double b, out int info);

public static double pchia(int n, double[] x, double[] f, double[] d, double a, double b, out int info, int incfd = 1, int skip = 0)
{
	return _pchia(n, x, f, d, incfd, skip, a, b, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_rpzero2")]
private extern static void _rpzero2(int n, double[] a, double[] rr, double[] ri, int iflag, int maxiter, out int iter, double[] s, double[] work, out int info);

public static void rpzero2(int n, double[] a, double[] rr, double[] ri, out int iter, double[] s, out int info, int iflag = 0, int maxiter = -1)
{
	if (maxiter <= 0)
		maxiter = 25*n;
	int lwork = Max(1, 8*n + 6);
	double[] work = new double[lwork];
	_rpzero2(n, a, rr, ri, iflag, maxiter, out iter, s, work, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dfzero_r")]
private extern static void _dfzero_r(out double b, out double c, double r, double re, double ae, out int info, out double xx, double yy, out int irev);

public delegate double DfzeroFunc(double x);

public static void dfzero(DfzeroFunc f, out double b, out double c, double r, out int info, double re = 1.0e-10, double ae = 1.0e-10)
{
	int irev;
	double xx, yy;

	irev = 0;
	yy = 0;
	do {
		_dfzero_r(out b, out c, r, re, ae, out info, out xx, yy, out irev);
		if (irev != 0)
			yy = f(xx);
	} while (irev != 0);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_hybrd1_r")]
private extern static void _hybrd1_r(int n, double[] x, double[] fvec, double tol, double[] work, int lwork, out int info, double[] xx, double[] yy, out int irev);

public delegate void Hybrd1Proc(int n, double[] x, double[] fvec, ref int iflag);

public static void hybrd1(Hybrd1Proc fcn, int n, double[] x, double[] fvec, out int info, double tol = 1.0e-10)
{
	int iflag, irev;
	double[] xx, yy;
	int lwork = n*(3*n + 13)/2;
	double[] work;

	if (n < 1) {
		info = -2;
		return;
	}
	if (tol < 0) {
		info = -5;
		return;
	}
	work = new double[lwork];
	xx = new double[n];
	yy = new double[n];
	irev = 0;
	do {
		_hybrd1_r(n, x, fvec, tol, work, lwork, out info, xx, yy, out irev);
		if (irev == 1 || irev == 2) {
			iflag = 1;
			fcn(n, xx, yy, ref iflag);
		} else if (irev == 3 || irev == 4) {
			iflag = 2;
			fcn(n, xx, yy, ref iflag);
		}
	} while (irev != 0);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dfmin_r")]
private extern static void _dfmin_r(double a, double b, double tol, out double xx, double yy, out int irev);

public delegate double DfminFunc(double x);

public static double dfmin(double a, double b, DfminFunc f, double tol = 1.0e-10)
{
	int irev;
	double xx, yy;

	irev = 0;
	yy = 0;
	do {
		_dfmin_r(a, b, tol, out xx, yy, out irev);
		if (irev != 0)
			yy = f(xx);
	} while (irev != 0);
	return xx;
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_optif0_r")]
private extern static void _optif0_r(int n, double[] x, double[] xpls, out double fpls, double[] work, int lwork, out int info, double[] xx, double yy, out int irev);

public delegate void Optif0Proc(int n, double[] x, out double fval);

public static void optif0(int n, double[] x, Optif0Proc fcn, double[] xpls, out double fpls, out int info)
{
	double[] xx;
	double yy;
	int irev;
	int lwork = n*(n + 10);
	double[] work;

	if (n < 1) {
		info = -1;
		fpls = 0;
		return;
	}
	work = new double[lwork];
	xx = new double[n];
	irev = 0;
	yy = 0;
	do {
		_optif0_r(n, x, xpls, out fpls, work, lwork, out info, xx, yy, out irev);
		if (irev >= 1 && irev <= 20)
			fcn(n, xx, out yy);
	} while (irev != 0);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_qk15_r")]
private extern static void _qk15_r(double a, double b, out double result, out double abserr, out double resabs, out double resasc, out double xx, double yy, out int irev);

public delegate double Qk15Func(double x);

public static void qk15(Qk15Func f, double a, double b, out double result, out double abserr)
{
	double resabs, resasc, xx, yy;
	int irev;

	irev = 0;
	yy = 0;
	do {
		_qk15_r(a, b, out result, out abserr, out resabs, out resasc, out xx, yy, out irev);
		if (irev >= 1 && irev <= 5)
			yy = f(xx);
	} while (irev != 0);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_qag_r")]
private extern static void _qag_r(double a, double b, double epsabs, double epsrel, int key, int limit, out double result, out double abserr, out int neval, out int last, double[] work, int lwork, int[] iwork, out int info, out double xx, double yy, out int irev);

public delegate double QagFunc(double x);

public static void qag(QagFunc f, double a, double b, out double result, out double abserr, out int neval, out int last, out int info, double epsabs = 1.0e-10, double epsrel = 1.0e-10, int key = 1, int limit = 100)
{
	double xx, yy;
	int irev;
	int lwork = 4*limit, liwork = limit;
	double[] work;
	int[] iwork;

	if (limit < 1) {
		info = -7;
		result = 0;
		abserr = 0;
		neval = 0;
		last = 0;
		return;
	}
	work = new double[lwork];
	iwork = new int[liwork];
	irev = 0;
	yy = 0;
	do {
		_qag_r(a, b, epsabs, epsrel, key, limit, out result, out abserr, out neval, out last, work, lwork, iwork, out info, out xx, yy, out irev);
		if (irev >= 1 && irev <= 15)
			yy = f(xx);
	} while (irev != 0);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_qagi_r")]
private extern static void _qagi_r(double bound, int inf, double epsabs, double epsrel, int limit, out double result, out double abserr, out int neval, out int last, double[] work, int lwork, int[] iwork, out int info, out double xx, double yy, out int irev);

public delegate double QagiFunc(double x);

public static void qagi(QagiFunc f, double bound, int inf, out double result, out double abserr, out int neval, out int last, out int info, double epsabs = 1.0e-10, double epsrel = 1.0e-10, int limit = 100)
{
	double xx, yy;
	int irev;
	int lwork = 4*limit, liwork = limit;
	double[] work;
	int[] iwork;

	if (limit < 1) {
		info = -6;
		result = 0;
		abserr = 0;
		neval = 0;
		last = 0;
		return;
	}
	work = new double[lwork];
	iwork = new int[liwork];
	irev = 0;
	yy = 0;
	do {
		_qagi_r(bound, inf, epsabs, epsrel, limit, out result, out abserr, out neval, out last, work, lwork, iwork, out info, out xx, yy, out irev);
		if (irev >= 1 && irev <= 18)
			yy = f(xx);
	} while (irev != 0);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_derkfa_r")]
private extern static void _derkfa_r(int n, ref double t, double[] y, double tout, double tend, ref double rtol, ref double atol, int itol, int mode, double[] work, int lwork, int[] iwork, int liwork, out int info, out double tt, double[] yy, double[] yyp, out int irev);

public delegate void DerkfaProc(int n, double t, double[] y, double[] yp);

public static void derkfa(int n, DerkfaProc fcn, ref double t, double[] y, double tout, double tend, int mode, double[] work, int[] iwork, out int info, double rtol = 1.0e-10, double atol = 0.0)
{
	const int itol = 0;
	int lwork = work.Length;
	int liwork = iwork.Length;
	double tt;
	double[] yy, yyp;
	int irev;

	if (n < 1) {
		info = -1;
		return;
	}
	yy = new double[n];
	yyp = new double[n];
	irev = 0;
	do {
		_derkfa_r(n, ref t, y, tout, tend, ref rtol, ref atol, itol, mode, work, lwork, iwork, liwork, out info, out tt, yy, yyp, out irev);
		if (irev != 0)
			fcn(n, tt, yy, yyp);
	} while (irev != 0);
	if (info == -2 || info == -5)
		info = info - 1;
	else if (info == -6 || info == -7)
		info = info - 4;
	else if (info == -9)
		info = -12;
	else if (info == -11)
		info = -7;
	else if (info == -13)
		info = -8;
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dopn43_r")]
private extern static void _dopn43_r(int n, ref double t, double[] y, double[] yp, double tout, double tend, ref double rtol, ref double atol, int itol, int mode, double[] work, int lwork, int[] iwork, int liwork, out int info, out double tt, double[] yy, double[] yypp, out int irev);

public delegate void Dopn43Proc(int n, double t, double[] y, double[] ypp);

public static void dopn43(int n, Dopn43Proc fcn, ref double t, double[] y, double[] yp, double tout, double tend, int mode, double[] work, int[] iwork, out int info, double rtol = 1.0e-10, double atol = 0.0)
{
	const int itol = 0;
	int lwork = work.Length;
	int liwork = iwork.Length;
	double tt;
	double[] yy, yypp;
	int irev;

	if (n < 1) {
		info = -1;
		return;
	}
	yy = new double[n];
	yypp = new double[n];
	irev = 0;
	do {
		_dopn43_r(n, ref t, y, yp, tout, tend, ref rtol, ref atol, itol, mode, work, lwork, iwork, liwork, out info, out tt, yy, yypp, out irev);
		if (irev >= 1 && irev <= 2)
			fcn(n, tt, yy, yypp);
	} while (irev != 0);
	if (info == -2 || info == -6)
		info = info - 1;
	else if (info == -7 || info == -8)
		info = info - 4;
	else if (info == -10)
		info = -13;
	else if (info == -12)
		info = -8;
	else if (info == -14)
		info = -9;
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_rfft1f")]
private extern static void _rfft1f(int n, int inc, double[] r, int lr, double[] wsave, int lwsave, double[] work, int lwork, out int info);

public static void rfft1f(int n, double[] r, double[] wsave, out int info, int inc = 1)
{
	int lr = r.Length;
	int lwsave = wsave.Length;
	int lwork = n;
	double[] work;

	if (n < 1) {
		info = -1;
		return;
	}
	work = new double[lwork];
	_rfft1f(n, inc, r, lr, wsave, lwsave, work, lwork, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_rfft1b")]
private extern static void _rfft1b(int n, int inc, double[] r, int lr, double[] wsave, int lwsave, double[] work, int lwork, out int info);

public static void rfft1b(int n, double[] r, double[] wsave, out int info, int inc = 1)
{
	int lr = r.Length;
	int lwsave = wsave.Length;
	int lwork = n;
	double[] work;

	if (n < 1) {
		info = -1;
		return;
	}
	work = new double[lwork];
	_rfft1b(n, inc, r, lr, wsave, lwsave, work, lwork, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_rfft1i")]
public extern static void _rfft1i(int n, double[] wsave, int lwsave, out int info);

public static void rfft1i(int n, double[] wsave, out int info)
{
	int lwsave = wsave.Length;

	if (n < 1) {
		info = -1;
		return;
	}
	_rfft1i(n, wsave, lwsave, out info);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_lmdif1_r")]
private extern static void _lmdif1_r(int m, int n, double[] x, double[] fvec, double tol, double[] work, int lwork, int[] iwork, out int info, double[] xx, double[] yy, out int irev);

public delegate void Lmdif1Proc(int m, int n, double[] x, double[] fvec, ref int iflag);

public static void lmdif1(Lmdif1Proc fcn, int m, int n, double[] x, double[] fvec, out int info, double tol = 1.0e-10)
{
	int iflag, irev;
	double[] xx, yy;
	int lwork = n*(m + 6) + 2*m, liwork = n;
	double[] work;
	int[] iwork;

	if (m < n) {
		info = -2;
		return;
	}
	if (n < 1) {
		info = -3;
		return;
	}
	if (tol < 0) {
		info = -6;
		return;
	}
	work = new double[lwork];
	iwork = new int[liwork];
	xx = new double[n];
	yy = new double[m];
	irev = 0;
	do {
		_lmdif1_r(m, n, x, fvec, tol, work, lwork, iwork, out info, xx, yy, out irev);
		if (irev == 1 || irev == 2) {
			iflag = 1;
			fcn(m, n, xx, yy, ref iflag);
		} else if (irev == 3) {
			iflag = 2;
			fcn(m, n, xx, yy, ref iflag);
		}
	} while (irev != 0);
}

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_init_genrand")]
public extern static void init_genrand(uint s);
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_genrand_int32")]
public extern static uint genrand_int32();
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_genrand_int31")]
public extern static int genrand_int31();
[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_genrand_res53")]
public extern static double genrand_res53();

[DllImport(DLL, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dlamch")]
public extern static double dlamch(char cmach);

//---

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_csr_dusmv")]
public extern static void csr_dusmv(char trans, int m, int n, double alpha, double[] val, int[] rowptr, int[] colind, int ibase, double[] x, int incx, double beta, double[] y, int incy, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_csr_dussv")]
public extern static void csr_dussv(char uplo, char trans, char diag, int n, double[] val, int[] rowptr, int[] colind, int ibase, double[] x, int incx, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_csr_dusmm")]
public extern static void csr_dusmm(char trans, char order, int m, int n, int nrhs, double alpha, double[] val, int[] rowptr, int[] colind, int ibase, int ldb, double[] b, double beta, int ldc, double[] c, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_csr_dussm")]
public extern static void csr_dussm(char uplo, char trans, char diag, char order, int n, int nrhs, double[] val, int[] rowptr, int[] colind, int ibase, int idx, double[] x, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_ssr_dusmv")]
public extern static void ssr_dusmv(char uplo, int n, double alpha, double[] val, int[] rowptr, int[] colind, int ibase, double[] x, int incx, double beta, double[] y, int incy, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_csc_csr")]
public extern static void csc_csr(int m, int n, double[] val, int[] ptr, int[] ind, int ibase, double[] val2, int[] ptr2, int[] ind2, int ibase2, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_coo_csr")]
public extern static void coo_csr(int m, int n, int nnz, double[] val, int[] rowind, int[] colind, int ibase, double[] val2, int[] rowptr2, int[] colind2, int ibase2, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_csr_coo")]
public extern static void csr_coo(int m, int n, double[] val, int[] rowptr, int[] colind, int ibase, double[] val2, int[] rowind2, int[] colind2, int ibase2, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_csr_ssr")]
public extern static void csr_ssr(char uplo, int n, double[] val, int[] rowptr, int[] colind, int ibase, int maxnnz2, double[] val2, int[] rowptr2, int[] colind2, int ibase2, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_ssr_csr")]
public extern static void ssr_csr(char uplo, int n, double[] val, int[] rowptr, int[] colind, int ibase, int maxnnz2, double[] val2, int[] rowptr2, int[] colind2, int ibase2, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_dense_csr")]
public extern static void dense_csr(int m, int n, int lda, double[] a, int maxnnz, double[] val, int[] rowptr, int[] colind, int ibase, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_csr_dense")]
public extern static void csr_dense(int m, int n, double[] val, int[] rowptr, int[] colind, int ibase, int lda, double[] a, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_bicg1")]
public extern static void _bicg1(int n, double[] val, int[] rowptr, int[] colind, double[] b, double[] x, double tol, int maxiter, out int iter, out double res, int lwork, double[] work, out int info);

public static void bicg1(int n, double[] val, int[] rowptr, int[] colind, double[] b, double[] x, double tol, int maxiter, out int iter, out double res, out int info)
{
	int lwork = 8*n;
	double[] work;
	work = new double[lwork];
	_bicg1(n, val, rowptr, colind, b, x, tol, maxiter, out iter, out res, lwork, work, out info);
}

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_cg1")]
public extern static void _cg1(char uplo, int n, double[] val, int[] rowptr, int[] colind, double[] b, double[] x, double tol, int maxiter, out int iter, out double res, int lwork, double[] work, out int info);

public static void cg1(char uplo, int n, double[] val, int[] rowptr, int[] colind, double[] b, double[] x, double tol, int maxiter, out int iter, out double res, out int info)
{
	int lwork = 5*n;
	double[] work;
	work = new double[lwork];
	_cg1(uplo, n, val, rowptr, colind, b, x, tol, maxiter, out iter, out res, lwork, work, out info);
}

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_mm_read")]
public extern static void _mm_read(string fname, byte[] matcode, out int nrow, out int ncol, out int nnz, double[] val, int lval, int[] rowind, int lrowind, int[] colind, int lcolind, int skip, int ibase, int format, int sort, out int info);

public static void mm_read(string fname, out string matcode, out int nrow, out int ncol, out int nnz, double[] val, int lval, int[] ptr, int lptr, int[] ind, int lind, int skip, int ibase, int format, int sort, out int info)
{
	var encoding = Encoding.GetEncoding("UTF-8");
	byte[] bmatcode = new byte[4];
	_mm_read(fname, bmatcode, out nrow, out ncol, out nnz, val, lval, ptr, lptr, ind, lind, skip, ibase, format, sort, out info);
	matcode = encoding.GetString(bmatcode);
}

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_mm_write")]
public extern static void _mm_write(string fname, char[] matcode, int nrow, int ncol, int nnz, double[] val, int lval, int[] ptr, int lptr, int[] ind, int lind, int ibase, int fchk, int format, out int info);

public static void mm_write(string fname, string matcode, int nrow, int ncol, int nnz, double[] val, int lval, int[] ptr, int lptr, int[] ind, int lind, int ibase, int fchk, int format, out int info)
{
	_mm_write(fname, matcode.ToCharArray(), nrow, ncol, nnz, val, lval, ptr, lptr, ind, lind, ibase, fchk, format, out info);
}

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_csr_check")]
public extern static void csr_check(int m, int n, double[] val, int[] rowptr, int[] colind, int[] result, out int info);

[DllImport(DLLSP, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_csx_check_sym")]
public extern static void csx_check_sym(int n, double[] val, int[] rowptr, int[] colind, out int info);

[DllImport(DLLPDE, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_fem2p")]
public extern static void _fem2p(int n, int ne, double[] x, double[] y, int ldknc, int[] knc, double[] p, double[] q, double[] f, int nb1, int[] ib, double[] bv, int nb2, int ldks2, int[] ks2, int ldalpha, double[] alpha, int ldbeta, double[] beta, int mxnnz, double[] a, int[] ia, int[] ja, int index, double[] b, int[] iwork, out int info);

public static void fem2p(int n, int ne, double[] x, double[] y, int ldknc, int[] knc, double[] p, double[] q, double[] f, int nb1, int[] ib, double[] bv, int nb2, int ldks2, int[] ks2, int ldalpha, double[] alpha, int ldbeta, double[] beta, int mxnnz, double[] a, int[] ia, int[] ja, int index, double[] b, out int info)
{
	int[] iwork;
	iwork = new int[2*n];
	_fem2p(n, ne, x, y, ldknc, knc, p, q, f, nb1, ib, bv, nb2, ldks2, ks2, ldalpha, alpha, ldbeta, beta, mxnnz, a, ia, ja, index, b, iwork, out info);
}

[DllImport(DLLPDE, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_mesh23")]
public extern static void mesh23(int nx, int ny, out int n, double sclx, double scly, double[] x, double[] y, out int ne, int ldknc, int[] knc, int ldks, int[] ks, int[] lb, out int nb, out int info);

[DllImport(DLLPDE, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_readgmsh22")]
public extern static void readgmsh22(string fname, out int n, double[] x, double[] y, double[] z, out int ne, int ldkc, int[] kc, int ldlb, int[] lb, out int info);

[DllImport(DLLPDE, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_writegmsh22")]
public extern static void writegmsh22(string fname, int n, double[] x, double[] y, double[] z, int ne, int ldkc, int[] kc, int ldlb, int[] lb, out int info);

[DllImport(DLLPDE, CallingConvention = CallingConvention.Cdecl, EntryPoint = "_writevtkug")]
public extern static void writevtkug(string fname, int n, double[] x, double[] y, double[] z, int ne, int ldkc, int[] kc, double[] u, out int info);

}
