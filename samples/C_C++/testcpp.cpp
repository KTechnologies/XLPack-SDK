/****************************************
 *                                      *
 *  C/C++ interface to XLPack           *
 *  Test program                        *
 *  Version 6.0 (June 14, 2022)         *
 *  (C) 2014-2022  K Technologies       *
 *                                      *
 ****************************************/

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <cerrno>
using namespace std;

#include "cnumlib"
#include "lapacke.h"

void test_sf()
{
	double x, x2, y, nu;

	cout << "** d1num\n";
	cout << "i = 1: " << d1num(1) << endl;
	cout << "i = 2: " << d1num(2) << endl;
	cout << "i = 3: " << d1num(3) << endl;
	cout << "i = 4: " << d1num(4) << endl;

	cout << "** factorial\n";
	x = 10;
	errno = 0;
	y = factorial(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** li\n";
	x = 2;
	errno = 0;
	y = li(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** ei\n";
	x = 1;
	errno = 0;
	y = ei(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** e1\n";
	x = 1;
	errno = 0;
	y = e1(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** digamma\n";
	x = 3.5;
	errno = 0;
	y = digamma(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besj0\n";
	x = 1;
	errno = 0;
	y = besj0(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besj1\n";
	x = 1;
	errno = 0;
	y = besj1(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besjnu\n";
	nu = 1; x = 1;
	errno = 0;
	y = besjnu(nu, x);
	cout << "nu = " << nu << ", x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besy0\n";
	x = 1;
	errno = 0;
	y = besy0(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besy1\n";
	x = 1;
	errno = 0;
	y = besy1(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besynu\n";
	nu = 1; x = 1;
	errno = 0;
	y = besynu(nu, x);
	cout << "nu = " << nu << ", x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besi0\n";
	x = 1;
	errno = 0;
	y = besi0(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;
	errno = 0;
	y = besi0e(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besi1\n";
	x = 1;
	errno = 0;
	y = besi1(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;
	errno = 0;
	y = besi1e(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besinu\n";
	nu = 1; x = 1;
	errno = 0;
	y = besinu(nu, x);
	cout << "nu = " << nu << ", x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besk0\n";
	x = 1;
	y = besk0(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;
	y = besk0e(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besk1\n";
	x = 1;
	errno = 0;
	y = besk1(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;
	errno = 0;
	y = besk1e(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** besknu\n";
	nu = 1; x = 1;
	errno = 0;
	y = besknu(nu, x);
	cout << "nu = " << nu << ", x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** celli1\n";
	x = 0.5;
	errno = 0;
	y = celli1(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** celli2\n";
	x = 0.5;
	errno = 0;
	y = celli2(x);
	cout << "x = " << x << ", y = " << y << ", errno = " << errno << endl;

	cout << "** celli3\n";
	x = 0.7; x2 = 0.5;
	errno = 0;
	y = celli3(x, x2);
	cout << "n = " << x << ", k = " << x2 << ", y = " << y << ", errno = " << errno << endl;
}

void test_dconst()
{
	double x;

	cout << "** dconst\n";
	for (int i = 0; i <= 35; i++) {
		x = dconst(i);
		cout << i << ": " << x << endl;
	}
}

void test_dgesv()
{
	const int n = 3;
	double a[n][n] = {
		{ 0.2, -0.32, -0.8 },
		{ -0.11, 0.81, -0.92 },
		{ -0.93, 0.37, -0.29 }
	};
	double b[] = { -0.3727, 0.4319, -1.4247 };
	double anorm, rcond = 0;
	int nrhs, lda, ldb, info;
	int ipiv[n];

	nrhs = 1; lda = n; ldb = n;
	anorm = LAPACKE_dlange(LAPACK_COL_MAJOR, '1', n, n, (double *)a, lda);
	info = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, (double *)a, lda, ipiv, b, ldb);
	if (info == 0)
		info = LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', n, (double *)a, lda, anorm, &rcond);
	cout << "** dgesv\n";
	cout << "x = " << b[0] << ", " << b[1] << ", " << b[2] << endl;
	cout << "rcond = " << rcond << ", info = " << info << endl;
}

void test_dposv()
{
	const int n = 3;
	double a[n][n] = {
		{ 2.2, 0.0, 0.0 },
		{ -0.11, 2.93, 0.0 },
		{ -0.32, 0.81, 2.37 }
	};
	double b[] = { -1.566, -2.8425, -1.1765 };
	double anorm, rcond = 0;
	int nrhs, lda, ldb, info;

	nrhs = 1; lda = n; ldb = n;
	anorm = LAPACKE_dlansy(LAPACK_COL_MAJOR, '1', 'U', n, (double *)a, lda);
	info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', n, nrhs, (double *)a, lda, b, ldb);
	if (info == 0)
		info = LAPACKE_dpocon(LAPACK_COL_MAJOR, 'U', n, (double *)a, lda, anorm, &rcond);
	cout << "** dposv\n";
	cout << "x = " << b[0] << ", " << b[1] << ", " << b[2] << endl;
	cout << "rcond = " << rcond << ", info = " << info << endl;
}

void test_dsyev()
{
	const int n = 3;
	double a[n][n] = {
		{ 2.2, 0.0, 0.0 },
		{ -0.11, 2.93, 0.0 },
		{ -0.32, 0.81, 2.37 }
	};
	double w[n];
	int lda, info;

	lda = n;
	info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, (double *)a, lda, w);
	cout << "** dsyev\n";
	cout << "Eigenvalues =\n";
	cout << w[0] << ", " << w[1] << ", " << w[2] << endl;
	cout << "Eigenvectors =\n";
	cout << a[0][0] << ", " << a[1][0] << ", " << a[2][0] << endl;
	cout << a[0][1] << ", " << a[1][1] << ", " << a[2][1] << endl;
	cout << a[0][2] << ", " << a[1][2] << ", " << a[2][2] << endl;
	cout << "info = " << info << endl;
}

void test_dgels()
{
	const int m = 4, n = 2;
	double x[] = { 0.2, 118.2, 337.4, 884.6 };
	double y[] = { 0.1, 118.1, 338.8, 888.0 };
	double a[n][m], ci[n];
	double s;
	int nrhs, lda, ldy, info;

	nrhs = 1; lda = m; ldy = m;
	for (int i = 0; i < m; i++) {
		a[0][i] = 1;
		a[1][i] = x[i];
	}
	info = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', m, n, nrhs, (double *)a, lda, y, ldy);
	cout << "** dgels\n";
	cout << "a0 = " << y[0] << ", a1 = " << y[1] << endl;
	cout << "info = " << info << endl;
	if (info == 0) {
		dgecov(0, n, lda, (double *)a, ci, info);
		s = 0.0;
		for (int i = n; i < m; i++)
			s = s + y[i]*y[i];
		s = s / (m - n);
		cout << "Std. dev. = " << sqrt(s*ci[0]) << ", " << sqrt(s*ci[1]) << endl;
		cout << "info = " << info << endl;
	}
}

void test_pchse()
{
	const int n = 4, ne = 2;
	double x[] = { 0.1, 0.11, 0.12, 0.13 };
	double y[] = { 2.3026, 2.2073, 2.1203, 2.0402 };
	double d[n], xe[ne], ye[ne], work[2*n];
	int incfd, skip, lwork = 2*n, info;

	incfd = 1; skip = 0;
	pchse(n, x, y, d, incfd, work, lwork, info);
	cout << "** pchse\n";
	cout << "info = " << info << endl;
	if (info == 0) {
		xe[0] = 0.115; xe[1] = 0.125;
		pchfe(n, x, y, d, incfd, skip, ne, xe, ye, info);
		cout << "ln(" << xe[0] << ") = " << ye[0] << endl;
		cout << "ln(" << xe[1] << ") = " << ye[1] << endl;
		cout << "info = " << info << endl;
	}
}

void test_pchia()
{
	const int n = 7;
	double x[] = { -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
	double y[] = { 0.5, 1.0, 0.5, 0.2, 0.1, 0.05882, 0.03846 };
	double d[n], work[2*n];
	double a, b, s;
	int incfd, skip, lwork = 2*n, info;

	incfd = 1; skip = 0;
	pchse(n, x, y, d, incfd, work, lwork, info);
	cout << "** pchse\n";
	cout << "info = " << info << endl;
	if (info == 0) {
		a = 0; b = 4;
		s = pchia(n, x, y, d, incfd, skip, a, b, info);
		cout << "** pchia\n";
		cout << "S = " << s << ", S(true) = " << atan(4.0) << endl;
		cout << "info = " << info << endl;
	}
}

void test_rpzero2()
{
	const int n = 5;
	double a[] = { 1.0, 0.0, 2.0, 2.0, -15.0, 10.0 };
	double zr[n], zi[n], s[n], work[8*n + 6];
	int iflag, maxiter, iter, info;

	iflag = 0;
	maxiter = 100;
	rpzero2(n, a, zr, zi, iflag, maxiter, iter, s, work, info);
	cout << "** rpzero2\n";
	for (int i = 0; i < n; i++)
		cout << zr[i] << "  " << zi[i] << "  " << s[i] << endl;
	cout << "iter = " << iter << ", info = " << info << endl;
}

double f_dfzero(double x)
{
	return (x*x - 2)*x - 5;
}

void test_dfzero()
{
	double b, c, r, re, ae;
	int info;

	b = 1.0; c = 3.0; r = b;
	re = 0.0; ae = 0.0;
	dfzero(f_dfzero, b, c, r, re, ae, info);
	cout << "** dfzero\n";
	cout << "x = " << b << ", info = " << info << endl;
}

void test_dfzero_r()
{
	double b, c, r, re, ae;
	double xx, yy = 0;
	int info, irev;

	b = 1.0; c = 3.0; r = b;
	re = 0.0; ae = 0.0;
	irev = 0;
	do {
		dfzero_r(b, c, r, re, ae, info, xx, yy, irev);
		if (irev != 0)
			yy = f_dfzero(xx);
	} while (irev != 0);
	cout << "** dfzero_r\n";
	cout << "x = " << b << ", info = " << info << endl;
}

void f_hybrd1(int n, double x[], double fvec[], int& iflag)
{
	fvec[0] = 4*x[0]*x[0] + x[1]*x[1] - 16;
	fvec[1] = x[0]*x[0] + x[1]*x[1] - 9;
}

void test_hybrd1()
{
	const int n = 2, lwork = n*(3*n + 17)/2;
	double x[n], fvec[n], work[lwork];
	double xtol = 1.0e-10;
	int info;

	x[0] = 1; x[1] = 2;
	hybrd1(f_hybrd1, n, x, fvec, xtol, work, lwork, info);
	cout << "** hybrd1\n";
	cout << "x[1] = " << x[0] << ", x[2] =" << x[1] << endl;
	cout << "info = " << info << endl;
}

void test_hybrd1_r()
{
	const int n = 2, lwork = n*(3*n + 17)/2;
	double x[n], fvec[n], work[lwork];
	double xtol = 1.0e-10;
	int info;
	double xx[n], yy[n];
	int iflag, irev;

	x[0] = 1; x[1] = 2;
	irev = 0;
	do {
		hybrd1_r(n, x, fvec, xtol, work, lwork, info, xx, yy, irev);
		if (irev >= 1 && irev <= 4) {
			iflag = 1;
			f_hybrd1(n, xx, yy, iflag);
		}
	} while (irev != 0);
	cout << "** hybrd1_r\n";
	cout << "x[1] = " << x[0] << ", x[2] =" << x[1] << endl;
	cout << "info = " << info << endl;
}

double f_dfmin(double x)
{
	return (x*x - 2)*x - 5;
}

void test_dfmin()
{
	double a, b, tol, x;

	a = 0; b = 1;
	tol = 1.0e-8;
	x = dfmin(a, b, f_dfmin, tol);
	cout << "** dfmin\n";
	cout << "x = " << x << endl;
}

void test_dfmin_r()
{
	double a, b, tol;
	double xx, yy = 0;
	int irev;

	a = 0; b = 1;
	tol = 1.0e-8;
	irev = 0;
	do {
		dfmin_r(a, b, tol, xx, yy, irev);
		if (irev != 0)
			yy = f_dfmin(xx);
	} while (irev != 0);
	cout << "** dfmin_r\n";
	cout << "x = " << xx << endl;
}

void f_optif0(int n, double x[], double& fval)
{
	// Rosenbrock function
	fval = 100*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]) + (1 - x[0])*(1 - x[0]);
}

void test_optif0()
{
	const int n = 2, lwork = n*(n + 11);
	double x[n], xpls[n], work[lwork];
	double fpls;
	int info;

	x[0] = -1.2; x[1] = 1.0;
	optif0(n, x, f_optif0, xpls, fpls, work, lwork, info);
	cout << "** optif0\n";
	cout << "xpls[1] = " << xpls[0] << ", xpls[2] = " << xpls[1] << endl;
	cout << "fpls = " << fpls << endl;
	cout << "info = " << info << endl;
}

void test_optif0_r()
{
	const int n = 2, lwork = n*(n + 11);
	double x[n], xpls[n], work[lwork];
	double fpls;
	int info;
	double xx[n], yy = 0;
	int irev;

	x[0] = -1.2; x[1] = 1.0;
	irev = 0;
	do {
		optif0_r(n, x, xpls, fpls, work, lwork, info, xx, yy, irev);
		if (irev >= 1 && irev <= 20)
			f_optif0(n, xx, yy);
	} while (irev != 0);
	cout << "** optif0_r\n";
	cout << "xpls[1] = " << xpls[0] << ", xpls[2] = " << xpls[1] << endl;
	cout << "fpls = " << fpls << endl;
	cout << "info = " << info << endl;
}

double f_qk15(double x)
{
	return 1/(1 + x*x);
}

void test_qk15()
{
	double a, b, result, abserr, resabs, resasc;

	// int(1/(1+x**2)) [0, 4] = atan(4)
	a = 0; b = 4;
	qk15(f_qk15, a, b, result, abserr, resabs, resasc);
	cout << "** qk15\n";
	cout << "result = " << result << ", abserr = " << abserr << endl;
}

void test_qk15_r()
{
	double a, b, result, abserr, resabs, resasc;
	double xx, yy = 0; int irev;

	// int(1/(1+x**2)) [0, 4] = atan(4)
	a = 0; b = 4;
	irev = 0;
	do {
		qk15_r(a, b, result, abserr, resabs, resasc, xx, yy, irev);
		if (irev != 0)
			yy = f_qk15(xx);
	} while (irev != 0);
	cout << "** qk15_r\n";
	cout << "result = " << result << ", abserr = " << abserr << endl;
}

double f_qag(double x)
{
	return 1/(1 + x*x);
}

void test_qag()
{
	const int limit = 100, lwork = 4*limit, liwork = limit;
	double a, b, epsabs, epsrel, result, abserr, work[lwork];
	int key, neval, last, iwork[liwork], info;

	// int(1/(1+x**2)) [0, 4] = atan(4)
	a = 0; b = 4;
	epsabs = 1.0e-10; epsrel = 1.0e-10;
	for (int i = 1; i <= 6; i++) {
		key = i;
		qag(f_qag, a, b, epsabs, epsrel, key, limit, result, abserr, neval, last, work, lwork, iwork, info);
		cout << "** qag (key = " << key << ")\n";
		cout << "result = " << result << ", abserr = " << abserr << ", neval = " << neval << ", last = " << last << endl;
	}
}

void test_qag_r()
{
	const int limit = 100, lwork = 4*limit, liwork = limit;
	double a, b, epsabs, epsrel, result, abserr, work[lwork];
	int key, neval, last, iwork[liwork], info;
	double xx, yy = 0; int irev;

	// int(1/(1+x**2)) [0, 4] = atan(4)
	a = 0; b = 4;
	epsabs = 1.0e-10; epsrel = 1.0e-10;
	for (int i = 1; i <= 6; i++) {
		key = i;
		irev = 0;
		do {
			qag_r(a, b, epsabs, epsrel, key, limit, result, abserr, neval, last, work, lwork, iwork, info, xx, yy, irev);
			if (irev != 0)
				yy = f_qag(xx);
		} while (irev != 0);
		cout << "** qag_r (key = " << key << ")\n";
		cout << "result = " << result << ", abserr = " << abserr << ", neval = " << neval << ", last = " << last << endl;
	}
}

double f_qagi(double x)
{
	return 1/(1 + x*x);
}

void test_qagi()
{
	const int limit = 100, lwork = 4*limit, liwork = limit;
	double bound, epsabs, epsrel, result, abserr, work[lwork];
	int inf, neval, last, iwork[liwork], info;

	// int(1/(1+x**2)) [-inf, +inf] = pi
	// int(1/(1+x**2)) [0, +inf] = pi / 2
	// int(1/(1+x**2)) [-inf, 0] = pi / 2
	bound = 0;
	epsabs = 1.0e-10; epsrel = 1.0e-10;
	inf = 2;
	qagi(f_qagi, bound, inf, epsabs, epsrel, limit, result, abserr, neval, last, work, lwork, iwork, info);
	cout << "** qagi [-inf, +inf]\n";
	cout << "result = " << result << ", abserr = " << abserr << ", neval = " << neval << ", last = " << last << endl;
	inf = 1;
	qagi(f_qagi, bound, inf, epsabs, epsrel, limit, result, abserr, neval, last, work, lwork, iwork, info);
	cout << "** qagi [0, +inf]\n";
	cout << "result = " << result << ", abserr = " << abserr << ", neval = " << neval << ", last = " << last << endl;
	inf = -1;
	qagi(f_qagi, bound, inf, epsabs, epsrel, limit, result, abserr, neval, last, work, lwork, iwork, info);
	cout << "** qagi [-inf, 0]\n";
	cout << "result = " << result << ", abserr = " << abserr << ", neval = " << neval << ", last = " << last << endl;
}

void test_qagi_r()
{
	const int limit = 100, lwork = 4*limit, liwork = limit;
	double bound, epsabs, epsrel, result, abserr, work[lwork];
	int inf, neval, last, iwork[liwork], info;
	double xx, yy = 0; int irev;

	// int(1/(1+x**2)) [-inf, +inf] = pi
	// int(1/(1+x**2)) [0, +inf] = pi / 2
	// int(1/(1+x**2)) [-inf, 0] = pi / 2
	bound = 0;
	epsabs = 1.0e-10; epsrel = 1.0e-10;
	inf = 2;
	irev = 0;
	do {
		qagi_r(bound, inf, epsabs, epsrel, limit, result, abserr, neval, last, work, lwork, iwork, info, xx, yy, irev);
		if (irev != 0)
			yy = f_qagi(xx);
	} while (irev != 0);
	cout << "** qagi_r [-inf, +inf]\n";
	cout << "result = " << result << ", abserr = " << abserr << ", neval = " << neval << ", last = " << last << endl;
	inf = 1;
	irev = 0;
	do {
		qagi_r(bound, inf, epsabs, epsrel, limit, result, abserr, neval, last, work, lwork, iwork, info, xx, yy, irev);
		if (irev != 0)
			yy = f_qagi(xx);
	} while (irev != 0);
	cout << "** qagi_r [0, +inf]\n";
	cout << "result = " << result << ", abserr = " << abserr << ", neval = " << neval << ", last = " << last << endl;
	inf = -1;
	irev = 0;
	do {
		qagi_r(bound, inf, epsabs, epsrel, limit, result, abserr, neval, last, work, lwork, iwork, info, xx, yy, irev);
		if (irev != 0)
			yy = f_qagi(xx);
	} while (irev != 0);
	cout << "** qagi_r [-inf, 0]\n";
	cout << "result = " << result << ", abserr = " << abserr << ", neval = " << neval << ", last = " << last << endl;
}

static double alfasq;
static int neval;

void f_derkf(int n, double t, double y[], double yp[])
{
	double r;

	r = y[0]*y[0] + y[1]*y[1];
	r = r*sqrt(r)/alfasq;
	yp[0] = y[2];
	yp[1] = y[3];
	yp[2] = -y[0]/r;
	yp[3] = -y[1]/r;
	neval += 1;
}

void test_derkf()
{
	const int n = 4, lwork = 9*n + 20, liwork = 20;
	double t, tout, y[n];
	double work[lwork]; int iwork[liwork];
	double tfinal, tprint;
	int info;

	double ecc = 0.25, alfa = M_PI/4;
	double rtol = 1e-10, atol = rtol; int itol = 0;
	int mode = 0; /* Interval mode */

	alfasq = alfa*alfa;
	t = 0;
	y[0] = 1 - ecc;
	y[1] = 0;
	y[2] = 0;
	y[3] = alfa*sqrt((1 + ecc)/(1 - ecc));
	tfinal = 12;
	tprint = 1;
	cout << "** derkf\n";
	neval = 0;
	info = 0;
	do {
		tout = t + tprint;
		derkf(n, f_derkf, t, y, tout, &rtol, &atol, itol, mode, work, lwork, iwork, liwork, info);
		if (info != 1)
			break;
		cout << "t = " << t << ", y[1] = " << y[0] << ", y[2] = " << y[1] << ", y[3] = " << y[2] << ", y[4] = " << y[3] << ", neval = " << neval << ", info = " << info << endl;
	} while (t < tfinal);
	printf("info = %d\n", info);
}

void test_derkf_2()
{
	const int n = 4, lwork = 11*n + 20, liwork = 20;
	double t, tout, y[n];
	double work[lwork]; int iwork[liwork];
	double tfinal, tprint;
	int info;

	double ecc = 0.25, alfa = M_PI/4;
	double rtol = 1e-10, atol = rtol; int itol = 0;
	int mode = 2; // Step mode (dense output)

	alfasq = alfa*alfa;
	t = 0;
	y[0] = 1 - ecc;
	y[1] = 0;
	y[2] = 0;
	y[3] = alfa*sqrt((1 + ecc)/(1 - ecc));
	tfinal = 12;
	tprint = 1;
	cout << "** derkf (2)\n";
	neval = 0;
	tout = t + tprint;
	info = 0;
	do {
		derkf(n, f_derkf, t, y, tfinal, &rtol, &atol, itol, mode, work, lwork, iwork, liwork, info);
		if (info == 1 || info == 2) {
			while (t >= tout) {
				double y1[n];
				derkf_int(n, tout, y1, work);
		cout << "t = " << tout << ", y[1] = " << y1[0] << ", y[2] = " << y1[1] << ", y[3] = " << y1[2] << ", y[4] = " << y1[3] << ", neval = " << neval << ", info = " << info << endl;
				tout = tout + tprint;
			}
		} else
			break;
	} while (t < tfinal);
	printf("info = %d\n", info);
}

void test_derkf_r()
{
	const int n = 4, lwork = 7*n + 20, liwork = 20;
	double t, tout, y[n];
	double work[lwork]; int iwork[liwork];
	double tt, yy[n], yyp[n]; int irev;
	double tfinal, tprint;
	int info;

	double ecc = 0.25, alfa = M_PI/4;
	double rtol = 1e-10, atol = rtol; int itol = 0;
	int mode = 0; /* Interval mode */

	alfasq = alfa*alfa;
	t = 0;
	y[0] = 1 - ecc;
	y[1] = 0;
	y[2] = 0;
	y[3] = alfa*sqrt((1 + ecc)/(1 - ecc));
	tfinal = 12;
	tprint = 1;
	cout << "** derkf_r\n";
	neval = 0;
	info = 0;
	do {
		tout = t + tprint;
		irev = 0;
		do {
			derkf_r(n, t, y, tout, &rtol, &atol, itol, mode, work, lwork, iwork, liwork, info, tt, yy, yyp, irev);
			if (irev != 0)
				f_derkf(n, tt, yy, yyp);
		} while (irev != 0);
		if (info != 1)
			break;
		cout << "t = " << t << ", y[1] = " << y[0] << ", y[2] = " << y[1] << ", y[3] = " << y[2] << ", y[4] = " << y[3] << ", neval = " << neval << ", info = " << info << endl;
	} while (t < tfinal);
	printf("info = %d\n", info);
}

void test_derkf_r_2()
{
	const int n = 4, lwork = 11*n + 20, liwork = 20;
	double t, tout, y[n];
	double work[lwork]; int iwork[liwork];
	double tt, yy[n], yyp[n]; int irev;
	double tfinal, tprint;
	int info;

	double ecc = 0.25, alfa = M_PI/4;
	double rtol = 1e-10, atol = rtol; int itol = 0;
	int mode = 2; // Step mode (dense output)

	alfasq = alfa*alfa;
	t = 0;
	y[0] = 1 - ecc;
	y[1] = 0;
	y[2] = 0;
	y[3] = alfa*sqrt((1 + ecc)/(1 - ecc));
	tfinal = 12;
	tprint = 1;
	cout << "** derkf_r (2)\n";
	neval = 0;
	tout = t + tprint;
	info = 0;
	do {
		irev = 0;
		do {
			derkf_r(n, t, y, tfinal, &rtol, &atol, itol, mode, work, lwork, iwork, liwork, info, tt, yy, yyp, irev);
			if (irev != 0)
				f_derkf(n, tt, yy, yyp);
		} while (irev != 0);
		if (info == 1 || info == 2) {
			while (t >= tout) {
				double y1[n];
				derkf_int(n, tout, y1, work);
		cout << "t = " << tout << ", y[1] = " << y1[0] << ", y[2] = " << y1[1] << ", y[3] = " << y1[2] << ", y[4] = " << y1[3] << ", neval = " << neval << ", info = " << info << endl;
				tout = tout + tprint;
			}
		} else
			break;
	} while (t < tfinal);
	printf("info = %d\n", info);
}

void test_rfft1()
{
	const int n = 10, lwsave = 17, inc = 1, lr = inc*(n - 1) + 1, lwork = n;
	double wsave[lwsave], r[lr], rcopy[lr], diff, work[lwork];
	int info, i, k;
	unsigned int seed;

	// Initialization
	cout << "** rfft1\n";
	seed = 13;
	init_genrand(seed);
	rfft1i(n, wsave, lwsave, info);
	if (info != 0) {
		cout << "Error during initialization\n";
		return;
	}
	// Generate test data
	for (i = 0; i < n; i++) {
		k = inc*i;
		r[k] = genrand_res53();
		rcopy[k] = r[k];
	}
	// Forward transform
	rfft1f(n, inc, r, lr, wsave, lwsave, work, lwork, info);
	if (info != 0) {
		cout << "Error in rfft1f\n";
		return;
	}
	// Backward transform
	rfft1b(n, inc, r, lr, wsave, lwsave, work, lwork, info);
	if (info != 0) {
		cout << "Error in rfft1b\n";
		return;
	}
	// Check results
	diff = 0;
	for (i = 0; i < n; i++) {
		k = inc*i;
		if (fabs(r[k] - rcopy[k]) > diff)
			diff = fabs(r[k] - rcopy[k]);
		cout << rcopy[k] << ", " << r[k] << ", " << fabs(r[k] - rcopy[k]) << endl;
	}
	cout << "diff(max) = " << diff << endl;
}

void f_lmdif1(int m, int n, double x[], double fvec[], int& iflag)
{
	double xdata[] = { 77.6, 114.9, 141.1, 190.8, 239.9, 289.0, 332.8, 378.4,
		434.8, 477.3, 536.8, 593.1, 689.1, 760.0 };
	double ydata[] = { 10.07, 14.73, 17.94, 23.93, 29.61, 35.18, 40.02, 44.82,
		50.76, 55.05, 61.01,  66.4, 75.47, 81.78 };
	for (int i = 0; i < m; i++)
		fvec[i] = ydata[i] - x[0]*(1 - exp(-xdata[i]*x[1]));
}

void test_lmdif1()
{
	const int m = 14, n = 2, lwork = n*(m + 6) + 2*m, liwork = n;
	double x[n], fvec[m], work[lwork];
	double tol = 1.0e-10;
	int iwork[liwork], info;

	x[0] = 500; x[1] = 0.0001;
	lmdif1(f_lmdif1, m, n, x, fvec, tol, work, lwork, iwork, info);
	cout << "** lmdif1\n";
	cout << "x[1] = " << x[0] << ", x[2] = " << x[1] << endl;
	cout << "info = " << info << endl;
}

void test_lmdif1_r()
{
	const int m = 14, n = 2, lwork = n*(m + 6) + 2*m, liwork = n;
	double x[n], fvec[m], work[lwork];
	double tol = 1.0e-10;
	int iwork[liwork], info;
	double xx[n], yy[m]; int iflag, irev;

	x[0] = 500; x[1] = 0.0001;
	irev = 0;
	do {
		lmdif1_r(m, n, x, fvec, tol, work, lwork, iwork, info, xx, yy, irev);
		if (irev >= 1 && irev <= 3) {
			iflag = 1;
			f_lmdif1(m, n, xx, yy, iflag);
		}
	} while (irev != 0);
	cout << "** lmdif1_r\n";
	cout << "x[1] = " << x[0] << ", x[2] = " << x[1] << endl;
	cout << "info = " << info << endl;
}

void test_rand()
{
	unsigned int seed, r32;
	int r31;
	double r53;

	seed = 11;
	cout << "** Random numbers: seed = " << seed <<  endl;
	init_genrand(seed);
	for (int i = 0; i < 10; i++) {
		r32 = genrand_int32();
		r31 = genrand_int31();
		r53 = genrand_res53();
		cout << r32 << ", " << r31 << ", " << r53 << endl;
	}
}

void test_dlamch()
{
	cout << "** dlamch\n";
	cout << "e: " << dlamch('e') << endl;
	cout << "s: " << dlamch('s') << endl;
	cout << "b: " << dlamch('b') << endl;
	cout << "p: " << dlamch('p') << endl;
	cout << "n: " << dlamch('n') << endl;
	cout << "r: " << dlamch('r') << endl;
	cout << "m: " << dlamch('m') << endl;
	cout << "u: " << dlamch('u') << endl;
	cout << "l: " << dlamch('l') << endl;
	cout << "o: " << dlamch('o') << endl;
}

int main()
{

	test_sf();
	test_dconst();
	test_dgesv();
	test_dposv();
	test_dsyev();
	test_dgels();
	test_pchse();
	test_pchia();
	test_rpzero2();
	test_dfzero();
	test_dfzero_r();
	test_hybrd1();
	test_hybrd1_r();
	test_dfmin();
	test_dfmin_r();
	test_optif0();
	test_optif0_r();
	test_qk15();
	test_qk15_r();
	test_qag();
	test_qag_r();
	test_qagi();
	test_qagi_r();
	test_derkf();
	test_derkf_r();
	test_derkf_2();
	test_derkf_r_2();
	test_rfft1();
	test_lmdif1();
	test_lmdif1_r();
	test_rand(); test_rand();
	test_dlamch();
	return 0;
}
