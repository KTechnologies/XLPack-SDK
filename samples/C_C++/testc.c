/****************************************
 *                                      *
 *  C/C++ interface to XLPack           *
 *  Test program                        *
 *  Version 7.0 (January 31, 2023)      *
 *  (C) 2014-2023  K Technologies       *
 *                                      *
 ****************************************/

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <errno.h>

#include "cnumlib.h"
#include "lapacke.h"

void test_sf()
{
	double x, x2, y, nu;

	printf("** d1num\n");
	printf("i = 1: %g\n", d1num(1));
	printf("i = 2: %g\n", d1num(2));
	printf("i = 3: %g\n", d1num(3));
	printf("i = 4: %g\n", d1num(4));

	printf("** factorial\n");
	x = 10;
	errno = 0;
	y = factorial(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** li\n");
	x = 2;
	errno = 0;
	y = li(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** ei\n");
	x = 1;
	errno = 0;
	y = ei(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** e1\n");
	x = 1;
	errno = 0;
	y = e1(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** digamma\n");
	x = 3.5;
	errno = 0;
	y = digamma(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** besj0\n");
	x = 1;
	errno = 0;
	y = besj0(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** besj1\n");
	x = 1;
	errno = 0;
	y = besj1(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** besjnu\n");
	nu = 1; x = 1;
	errno = 0;
	y = besjnu(nu, x);
	printf("nu = %g, x = %g, y = %g, errno = %d\n", nu, x, y, errno);

	printf("** besy0\n");
	x = 1;
	errno = 0;
	y = besy0(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** besy1\n");
	x = 1;
	errno = 0;
	y = besy1(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** besynu\n");
	nu = 1; x = 1;
	errno = 0;
	y = besynu(nu, x);
	printf("nu = %g, x = %g, y = %g, errno = %d\n", nu, x, y, errno);

	printf("** besi0\n");
	x = 1;
	errno = 0;
	y = besi0(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);
	errno = 0;
	y = besi0e(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** besi1\n");
	x = 1;
	errno = 0;
	y = besi1(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);
	errno = 0;
	y = besi1e(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** besinu\n");
	nu = 1; x = 1;
	errno = 0;
	y = besinu(nu, x);
	printf("nu = %g, x = %g, y = %g, errno = %d\n", nu, x, y, errno);

	printf("** besk0\n");
	x = 1;
	errno = 0;
	y = besk0(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);
	errno = 0;
	y = besk0e(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** besk1\n");
	x = 1;
	errno = 0;
	y = besk1(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);
	errno = 0;
	y = besk1e(x);
	printf("x = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** besknu\n");
	nu = 1; x = 1;
	errno = 0;
	y = besknu(nu, x);
	printf("nu = %g, x = %g, y = %g, errno = %d\n", nu, x, y, errno);

	printf("** celli1\n");
	x = 0.5;
	errno = 0;
	y = celli1(x);
	printf("k = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** celli2\n");
	x = 0.5;
	errno = 0;
	y = celli2(x);
	printf("k = %g, y = %g, errno = %d\n", x, y, errno);

	printf("** celli3\n");
	x = 0.7; x2 = 0.5;
	errno = 0;
	y = celli3(x, x2);
	printf("n = %g, k = %g, y = %g, errno = %d\n", x, x2, y, errno);
}

void test_dconst()
{
	double x;

	printf("** dconst\n");
	for (int i = 0; i <= 35; i++) {
		x = dconst(i);
		printf("%d: %g\n", i, x);
	}
}

/* Parameters for test_dgesv */
#define N	3

void test_dgesv()
{
	double a[N][N] = {
		{ 0.2, -0.32, -0.8 },
		{ -0.11, 0.81, -0.92 },
		{ -0.93, 0.37, -0.29 }
	};
	double b[] = { -0.3727, 0.4319, -1.4247 };
	double anorm, rcond = 0;
	int nrhs, lda, ldb, info;
	int ipiv[N];

	nrhs = 1; lda = N; ldb = N;
	anorm = LAPACKE_dlange(LAPACK_COL_MAJOR, '1', N, N, (double *)a, lda);
	info = LAPACKE_dgesv(LAPACK_COL_MAJOR, N, nrhs, (double *)a, lda, ipiv, b, ldb);
	if (info == 0)
		info = LAPACKE_dgecon(LAPACK_COL_MAJOR, '1', N, (double *)a, lda, anorm, &rcond);
	printf("** dgesv\n");
	printf("x = %g  %g  %g\n", b[0], b[1], b[2]);
	printf("rcond = %g, info = %d\n", rcond, info);
}

#undef N

/* Parameters for test_dposv */
#define N	3

void test_dposv()
{
	double a[N][N] = {
		{ 2.2, 0.0, 0.0 },
		{ -0.11, 2.93, 0.0 },
		{ -0.32, 0.81, 2.37 }
	};
	double b[] = { -1.566, -2.8425, -1.1765 };
	double anorm, rcond = 0;
	int nrhs, lda, ldb, info;

	nrhs = 1; lda = N; ldb = N;
	anorm = LAPACKE_dlansy(LAPACK_COL_MAJOR, '1', 'U', N, (double *)a, lda);
	info = LAPACKE_dposv(LAPACK_COL_MAJOR, 'U', N, nrhs, (double *)a, lda, b, ldb);
	if (info == 0)
		info = LAPACKE_dpocon(LAPACK_COL_MAJOR, 'U', N, (double *)a, lda, anorm, &rcond);
	printf("** dposv\n");
	printf("x =\n");
	printf("%g  %g  %g\n", b[0], b[1], b[2]);
	printf("rcond = %g, info = %d\n", rcond, info);
}

#undef N

/* Parameters for test_dsyev */
#define N	3

void test_dsyev()
{
	double a[N][N] = {
		{ 2.2, 0.0, 0.0 },
		{ -0.11, 2.93, 0.0 },
		{ -0.32, 0.81, 2.37 }
	};
	double w[N];
	int n, lda, info;

	lda = N;
	info = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', N, (double *)a, lda, w);
	printf("** dsyev\n");
	printf("Eigenvalues =\n");
	printf("%g  %g  %g\n", w[0], w[1], w[2]);
	printf("Eigenvectors =\n");
	printf("%g  %g  %g\n", a[0][0], a[1][0], a[2][0]);
	printf("%g  %g  %g\n", a[0][1], a[1][1], a[2][1]);
	printf("%g  %g  %g\n", a[0][2], a[1][2], a[2][2]);
    printf("info = %d\n", info);
}

#undef N

/* Parameters for test_dgels */
#define M	4
#define N	2

void test_dgels()
{
	double x[] = { 0.2, 118.2, 337.4, 884.6 };
	double y[] = { 0.1, 118.1, 338.8, 888.0 };
	double a[N][M], ci[N];
	double s;
	int nrhs, lda, ldy, info;

	nrhs = 1; lda = M; ldy = M;
	for (int i = 0; i < M; i++) {
		a[0][i] = 1;
		a[1][i] = x[i];
	}
	info = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', M, N, nrhs, (double *)a, lda, y, ldy);
	printf("** dgels\n");
	printf("a0 = %g, a1 = %g\n", y[0], y[1]);
	printf("info = %d\n", info);
	if (info == 0) {
#if defined(_WIN32)
		dgecov(0, N, lda, (double *)a, ci, &info);
#else
		dgecov(0, N, lda, a, ci, &info);
#endif
		s = 0.0;
		for (int i = N; i < M; i++)
			s = s + y[i]*y[i];
		s = s / (M - N);
		printf("Std. dev. = %g, %g\n", sqrt(s*ci[0]), sqrt(s*ci[1]));
		printf("info = %d\n", info);
	}
}

#undef M
#undef N

/* Parameters for test_pchse */
#define N	4
#define NE	2
#define LWORK	(2*N)

void test_pchse()
{
	double x[] = { 0.1, 0.11, 0.12, 0.13 };
	double y[] = { 2.3026, 2.2073, 2.1203, 2.0402 };
	double d[N], xe[NE], ye[NE], work[LWORK];
	int incfd, skip, lwork = LWORK, info;

	incfd = 1; skip = 0;
	pchse(N, x, y, d, incfd, work, lwork, &info);
	printf("** pchse\n");
	printf("info = %d\n", info);
	if (info == 0) {
		xe[0] = 0.115; xe[1] = 0.125;
		pchfe(N, x, y, d, incfd, skip, NE, xe, ye, &info);
		printf("ln(%g) = %g\n", xe[0], ye[0]);
		printf("ln(%g) = %g\n", xe[1], ye[1]);
		printf("info = %d\n", info);
	}
}

#undef N
#undef NE
#undef LWORK

/* Parameters for test_pchia */
#define N	7
#define LWORK	(2*N)

void test_pchia()
{
	double x[] = { -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
	double y[] = { 0.5, 1.0, 0.5, 0.2, 0.1, 0.05882, 0.03846 };
	double d[N], work[LWORK];
	double a, b, s;
	int incfd, skip, lwork = LWORK, info;

	incfd = 1; skip = 0;
	pchse(N, x, y, d, incfd, work, lwork, &info);
	printf("** pchse\n");
	printf("info = %d\n", info);
	if (info == 0) {
		a = 0; b = 4;
		s = pchia(N, x, y, d, incfd, skip, a, b, &info);
		printf("** pchia\n");
		printf("S = %g, S(true) = %g\n", s, atan(4.0));
		printf("info = %d\n", info);
	}
}

#undef N
#undef LWORK

/* Parameters for test_rpzero2 */
#define N	5
#define LWORK	(8*N + 6)

void test_rpzero2()
{
	double a[] = { 1.0, 0.0, 2.0, 2.0, -15.0, 10.0 };
	double zr[N], zi[N], s[N], work[LWORK];
	int n, iflag, maxiter, iter, info;

	iflag = 0;
	maxiter = 100;
	rpzero2(N, a, zr, zi, iflag, maxiter, &iter, s, work, &info);
	printf("** rpzero2\n");
	for (int i = 0; i < N; i++)
		printf("%g  %g  %g\n", zr[i], zi[i], s[i]);
	printf("iter = %d, info = %d\n", iter, info);
}

#undef N
#undef LWORK

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
	dfzero(f_dfzero, &b, &c, r, re, ae, &info);
	printf("** dfzero\n");
	printf("x = %g, info = %d\n", b, info);
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
		dfzero_r(&b, &c, r, re, ae, &info, &xx, yy, &irev);
		if (irev != 0)
			yy = f_dfzero(xx);
	} while (irev != 0);
	printf("** dfzero_r\n");
	printf("x = %g, info = %d\n", b, info);
}

/* Parameters for test_hybrd1 */
#define N	2
#define LWORK	(N*(3*N + 17)/2)

void f_hybrd1(int n, double x[], double fvec[], int *iflag)
{
	fvec[0] = 4*x[0]*x[0] + x[1]*x[1] - 16;
	fvec[1] = x[0]*x[0] + x[1]*x[1] - 9;
}

void test_hybrd1()
{
	double x[N], fvec[N], work[LWORK];
	double xtol = 1.0e-10;
	int n = N, lwork = LWORK, info;

	x[0] = 1; x[1] = 2;
	hybrd1(f_hybrd1, n, x, fvec, xtol, work, lwork, &info);
	printf("** hybrd1\n");
	printf("x[1] = %g, x[2] = %g\n", x[0], x[1]);
	printf("info = %d\n", info);
}

void test_hybrd1_r()
{
	double x[N], fvec[N], work[LWORK];
	double xtol = 1.0e-10;
	int n = N, lwork = LWORK, info;
	double xx[N], yy[N];
	int iflag, irev;

	x[0] = 1; x[1] = 2;
	irev = 0;
	do {
		hybrd1_r(n, x, fvec, xtol, work, lwork, &info, xx, yy, &irev);
		if (irev >= 1 && irev <= 4) {
			iflag = 1;
			f_hybrd1(n, xx, yy, &iflag);
		}
	} while (irev != 0);
	hybrd1(f_hybrd1, n, x, fvec, xtol, work, lwork, &info);
	printf("** hybrd1_r\n");
	printf("x[1] = %g, x[2] = %g\n", x[0], x[1]);
	printf("info = %d\n", info);
}

#undef N
#undef LWORK

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
	printf("** dfmin\n");
	printf("x = %g\n", x);
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
		dfmin_r(a, b, tol, &xx, yy, &irev);
		if (irev != 0)
			yy = f_dfmin(xx);
	} while (irev != 0);
	printf("** dfmin_r\n");
	printf("x = %g\n", xx);
}

/* Parameters for test_optif0 */
#define N	2
#define LWORK	(N*(N + 11))

void f_optif0(int n, double x[], double *fval)
{
	/* Rosenbrock function */
	*fval = 100*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]) + (1 - x[0])*(1 - x[0]);
}

void test_optif0()
{
	double x[N], xpls[N], work[LWORK];
	double fpls;
	int n, lwork = LWORK, info;

	n = N;
	x[0] = -1.2; x[1] = 1.0;
	optif0(n, x, f_optif0, xpls, &fpls, work, lwork, &info);
	printf("** optif0\n");
	printf("xpls[1] = %g, xpls[2] = %g\n", xpls[0], xpls[1]);
	printf("fpls = %g\n", fpls);
	printf("info = %d\n", info);
}

void test_optif0_r()
{
//	double x[N], xpls[N], work[LWORK];
	double x[2], xpls[2], work[100];
	double fpls;
	int n, lwork = LWORK, info;
	double xx[N], yy = 0;
	int irev;

	n = N;
	x[0] = -1.2; x[1] = 1.0;
	irev = 0;
	do {
		optif0_r(n, x, xpls, &fpls, work, lwork, &info, xx, yy, &irev);
		if (irev >= 1 && irev <= 20)
			f_optif0(n, xx, &yy);
	} while (irev != 0);
	printf("** optif0_r\n");
	printf("xpls[1] = %g, xpls[2] = %g\n", xpls[0], xpls[1]);
	printf("fpls = %g\n", fpls);
	printf("info = %d\n", info);
}

#undef N
#undef LWORK

double f_qk15(double x)
{
	return 1/(1 + x*x);
}

void test_qk15()
{
	double a, b, result, abserr, resabs, resasc;

	/* int(1/(1+x**2)) [0, 4] = atan(4) */
	a = 0; b = 4;
	qk15(f_qk15, a, b, &result, &abserr, &resabs, &resasc);
	printf("** qk15\n");
	printf("result = %g, abserr = %g\n", result, abserr);
}

void test_qk15_r()
{
	double a, b, result, abserr, resabs, resasc;
	double xx, yy = 0; int irev;

	/* int(1/(1+x**2)) [0, 4] = atan(4) */
	a = 0; b = 4;
	irev = 0;
	do {
		qk15_r(a, b, &result, &abserr, &resabs, &resasc, &xx, yy, &irev);
		if (irev != 0)
			yy = f_qk15(xx);
	} while (irev != 0);
	printf("** qk15_r\n");
	printf("result = %g, abserr = %g\n", result, abserr);
}

/* Parameters for test_qag */
#define LIMIT	100
#define LWORK	(4*LIMIT)
#define LIWORK	LIMIT

double f_qag(double x)
{
	return 1/(1 + x*x);
}

void test_qag()
{
	double a, b, epsabs, epsrel, result, abserr, work[LWORK];
	int key, limit, neval, last, iwork[LIWORK], lwork = LWORK, info;

	/* int(1/(1+x**2)) [0, 4] = atan(4) */
	a = 0; b = 4;
	epsabs = 1.0e-10; epsrel = 1.0e-10;
	limit = LIMIT;
	for (int i = 1; i <= 6; i++) {
		key = i;
		qag(f_qag, a, b, epsabs, epsrel, key, limit, &result, &abserr, &neval, &last, work, lwork, iwork, &info);
		printf("** qag (key = %d)\n", key);
		printf("result = %g, abserr = %g, neval = %d, last = %d\n", result, abserr, neval, last);
	}
}

void test_qag_r()
{
	double a, b, epsabs, epsrel, result, abserr, work[LWORK];
	int key, limit, neval, last, iwork[LIWORK], lwork = LWORK, info;
	double xx, yy = 0; int irev;

	/* int(1/(1+x**2)) [0, 4] = atan(4) */
	a = 0; b = 4;
	epsabs = 1.0e-10; epsrel = 1.0e-10;
	limit = LIMIT;
	for (int i = 1; i <= 6; i++) {
		key = i;
		irev = 0;
		do {
			qag_r(a, b, epsabs, epsrel, key, limit, &result, &abserr, &neval, &last, work, lwork, iwork, &info, &xx, yy, &irev);
			if (irev != 0)
				yy = f_qag(xx);
		} while (irev != 0);
		printf("** qag_r (key = %d)\n", key);
		printf("result = %g, abserr = %g, neval = %d, last = %d\n", result, abserr, neval, last);
	}
}

#undef N
#undef LWORK
#undef LIWORK

/* Parameters for test_qagi */
#define LIMIT	100
#define LWORK	(4*LIMIT)
#define LIWORK	LIMIT

double f_qagi(double x)
{
	return 1/(1 + x*x);
}

void test_qagi()
{
	double bound, epsabs, epsrel, result, abserr, work[LWORK];
	int inf, limit, neval, last, iwork[LIWORK], lwork = LWORK, info;

	/* int(1/(1+x**2)) [-inf, +inf] = pi */
	/* int(1/(1+x**2)) [0, +inf] = pi / 2 */
	/* int(1/(1+x**2)) [-inf, 0] = pi / 2 */
	bound = 0;
	epsabs = 1.0e-10; epsrel = 1.0e-10;
	limit = LIMIT;
	inf = 2;
	qagi(f_qagi, bound, inf, epsabs, epsrel, limit, &result, &abserr, &neval, &last, work, lwork, iwork, &info);
	printf("** qagi [-inf, +inf]\n");
	printf("result = %g, abserr = %g, neval = %d, last = %d\n", result, abserr, neval, last);
	inf = 1;
	qagi(f_qagi, bound, inf, epsabs, epsrel, limit, &result, &abserr, &neval, &last, work, lwork, iwork, &info);
	printf("** qagi [0, +inf]\n");
	printf("result = %g, abserr = %g, neval = %d, last = %d\n", result, abserr, neval, last);
	inf = -1;
	qagi(f_qagi, bound, inf, epsabs, epsrel, limit, &result, &abserr, &neval, &last, work, lwork, iwork, &info);
	printf("** qagi [-inf, 0]\n");
	printf("result = %g, abserr = %g, neval = %d, last = %d\n", result, abserr, neval, last);
}

void test_qagi_r()
{
	double bound, epsabs, epsrel, result, abserr, work[LWORK];
	int inf, limit, neval, last, iwork[LIWORK], lwork = LWORK, info;
	double xx, yy = 0; int irev;

	/* int(1/(1+x**2)) [-inf, +inf] = pi */
	/* int(1/(1+x**2)) [0, +inf] = pi / 2 */
	/* int(1/(1+x**2)) [-inf, 0] = pi / 2 */
	bound = 0;
	epsabs = 1.0e-10; epsrel = 1.0e-10;
	limit = LIMIT;
	inf = 2;
	irev = 0;
	do {
		qagi_r(bound, inf, epsabs, epsrel, limit, &result, &abserr, &neval, &last, work, lwork, iwork, &info, &xx, yy, &irev);
		if (irev != 0)
			yy = f_qagi(xx);
	} while (irev != 0);
	printf("** qagi_r [-inf, +inf]\n");
	printf("result = %g, abserr = %g, neval = %d, last = %d\n", result, abserr, neval, last);
	inf = 1;
	irev = 0;
	do {
		qagi_r(bound, inf, epsabs, epsrel, limit, &result, &abserr, &neval, &last, work, lwork, iwork, &info, &xx, yy, &irev);
		if (irev != 0)
			yy = f_qagi(xx);
	} while (irev != 0);
	printf("** qagi_r [0, +inf]\n");
	printf("result = %g, abserr = %g, neval = %d, last = %d\n", result, abserr, neval, last);
	inf = -1;
	irev = 0;
	do {
		qagi_r(bound, inf, epsabs, epsrel, limit, &result, &abserr, &neval, &last, work, lwork, iwork, &info, &xx, yy, &irev);
		if (irev != 0)
			yy = f_qagi(xx);
	} while (irev != 0);
	printf("** qagi_r [-inf, 0]\n");
	printf("result = %g, abserr = %g, neval = %d, last = %d\n", result, abserr, neval, last);
}

#undef N
#undef LWORK
#undef LIWORK

/* Parameters for test_derkfa */
#define N	4
#define LWORK	(11*N + 40)
#define LIWORK	40

void f_derkfa(int n, double t, double y[], double yp[])
{
	double alfa = M_PI_4, r;

	r = y[0]*y[0] + y[1]*y[1];
	r = r*sqrt(r)/(alfa*alfa);
	yp[0] = y[2];
	yp[1] = y[3];
	yp[2] = -y[0]/r;
	yp[3] = -y[1]/r;
}

void test_derkfa()
{
	int n = N;
	double t, y[N], tout, tend;
	double work[LWORK];
	int lwork = LWORK, iwork[LIWORK], liwork = LIWORK;
	int info;

	double ecc = 0.25, alfa = M_PI_4;
	double rtol = 1e-10, atol = rtol; int itol = 0;
	int mode = 2;

	printf("** derkfa\n");
	t = 0;
	y[0] = 1 - ecc;
	y[1] = 0;
	y[2] = 0;
	y[3] = alfa*sqrt((1 + ecc)/(1 - ecc));
	tend = 12;
	info = 0;
	for (int i = 1; i <= 12; i++) {
		tout = i;
		derkfa(n, f_derkfa, &t, y, tout, tend, &rtol, &atol, itol, mode, work, -lwork, iwork, -liwork, &info);
		if (info < 0 || info > 10)
			break;
		printf("%4.1f  %15.10f  %15.10f  %15.10f  %15.10f  %d  %d\n", t, y[0], y[1], y[2], y[3], iwork[13], info);
	}
	printf("info = %d\n", info);
}

#undef LWORK
#define LWORK	(9*N + 40)

void test_derkfa_r()
{
	int n = N;
	double t, y[N], tout, tend;
	double work[LWORK];
	int lwork = LWORK, iwork[LIWORK], liwork = LIWORK;
	double tt, yy[N], yyp[N]; int irev;
	int info;

	double ecc = 0.25, alfa = M_PI_4;
	double rtol = 1e-10, atol = rtol; int itol = 0;
	int mode = 2;

	printf("** derkfa_r\n");
	t = 0;
	y[0] = 1 - ecc;
	y[1] = 0;
	y[2] = 0;
	y[3] = alfa*sqrt((1 + ecc)/(1 - ecc));
	tend = 12;
	info = 0;
	for (int i = 1; i <= 12; i++) {
		tout = i;
		irev = 0;
		do {
			derkfa_r(n, &t, y, tout, tend, &rtol, &atol, itol, mode, work, -lwork, iwork, -liwork, &info, &tt, yy, yyp, &irev);
			if (irev != 0)
				f_derkfa(n, tt, yy, yyp);
		} while (irev != 0);
		if (info < 0 || info > 10)
			break;
		printf("%4.1f  %15.10f  %15.10f  %15.10f  %15.10f  %d  %d\n", t, y[0], y[1], y[2], y[3], iwork[13], info);
	}
	printf("info = %d\n", info);
}

#undef N
#undef LWORK
#undef LIWORK

/* Parameters for test_dopn43 */
#define N	2
#define LWORK	(10*N + 40)
#define LIWORK	30

void f_dopn43(int n, double t, double y[], double ypp[])
{
	double alfa = M_PI_4, r;

	r = y[0]*y[0] + y[1]*y[1];
	r = r*sqrt(r)/(alfa*alfa);
	ypp[0] = -y[0]/r;
	ypp[1] = -y[1]/r;
}

void test_dopn43()
{
	int n = N;
	double t, y[N], yp[N], tout, tend;
	double work[LWORK];
	int lwork = LWORK, iwork[LIWORK], liwork = LIWORK;
	int info;

	double ecc = 0.25, alfa = M_PI_4;
	double rtol = 1e-10, atol = rtol; int itol = 0;
	int mode = 2;

	printf("** dopn43\n");
	t = 0;
	y[0] = 1 - ecc;
	y[1] = 0;
	yp[0] = 0;
	yp[1] = alfa*sqrt((1 + ecc)/(1 - ecc));
	tend = 12;
	info = 0;
	for (int i = 1; i <= 12; i++) {
		tout = i;
		dopn43(n, f_dopn43, &t, y, yp, tout, tend, &rtol, &atol, itol, mode, work, -lwork, iwork, -liwork, &info);
		if (info < 0 || info > 10)
			break;
		printf("%4.1f  %15.10f  %15.10f  %15.10f  %15.10f  %d  %d\n", t, y[0], y[1], yp[0], yp[1], iwork[13], info);
	}
	printf("info = %d\n", info);
}

#undef LWORK
#define LWORK	(8*N + 40)

void test_dopn43_r()
{
	int n = N;
	double t, y[N], yp[N], tout, tend;
	double work[LWORK];
	int lwork = LWORK, iwork[LIWORK], liwork = LIWORK;
	double tt, yy[N], yypp[N]; int irev;
	int info;

	double ecc = 0.25, alfa = M_PI_4;
	double rtol = 1e-10, atol = rtol; int itol = 0;
	int mode = 2;

	printf("** dopn43_r\n");
	t = 0;
	y[0] = 1 - ecc;
	y[1] = 0;
	yp[0] = 0;
	yp[1] = alfa*sqrt((1 + ecc)/(1 - ecc));
	tend = 12;
	info = 0;
	for (int i = 1; i <= 12; i++) {
		tout = i;
		irev = 0;
		do {
			dopn43_r(n, &t, y, yp, tout, tend, &rtol, &atol, itol, mode, work, -lwork, iwork, -liwork, &info, &tt, yy, yypp, &irev);
			if (irev != 0)
				f_dopn43(n, tt, yy, yypp);
		} while (irev != 0);
		if (info < 0 || info > 10)
			break;
		printf("%4.1f  %15.10f  %15.10f  %15.10f  %15.10f  %d  %d\n", t, y[0], y[1], yp[0], yp[1], iwork[13], info);
	}
	printf("info = %d\n", info);
}

#undef N
#undef LWORK
#undef LIWORK

/* Parameters for test_rfft1 */
#define N	10
#define LWSAVE	17
#define INC	1
#define LR	(INC*(N - 1) + 1)
#define LWORK	N

void test_rfft1()
{
	double wsave[LWSAVE], r[LR], rcopy[LR], diff, work[LWORK];
	int n = N, lwsave = LWSAVE, lr = LR, inc = INC, lwork = LWORK, info, i, k;
	unsigned int seed;

	/* Initialization */
	printf("** rfft1\n");
	seed = 13;
	init_genrand(seed);
	rfft1i(n, wsave, lwsave, &info);
	if (info != 0) {
		printf("Error during initialization\n");
		return;
	}
	/* Generate test data */
	for (i = 0; i < n; i++) {
		k = inc*i;
		r[k] = genrand_res53();
		rcopy[k] = r[k];
	}
	/* Forward transform */
	rfft1f(n, inc, r, lr, wsave, lwsave, work, lwork, &info);
	if (info != 0) {
		printf("Error in rfft1f\n");
		return;
	}
	/* Backward transform */
	rfft1b(n, inc, r, lr, wsave, lwsave, work, lwork, &info);
	if (info != 0) {
		printf("Error in rfft1b\n");
		return;
	}
	/* Check results */
	diff = 0;
	for (i = 0; i < n; i++) {
		k = inc*i;
		if (fabs(r[k] - rcopy[k]) > diff)
			diff = fabs(r[k] - rcopy[k]);
		printf("%g  %g  %g\n", rcopy[k], r[k], fabs(r[k] - rcopy[k]));
	}
	printf("diff(max) = %g\n", diff);
}

#undef N
#undef LWSAVE
#undef INC
#undef LR
#undef LWORK

/* Parameters for test_lmdif1 */
#define M	14
#define N	2
#define LWORK	(N*(M + 6) + 2*M)
#define LIWORK	N

void f_lmdif1(int m, int n, double x[], double fvec[], int *iflag)
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
	double x[N], fvec[M], work[LWORK];
	double tol = 1.0e-10;
	int m = M, n = N, iwork[LIWORK], lwork = LWORK, info;

	x[0] = 500; x[1] = 0.0001;
	lmdif1(f_lmdif1, m, n, x, fvec, tol, work, lwork, iwork, &info);
	printf("** lmdif1\n");
	printf("x[1] = %g, x[2] = %g\n", x[0], x[1]);
	printf("info = %d\n", info);
}

void test_lmdif1_r()
{
	double x[N], fvec[M], work[LWORK];
	double tol = 1.0e-10;
	int m = M, n = N, iwork[LIWORK], lwork = LWORK, info;
	double xx[N], yy[M]; int iflag, irev;

	x[0] = 500; x[1] = 0.0001;
	irev = 0;
	do {
		lmdif1_r(m, n, x, fvec, tol, work, lwork, iwork, &info, xx, yy, &irev);
		if (irev >= 1 && irev <= 3) {
			iflag = 1;
			f_lmdif1(m, n, xx, yy, &iflag);
		}
	} while (irev != 0);
	printf("** lmdif1_r\n");
	printf("x[1] = %g, x[2] = %g\n", x[0], x[1]);
	printf("info = %d\n", info);
}

#undef M
#undef N
#undef LWORK
#undef LIWORK

void test_rand()
{
	unsigned int seed, r32;
	int r31;
	double r53;

	seed = 11;
	printf("** Random numbers: seed = %u\n", seed);
	init_genrand(seed);
	for (int i = 0; i < 10; i++) {
		r32 = genrand_int32();
		r31 = genrand_int31();
		r53 = genrand_res53();
		printf("%12u %12d %20.15g\n", r32, r31, r53);
	}
}

void test_dlamch()
{
	printf("** dlamch\n");
	printf("e: %g\n", dlamch('e'));
	printf("s: %g\n", dlamch('s'));
	printf("b: %g\n", dlamch('b'));
	printf("p: %g\n", dlamch('p'));
	printf("n: %g\n", dlamch('n'));
	printf("r: %g\n", dlamch('r'));
	printf("m: %g\n", dlamch('m'));
	printf("u: %g\n", dlamch('u'));
	printf("l: %g\n", dlamch('l'));
	printf("o: %g\n", dlamch('o'));
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
	test_derkfa();
	test_derkfa_r();
	test_dopn43();
	test_dopn43_r();
	test_rfft1();
	test_lmdif1();
	test_lmdif1_r();
	test_rand(); test_rand();
	test_dlamch();

	return 0;
}
