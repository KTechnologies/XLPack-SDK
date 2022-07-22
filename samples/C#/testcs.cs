/*****************************************
 *                                       *
 *  Experimental C# interface to XLPack  *
 *  Test program                         *
 *  Version 6.0 (June 14, 2022)          *
 *  (C) 2014-2022  K Technologies        *
 *                                       *
 *****************************************/

using System;
using System.IO;
using static System.Math;
using static XLPack;

public class Test
{

static void Test_sf()
{
	double x, x2, y, nu;
	uint ix;
	int errno;

	Console.Write("** d1num\n");
	Console.Write("i = 1: {0}\n", d1num(1));
	Console.Write("i = 2: {0}\n", d1num(2));
	Console.Write("i = 3: {0}\n", d1num(3));
	Console.Write("i = 4: {0}\n", d1num(4));

	Console.Write("** factorial\n");
	ix = 10;
	y = factorial(ix, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", ix, y, errno);

	Console.Write("** li\n");
	x = 2;
	y = li(x, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** ei\n");
	x = 1;
	y = ei(x, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** digamma\n");
	x = 3.5;
	y = digamma(x, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** besj0\n");
	x = 1;
	y = besj0(x, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** besj1\n");
	x = 1;
	y = besj1(x, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** besjnu\n");
	nu = 1; x = 1;
	y = besjnu(nu, x, out errno);
	Console.Write("nu = {0}, x = {1}, y = {2}, errno = {3}\n", nu, x, y, errno);

	Console.Write("** besy0\n");
	x = 1;
	y = besy0(x, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** besy1\n");
	x = 1;
	y = besy1(x, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** besynu\n");
	nu = 1; x = 1;
	y = besynu(nu, x, out errno);
	Console.Write("nu = {0}, x = {1}, y = {2}, errno = {3}\n", nu, x, y, errno);

	Console.Write("** besi0\n");
	x = 1;
	y = besi0(x, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);
	y = besi0(x, out errno, 2);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** besi1\n");
	x = 1;
	y = besi1(x, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);
	y = besi1(x, out errno, 2);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** besinu\n");
	nu = 1; x = 1;
	y = besinu(nu, x, out errno);
	Console.Write("nu = {0}, x = {1}, y = {2}, errno = {3}\n", nu, x, y, errno);

	Console.Write("** besk0\n");
	x = 1;
	y = besk0(x, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);
	y = besk0(x, out errno, 2);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** besk1\n");
	x = 1;
	y = besk1(x, out errno);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);
	y = besk1(x, out errno, 2);
	Console.Write("x = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** besknu\n");
	nu = 1; x = 1;
	y = besknu(nu, x, out errno);
	Console.Write("nu = {0}, x = {1}, y = {2}, errno = {3}\n", nu, x, y, errno);

	Console.Write("** celli1\n");
	x = 0.5;
	y = celli1(x, out errno);
	Console.Write("k = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** celli2\n");
	x = 0.5;
	y = celli2(x, out errno);
	Console.Write("k = {0}, y = {1}, errno = {2}\n", x, y, errno);

	Console.Write("** celli3\n");
	x = 0.7; x2 = 0.5;
	y = celli3(x, x2, out errno);
	Console.Write("n = {0}, k = {1}, y = {2}, errno = {3}\n", x, x2, y, errno);
}

static void Test_dconst()
{
	double x;

	Console.Write("** dconst\n");
	for (int i = 0; i <= 35; i++) {
		x = dconst(i);
		Console.Write("{0}: {1}\n", i, x);
	}
}

static void Test_dgesv()
{
	const int n = 3;
	double[,] a =
		{{ 0.2, -0.32, -0.8 },
		 { -0.11, 0.81, -0.92 },
		 { -0.93, 0.37, -0.29 }};
	double[] b = { -0.3727, 0.4319, -1.4247 };
	int[] ipiv = new int[n];
	double anorm, rcond = 0;
	int info;

	anorm = dlange('1', n, n, a);
	dgesv(n, a, ipiv, b, out info);
	if (info == 0)
		dgecon('1', n, a, anorm, out rcond, out info);
	Console.Write("** dgesv\n");
	Console.Write("x = {0}  {1}  {2}\n", b[0], b[1], b[2]);
	Console.Write("rcond = {0}, info = {1}\n", rcond, info);
}

static void Test_dposv()
{
	const int n = 3;
	double[,] a =
		{{ 2.2, 0.0, 0.0 },
		 { -0.11, 2.93, 0.0 },
		 { -0.32, 0.81, 2.37 }};
	double[] b = { -1.566, -2.8425, -1.1765 };
	double anorm, rcond = 0;
	int info;

	anorm = dlansy('1', 'U', n, a);
	dposv('U', n, a, b, out info);
	if (info == 0)
		dpocon('U', n, a, anorm, out rcond, out info);
	Console.Write("** dposv\n");
	Console.Write("x =\n");
	Console.Write("{0}  {1}  {2}\n", b[0], b[1], b[2]);
	Console.Write("rcond = {0}, info = {1}\n", rcond, info);
}

static void Test_dsyev()
{
	const int n = 3;
	double[,] a =
		{{ 2.2, 0.0, 0.0 },
		 { -0.11, 2.93, 0.0 },
		 { -0.32, 0.81, 2.37 }};
	double[] w = new double[n];
	int info;

	dsyev('V', 'U', n, a, w, out info);
	Console.Write("** dsyev\n");
	Console.Write("Eigenvalues =\n");
	Console.Write("{0}  {1}  {2}\n", w[0], w[1], w[2]);
	Console.Write("Eigenvectors =\n");
	Console.Write("{0,22}{1,22}{2,22}\n", a[0, 0], a[1, 0], a[2, 0]);
	Console.Write("{0,22}{1,22}{2,22}\n", a[0, 1], a[1, 1], a[2, 1]);
	Console.Write("{0,22}{1,22}{2,22}\n", a[0, 2], a[1, 2], a[2, 2]);
    Console.Write("info = {0}\n", info);
}

static void Test_dgels()
{
	const int m = 4, n = 2;
	double[] x = { 0.2, 118.2, 337.4, 884.6 };
	double[] y = { 0.1, 118.1, 338.8, 888.0 };
	double[,] a = new double[n, m];
	double[] ci = new double[n];
	double s;
	int info;

	for (int i = 0; i < m; i++) {
		a[0, i] = 1;
		a[1, i] = x[i];
	}
	dgels('n', m, n, a, y, out info);
	Console.Write("** dgels\n");
	Console.Write("a0 = {0}, a1 = {1}\n", y[0], y[1]);
	Console.Write("info = {0}\n", info);
	if (info == 0) {
		dgecov(0, n, a, ci, out info);
		s = 0.0;
		for (int i = n; i < m; i++)
			s = s + y[i]*y[i];
		s = s / (m - n);
		Console.Write("Std. dev. = {0}, {1}\n", Sqrt(s*ci[0]), Sqrt(s*ci[1]));
		Console.Write("info = {0}\n", info);
	}
}

static void Test_pchse()
{
	const int n = 4, ne = 2;
	double[] x = { 0.1, 0.11, 0.12, 0.13 };
	double[] y = { 2.3026, 2.2073, 2.1203, 2.0402 };
	double[] d = new double[n], xe = new double[2], ye = new double[2];
	int info;

	pchse(n, x, y, d, out info);
	Console.Write("** pchse\n");
	Console.Write("info = {0}\n", info);
	if (info == 0) {
		xe[0] = 0.115; xe[1] = 0.125;
		pchfe(n, x, y, d, ne, xe, ye, out info);
		Console.Write("ln({0}) = {1}\n", xe[0], ye[0]);
		Console.Write("ln({0}) = {1}\n", xe[1], ye[1]);
		Console.Write("info = {0}\n", info);
	}
}

static void Test_pchia()
{
	const int n = 7;
	double[] x = { -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
	double[] y = { 0.5, 1.0, 0.5, 0.2, 0.1, 0.05882, 0.03846 };
	double[] d = new double[n];
	double a, b, s;
	int info;

	pchse(n, x, y, d, out info);
	Console.Write("** pchse\n");
	Console.Write("info = {0}\n", info);
	if (info == 0) {
		a = 0; b = 4;
		s = pchia(n, x, y, d, a, b, out info);
		Console.Write("** pchia\n");
		Console.Write("S = {0}, S(true) = {1}\n", s, Atan(4.0));
		Console.Write("info = {0}\n", info);
	}
}

static void Test_rpzero2()
{
	const int n = 5;
	double[] a = { 1.0, 0.0, 2.0, 2.0, -15.0, 10.0 };
	double[] zr = new double[n], zi = new double[n], s = new double[n];
	int iter, info;

	rpzero2(n, a, zr, zi, out iter, s, out info);
	Console.Write("** rpzero2\n");
	for (int i = 0; i < n; i++)
		Console.Write("{0,22}{1,22}{2,22}\n", zr[i], zi[i], s[i]);
	Console.Write("iter = {0}, info = {1}\n", iter, info);
}

static double f_dfzero(double x)
{
	return (x*x - 2)*x - 5;
}

static void Test_dfzero()
{
	double b, c, r;
	int info;

	b = 1.0; c = 3.0; r = b;
	dfzero(f_dfzero, out b, out c, r, out info);
	Console.Write("** dfzero\n");
	Console.Write("x = {0}, info = {1}\n", b, info);
}

static void f_hybrd1(int n, double[] x, double[] fvec, ref int iflag)
{
	fvec[0] = 4*x[0]*x[0] + x[1]*x[1] - 16;
	fvec[1] = x[0]*x[0] + x[1]*x[1] - 9;
}

static void Test_hybrd1()
{
	const int n = 2;
	double[] x = new double[n], fvec = new double[n];
	int info;

	x[0] = 1; x[1] = 2;
	hybrd1(f_hybrd1, n, x, fvec, out info);
	Console.Write("** hybrd1\n");
	Console.Write("x[1] = {0}, x[2] = {1}\n", x[0], x[1]);
	Console.Write("info = {0}\n", info);
}

static double f_dfmin(double x)
{
	return (x*x - 2)*x - 5;
}

static void Test_dfmin()
{
	double a, b, x;

	a = 0; b = 1;
	x = dfmin(a, b, f_dfmin);
	Console.Write("** dfmin\n");
	Console.Write("x = {0}\n", x);
}

static void f_optif0(int n, double[] x, out double fval)
{
	// Rosenbrock function
	fval = 100*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]) + (1 - x[0])*(1 - x[0]);
}

static void Test_optif0()
{
	const int n = 2;
	double[] x = new double[n], xpls = new double[n];
	double fpls;
	int info;

	x[0] = -1.2; x[1] = 1.0;
	optif0(n, x, f_optif0, xpls, out fpls, out info);
	Console.Write("** optif0\n");
	Console.Write("xpls[1] = {0}, xpls[2] = {1}\n", xpls[0], xpls[1]);
	Console.Write("fpls = {0}\n", fpls);
	Console.Write("info = {0}\n", info);
}

static double f_qk15(double x)
{
	return 1/(1 + x*x);
}

static void Test_qk15()
{
	double a, b, result, abserr;

	// int(1/(1+x**2)) [0, 4] = atan(4)
	a = 0; b = 4;
	qk15(f_qk15, a, b, out result, out abserr);
	Console.Write("** qk15\n");
	Console.Write("result = {0}, abserr = {1}\n", result, abserr);
}

static double f_qag(double x)
{
	return 1/(1 + x*x);
}

static void Test_qag()
{
	double a, b, epsabs, epsrel, result, abserr;
	int key, neval, last, info;

	// int(1/(1+x**2)) [0, 4] = atan(4)
	a = 0; b = 4;
	epsabs = 1.0e-10; epsrel = 1.0e-10;
	for (int i = 1; i <= 6; i++) {
		key = i;
		qag(f_qag, a, b, out result, out abserr, out neval, out last, out info, epsabs, epsrel, key);
		Console.Write("** qag (key = {0})\n", key);
		Console.Write("result = {0}, abserr = {1}, neval = {2}, last = {3}\n", result, abserr, neval, last);
	}
}

static double f_qagi(double x)
{
	return 1/(1 + x*x);
}

static void Test_qagi()
{
	double bound, result, abserr;
	int inf, neval, last, info;

	// int(1/(1+x**2)) [-inf, +inf] = pi
	// int(1/(1+x**2)) [0, +inf] = pi / 2
	// int(1/(1+x**2)) [-inf, 0] = pi / 2
	bound = 0;
	inf = 2;
	qagi(f_qagi, bound, inf, out result, out abserr, out neval, out last, out info);
	Console.Write("** qagi [-inf, +inf]\n");
	Console.Write("result = {0}, abserr = {1}, neval = {2}, last = {3}\n", result, abserr, neval, last);
	inf = 1;
	qagi(f_qagi, bound, inf, out result, out abserr, out neval, out last, out info);
	Console.Write("** qagi [0, +inf]\n");
	Console.Write("result = {0}, abserr = {1}, neval = {2}, last = {3}\n", result, abserr, neval, last);
	inf = -1;
	qagi(f_qagi, bound, inf, out result, out abserr, out neval, out last, out info);
	Console.Write("** qagi [-inf, 0]\n");
	Console.Write("result = {0}, abserr = {1}, neval = {2}, last = {3}\n", result, abserr, neval, last);
}

static double alfasq;

static void f_derkf(int n, double t, double[] y, double[] yp)
{
	double r;

	r = y[0]*y[0] + y[1]*y[1];
	r = r*Sqrt(r)/alfasq;
	yp[0] = y[2];
	yp[1] = y[3];
	yp[2] = -y[0]/r;
	yp[3] = -y[1]/r;
}

static void Test_derkf()
{
	const int n = 4;
	const int lwork = 7*n + 25, liwork = 20;
	double t, tout;
	double[] y = new double[n];
	double[] work = new double[lwork];
	int[] iwork = new int[liwork];
	double tfinal, tprint;
	int info;
	double ecc = 0.25, alfa = Math.PI/4;

	alfasq = alfa*alfa;
	t = 0;
	y[0] = 1 - ecc;
	y[1] = 0;
	y[2] = 0;
	y[3] = alfa*Sqrt((1 + ecc)/(1 - ecc));
	tfinal = 12;
	tprint = 1;
	Console.Write("** derkf\n");
	info = 0;
	do {
		tout = t + tprint;
		derkf(n, f_derkf, ref t, y, tout, work, iwork, out info);
		if (info != 1)
			break;
		Console.Write("t = {0}, y[1] = {1}, y[2] = {2}, y[3] = {3}, y[4] = {4}\n", t, y[0], y[1], y[2], y[3]);
	} while (t < tfinal);
	Console.Write("info = {0}\n", info);
}

static void Test_derkf_2()
{
	const int n = 4;
	const int lwork = 9*n + 25, liwork = 20;
	double t, tout;
	double[] y = new double[n];
	double[] work = new double[lwork];
	int[] iwork = new int[liwork];
	double tfinal, tprint;
	int info;
	double ecc = 0.25, alfa = Math.PI/4;

	alfasq = alfa*alfa;
	t = 0;
	y[0] = 1 - ecc;
	y[1] = 0;
	y[2] = 0;
	y[3] = alfa*Sqrt((1 + ecc)/(1 - ecc));
	tfinal = 12;
	tprint = 1;
	Console.Write("** derkf (2) (dense output)\n");
	tout = t + tprint;
	info = 0;
	do {
		derkf(n, f_derkf, ref t, y, tfinal, work, iwork, out info, mode:2);
		if (info == 1 || info == 2) {
			while (t >= tout) {
				double[] y1 = new double[n];
				derkf_int(n, tout, y1, work);
				Console.Write("t = {0}, y[1] = {1}, y[2] = {2}, y[3] = {3}, y[4] = {4}\n", tout, y1[0], y1[1], y1[2], y1[3]);
				tout = tout + tprint;
			}
		} else
			break;
	} while (t < tfinal);
	Console.Write("info = {0}\n", info);
}

static void Test_rfft1()
{
	double[] wsave, r, rcopy;
	double diff;
	int n, lwsave, info, i;
	uint seed;

	// Initialization
	Console.Write("** rfft1\n");
	seed = 13;
	init_genrand(seed);
	n = 10;
	lwsave = n + (int)(Log(n)/Log(2)) + 4;
	wsave = new double[lwsave];
	rfft1i(n, wsave, out info);
	if (info != 0) {
		Console.Write("Error during initialization\n");
		return;
	}
	// Generate test data
	r = new double[n];
	rcopy = new double[n];
	for (i = 0; i < n; i++) {
		r[i] = genrand_res53();
		rcopy[i] = r[i];
	}
	// Forward transform
	rfft1f(n, r, wsave, out info);
	if (info != 0) {
		Console.Write("Error in rfft1f\n");
		return;
	}
	// Backward transform
	rfft1b(n, r, wsave, out info);
	if (info != 0) {
		Console.Write("Error in rfft1b\n");
		return;
	}
	// Check results
	diff = 0;
	for (i = 0; i < n; i++) {
		if (Abs(r[i] - rcopy[i]) > diff)
			diff = Abs(r[i] - rcopy[i]);
		Console.Write("{0,22}{1,22}{2,22}\n", rcopy[i], r[i], Abs(r[i] - rcopy[i]));
	}
	Console.Write("diff(max) = {0}\n", diff);
}

static void f_lmdif1(int m, int n, double[] x, double[] fvec, ref int iflag)
{
	double[] xdata = { 77.6, 114.9, 141.1, 190.8, 239.9, 289.0, 332.8, 378.4,
		434.8, 477.3, 536.8, 593.1, 689.1, 760.0 };
	double[] ydata = { 10.07, 14.73, 17.94, 23.93, 29.61, 35.18, 40.02, 44.82,
		50.76, 55.05, 61.01,  66.4, 75.47, 81.78 };

	for (int i = 0; i < m; i++)
		fvec[i] = ydata[i] - x[0]*(1 - Exp(-xdata[i]*x[1]));
}

static void Test_lmdif1()
{
	const int m = 14, n = 2;
	double[] x = new double[n], fvec = new double[m];
	int info;

	x[0] = 500; x[1] = 0.0001;
	lmdif1(f_lmdif1, m, n, x, fvec, out info);
	Console.Write("** lmdif1\n");
	Console.Write("x[1] = {0}, x[2] = {1}\n", x[0], x[1]);
	Console.Write("info = {0}\n", info);
}

static void Test_rand()
{
	uint seed, r32;
	int r31;
	double r53;

	seed = 11;
	Console.Write("** Random numbers: seed = {0}\n", seed);
	init_genrand(seed);
	for (int i = 0; i < 10; i++) {
		r32 = genrand_int32();
		r31 = genrand_int31();
		r53 = genrand_res53();
		Console.Write("{0,12}{1,12}{2,22}\n", r32, r31, r53);
	}
}

static void Test_dlamch()
{
	Console.Write("** dlamch\n");
	Console.Write("e: {0}\n", dlamch('e'));
	Console.Write("s: {0}\n", dlamch('s'));
	Console.Write("b: {0}\n", dlamch('b'));
	Console.Write("p: {0}\n", dlamch('p'));
	Console.Write("n: {0}\n", dlamch('n'));
	Console.Write("r: {0}\n", dlamch('r'));
	Console.Write("m: {0}\n", dlamch('m'));
	Console.Write("u: {0}\n", dlamch('u'));
	Console.Write("l: {0}\n", dlamch('l'));
	Console.Write("o: {0}\n", dlamch('o'));
}

public static int Main(string[] args)
{

	Test_sf();
	Test_dconst();
	Test_dgesv();
	Test_dposv();
	Test_dsyev();
	Test_dgels();
	Test_pchse();
	Test_pchia();
	Test_rpzero2();
	Test_dfzero();
	Test_hybrd1();
	Test_dfmin();
	Test_optif0();
	Test_qk15();
	Test_qag();
	Test_qagi();
	Test_derkf();
	Test_derkf_2();
	Test_rfft1();
	Test_lmdif1();
	Test_rand(); Test_rand();
	Test_dlamch();
	return 0;
}

}
