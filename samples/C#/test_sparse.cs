/*****************************************
 *                                       *
 *  Experimental C# interface to XLPack  *
 *  Test program                         *
 *  Version 7.0 (January 31, 2023)       *
 *  (C) 2014-2023  K Technologies        *
 *                                       *
 *****************************************/

using System;
using System.IO;
using static System.Math;
using static XLPack;

public class Test
{

// Test SPUtils

static void Test_csr_dense()
{
	double[] a = {
		0.2, 0.11, 0, -0.11, -0.8,
		0.11, 0.93, 0.81, 0.93, -0.92,
		0, -0.81, 0.37, 0.81, 0,
		0.8, 0.92, 0, 0.9, 0.86 };
	double[] a_val = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int[] a_ptr = { 0, 3, 7, 9, 13, 16 };
	int[] a_ind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	const int m = 5, n = 4, lda = 5, index = 0;
	double[] a2 = new double[m*n];
	int info;

	Console.WriteLine("** Test csr_dense");
	csr_dense(m, n, a_val, a_ptr, a_ind, index, lda, a2, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				if (a[j*lda + i] != a2[j*lda + i])
					Console.WriteLine("error: unmatched {0} {1}", i, j);
			}
		}
	}
}

static void Test_dense_csr()
{
	double[] a = {
		0.2, 0.11, 0, -0.11, -0.8,
		0.11, 0.93, 0.81, 0.93, -0.92,
		0, -0.81, 0.37, 0.81, 0,
		0.8, 0.92, 0, 0.9, 0.86 };
	double[] a_val = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int[] a_ptr = { 0, 3, 7, 9, 13, 16 };
	int[] a_ind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	const int m = 5, n = 4, lda = 5, nnz = 16, index = 0;
	double[] a2_val = new double[nnz];
	int[] a2_ptr = new int[m + 1], a2_ind = new int[nnz];
	int info;

	Console.WriteLine("** Test dense_csr");
	dense_csr(m, n, lda, a, nnz, a2_val, a2_ptr, a2_ind, index, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		for (int i = 0; i < m + 1; i++) {
			if (a2_ptr[i] != a_ptr[i])
				Console.WriteLine("error: unmatched ptr {0}", i);
		}
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != a_val[i])
				Console.WriteLine("error: unmatched val {0}", i);
			if (a2_ind[i] != a_ind[i])
				Console.WriteLine("error: unmatched ind {0}", i);
		}
	}
}

static void Test_csr_coo()
{
	double[] a_val = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int[] a_rowptr = { 0, 3, 7, 9, 13, 16 };
	int[] a_colind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	int[] a_rowind = { 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4 };
	int m = 5, n = 4, nnz = 16, index = 0;
	double[] a2_val = new double[nnz];
	int[] a2_rowind = new int[nnz], a2_colind = new int[nnz];
	int info;

	Console.WriteLine("** Test csr_coo");
	csr_coo(m, n, a_val, a_rowptr, a_colind, index, a2_val, a2_rowind, a2_colind, index, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != a_val[i])
				Console.WriteLine("error: unmatched val {0}", i);
			if (a2_rowind[i] != a_rowind[i])
				Console.WriteLine("error: unmatched rowind {0}", i);
			if (a2_colind[i] != a_colind[i])
				Console.WriteLine("error: unmatched colind {0}", i);
		}
	}
}

static void Test_coo_csr()
{
	double[] a_val = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int[] a_rowptr = { 0, 3, 7, 9, 13, 16 };
	int[] a_colind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	int[] a_rowind = { 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4 };
	int m = 5, n = 4, nnz = 16, index = 0;
	double[] a2_val = new double[nnz];
	int[] a2_rowptr = new int[m + 1];
	int[] a2_colind = new int[nnz];
	int info;

	Console.WriteLine("** Test coo_csr");
	coo_csr(m, n, nnz, a_val, a_rowind, a_colind, index, a2_val, a2_rowptr, a2_colind, index, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		for (int i = 0; i < m + 1; i++) {
			if (a2_rowptr[i] != a_rowptr[i])
				Console.WriteLine("error: unmatched rowptr {0}", i);
		}
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != a_val[i])
				Console.WriteLine("error: unmatched val {0}", i);
			if (a2_colind[i] != a_colind[i])
				Console.WriteLine("error: unmatched colind {0}", i);
		}
	}
}

static void Test_csr_ssr()
{
	double[] a_val = { 2.2, -0.11, -0.8, -0.11, 2.93, 0.81, -0.92, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int[] a_ptr = { 0, 3, 7, 9, 12 };
	int[] a_ind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 3 };
	double[] s_val = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int[] s_ptr = { 0, 1, 3, 5, 8 };
	int[] s_ind = { 0, 0, 1, 1, 2, 0, 1, 3 };
	int n = 4, nnz = 8, index = 0;
	double[] a2_val = new double[nnz];
	int[] a2_ptr = new int[n + 1];
	int[] a2_ind = new int[nnz];
	int info;

	Console.WriteLine("** Test csr_ssr");
	csr_ssr('L', n, a_val, a_ptr, a_ind, index, nnz, a2_val, a2_ptr, a2_ind, index, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		for (int i = 0; i < n + 1; i++) {
			if (a2_ptr[i] != s_ptr[i])
				Console.WriteLine("error: unmatched ptr {0}", i);
		}
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != s_val[i])
				Console.WriteLine("error: unmatched val {0}", i);
			if (a2_ind[i] != s_ind[i])
				Console.WriteLine("error: unmatched ind {0}", i);
		}
	}
}

static void Test_ssr_csr()
{
	double[] a_val = { 2.2, -0.11, -0.8, -0.11, 2.93, 0.81, -0.92, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int[] a_ptr = { 0, 3, 7, 9, 12 };
	int[] a_ind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 3 };
	double[] s_val = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int[] s_ptr = { 0, 1, 3, 5, 8 };
	int[] s_ind = { 0, 0, 1, 1, 2, 0, 1, 3 };
	int n = 4, nnz = 12, index = 0;
	double[] a2_val = new double[nnz];
	int[] a2_ptr = new int[n + 1];
	int[] a2_ind = new int[nnz];
	int info;

	Console.WriteLine("** Test ssr_csr");
	ssr_csr('L', n, s_val, s_ptr, s_ind, index, nnz, a2_val, a2_ptr, a2_ind, index, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		for (int i = 0; i < n + 1; i++) {
			if (a2_ptr[i] != a_ptr[i])
				Console.WriteLine("error: unmatched ptr {0}", i);
		}
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != a_val[i])
				Console.WriteLine("error: unmatched val {0}", i);
			if (a2_ind[i] != a_ind[i])
				Console.WriteLine("error: unmatched ind {0}", i);
		}
	}
}

static void Test_csc_csr()
{
	double[] a_val = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int[] a_ptr = { 0, 3, 7, 9, 13, 16 };
	int[] a_ind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	double[] at_val = { 0.2, 0.11, -0.11, -0.8, 0.11, 0.93, 0.81, 0.93, -0.92, -0.81, 0.37, 0.81, 0.8, 0.92, 0.9, 0.86 };
	int[] at_ptr = { 0, 4, 9, 12, 16 };
	int[] at_ind = { 0, 1, 3, 4, 0, 1, 2, 3, 4, 1, 2, 3, 0, 1, 3, 4 };
	int m = 5, n = 4, nnz = 16, index = 0;
	double[] a2_val = new double[nnz];
	int[] a2_ptr = new int[n + 1];
	int[] a2_ind = new int[nnz];
	int info;

	Console.WriteLine("** Test csc_csr");
	csc_csr(n, m, a_val, a_ptr, a_ind, index, a2_val, a2_ptr, a2_ind, index, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		for (int i = 0; i < n + 1; i++) {
			if (a2_ptr[i] != at_ptr[i])
				Console.WriteLine("error: unmatched ptr {0}", i);
		}
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != at_val[i])
				Console.WriteLine("error: unmatched val {0}", i);
			if (a2_ind[i] != at_ind[i])
				Console.WriteLine("error: unmatched ind {0}", i);
		}
	}
}

// Test SPBLAS

static void Test_csr_dusmv()
{
	double[] a = {
		0.2, 0.11, 0, -0.11, -0.8,
		0.11, 0.93, 0.81, 0.93, -0.92,
		0, -0.81, 0.37, 0.81, 0,
		0.8, 0.92, 0, 0.9, 0.86 };
	double[] a_val = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int[] a_ptr = { 0, 3, 7, 9, 13, 16 };
	int[] a_ind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	double[] x = { 0.73, 0.86, 0.44, 0.87 };
	double[] y = { 0.28, 0.30, 0.69, 0.50, 0.15 };
	double alpha = 2.0, beta = 3.0;
	int m = 5, n = 4, lda = 5, index = 0;
	double[] y_bup = new double[m];
	int info;

	Console.WriteLine("** Test csr_dusmv");
	for (int i = 0; i < m; i++)
		y_bup[i] = y[i];
	csr_dusmv('N', m, n, alpha, a_val, a_ptr, a_ind, index, x, 1, beta, y, 1, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		double err = 0;
		for (int i = 0; i < m; i++) {
			double s = 0;
			for (int j = 0; j < n; j++)
				s = s + a[lda*j + i]*x[j];
			err = Max(err, Abs(y[i] - (alpha*s + beta*y_bup[i])));
		}
		Console.WriteLine("err = {0}", err);
	}
}

static void Test_csr_dussv()
{
	double[] a = {
		2.2, -0.11, 0, -0.8,
		0, 2.93, 0.81, -0.92,
		0, 0, 2.37, 0,
		0, 0, 0, 2.86 };
	double[] s_val = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86 };
	int[] s_ptr = { 0, 1, 3, 5, 8 };
	int[] s_ind = { 0, 0, 1, 1, 2, 0, 1, 3 };
	double[] b = { 0.73, 0.86, 0.44, 0.87 };
	int n = 4, lda = 4, index = 0;
	double[] x = new double[n];
	int info;

	Console.WriteLine("** Test csr_dussv");
	for (int i = 0; i < n; i++)
		x[i] = b[i];
	csr_dussv('L', 'N', 'N', n, s_val, s_ptr, s_ind, index, x, 1, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		double err = 0;
		for (int i = 0; i < n; i++) {
			double s = 0;
			for (int j = 0; j <= i; j++)
				s = s + a[lda*j + i]*x[j];
			err = Max(err, Abs(b[i] - s));
		}
		Console.WriteLine("err = {0}", err);
	}
}

static void Test_csr_dusmm()
{
	double[] a = {
		0.2, 0.11, 0, -0.11, -0.8,
		0.11, 0.93, 0.81, 0.93, -0.92,
		0, -0.81, 0.37, 0.81, 0,
		0.8, 0.92, 0, 0.9, 0.86 };
	double[] a_val = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int[] a_ptr = { 0, 3, 7, 9, 13, 16 };
	int[] a_ind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	double[] b = { 0.73, 0.86, 0.44, 0.87,
		0.0057, 0.89, 0.75, 0.18,
		0.71, 0.46, 0.69, 0.64 };
	double[] c = { 0.28, 0.30, 0.69, 0.50, 0.15,
		0.99, 0.24, 0.4, 0.11, 0.48,
		0.0029, 0.67, 0.38, 0.46, 0.014 };
	double alpha = 2, beta = 3;
	int m = 5, n = 4, l = 3, lda = 5, ldb = 4, ldc = 5, index = 0;
	double[] c_bup = new double[m*l];
	int info;

	Console.WriteLine("** Test csr_dusmm");
	for (int i = 0; i < m*l; i++)
		c_bup[i] = c[i];
	csr_dusmm('N', 'C', m, n, l, alpha, a_val, a_ptr, a_ind, index, ldb, b, beta, ldc, c, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		double err = 0;
		for (int k = 0; k < l; k++) {
			for (int i = 0; i < m; i++) {
				double s = 0;
				for (int j = 0; j < n; j++)
					s = s + a[lda*j + i]*b[ldb*k + j];
				err = Max(err, Abs(c[ldc*k + i] - (alpha*s + beta*c_bup[ldc*k + i])));
			}
		}
		Console.WriteLine("err = {0}", err);
	}
}

static void Test_csr_dussm()
{
	double[] a = {
		2.2, -0.11, 0, -0.8,
		0, 2.93, 0.81, -0.92,
		0, 0, 2.37, 0,
		0, 0, 0, 2.86 };
	double[] s_val = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86 };
	int[] s_ptr = { 0, 1, 3, 5, 8 };
	int[] s_ind = { 0, 0, 1, 1, 2, 0, 1, 3 };
	double[] b = { 0.73, 0.86, 0.44, 0.87,
		0.0057, 0.89, 0.75, 0.18,
		0.71, 0.46, 0.69, 0.64 };
	int n = 4, nrhs = 3, lda = 4, ldb = 4, ldx = 4, index = 0;
	double[] x = new double[n*nrhs];
	int info;

	Console.WriteLine("** Test csr_dussm");
	for (int i = 0; i < n*nrhs; i++)
		x[i] = b[i];
	csr_dussm('L', 'N', 'N', 'C', n, nrhs, s_val, s_ptr, s_ind, index, ldx, x, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		double err = 0;
		for (int k = 0; k < nrhs; k++) {
			for (int i = 0; i < n; i++) {
				double s = 0;
				for (int j = 0; j <= i; j++)
					s = s + a[lda*j + i]*x[ldx*k + j];
				err = Max(err, Abs(b[ldb*k + i] - s));
			}
		}
		Console.WriteLine("err = {0}", err);
	}
}

static void Test_ssr_dusmv()
{
	double[] a = {
		2.2, -0.11, 0, -0.8,
		-0.11, 2.93, 0.81, -0.92,
		0, 0.81, 2.37, 0,
		-0.8, -0.92, 0, 2.86 };
	double[] s_val = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86 };
	int[] s_ptr = { 0, 1, 3, 5, 8 };
	int[] s_ind = { 0, 0, 1, 1, 2, 0, 1, 3 };
	double[] x = { 0.73, 0.86, 0.44, 0.87 };
	double[] y = { 0.28, 0.30, 0.69, 0.50 };
	double alpha = 2.0, beta = 3.0;
	int n = 4, lda = 4, index = 0;
	double[] y_bup = new double[n];
	int info;

	Console.WriteLine("** Test ssr_dusmv");
	for (int i = 0; i < n; i++)
		y_bup[i] = y[i];
	ssr_dusmv('L', n, alpha, s_val, s_ptr, s_ind, index, x, 1, beta, y, 1, out info);
	Console.WriteLine("info = {0}", info);
	if (info == 0) {
		double err = 0;
		for (int i = 0; i < n; i++) {
			double s = 0;
			for (int j = 0; j < n; j++)
				s = s + a[lda*j + i]*x[j];
			err = Max(err, Abs(y[i] - (alpha*s + beta*y_bup[i])));
		}
		Console.WriteLine("err = {0}", err);
	}
}

// Test iterative solver

static void Test_cg1()
{
	double[] s_val = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int[] s_ptr = { 0, 1, 3, 5, 8 };
	int[] s_ind = { 0, 0, 1, 1, 2, 0, 1, 3 };
	double[] b = { 1.433, 1.3137, 2.3799, -0.7992 };
	double[] x = { 0, 0, 0, 0 };
	double tol = 1.0e-10, res;
	int maxiter = 100, iter;
	int n = 4, info;

	Console.WriteLine("** Test_cg1");
	cg1('L', n, s_val, s_ptr, s_ind, b, x, tol, maxiter, out iter, out res, out info);
	Console.WriteLine("x =");
	Console.WriteLine("{0}  {1}  {2}  {3}", x[0], x[1], x[2], x[3]);
	Console.WriteLine("iter = {0}, res = {1}, info = {2}", iter, res, info);
}

static void Test_bicg1()
{
	double[] a_val = { 2.2, -0.11, -0.8, 0.11, 2.93, 0.81, -0.92, 0.81, 2.37, 0.8, 0.92, 0.86 };
	int[] a_ptr = { 0, 3, 7, 9, 12 };
	int[] a_ind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 3 };
	double[] b = { 1.433, 1.4325, 2.3799, 0.2488 };
	double[] x = { 0, 0, 0, 0 };
	double tol = 1.0e-10, res;
	int maxiter = 100, iter;
	int n = 4, info;

	Console.WriteLine("** Test_bicg1");
	bicg1(n, a_val, a_ptr, a_ind, b, x, tol, maxiter, out iter, out res, out info);
	Console.WriteLine("x =");
	Console.WriteLine("{0}  {1}  {2}  {3}", x[0], x[1], x[2], x[3]);
	Console.WriteLine("iter = {0}, res = {1}, info = {2}", iter, res, info);
}

// Test I/O routines

static void Test_mmio()
{
	string fname = "test.mtx";
	string matcode;
	double[] a_val = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int[] a_ptr = { 0, 3, 7, 9, 13, 16 };
	int[] a_ind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	int m = 5, n = 4, nnz = 16;
	double[] b_val = new double[nnz];
	int[] b_ptr = new int[m + 1];
	int[] b_ind = new int[nnz];
	int skip = 0, index = 0, fchk = 0, format = 0, sort = 0;
	int nrow, ncol, nzero, info;

	Console.WriteLine("** Test mmio");
	// Write t_csr to the file
	mm_write(fname, "MCRS", m, n, nnz, a_val, nnz, a_ptr, m + 1, a_ind, nnz, index, fchk, format, out info);
	Console.WriteLine("mm_write: info = {0}", info);
	if (info != 0)
		return;
	// Read to a_csr
	mm_read(fname, out matcode, out nrow, out ncol, out nzero, b_val, nnz, b_ptr, m + 1, b_ind, nnz, skip, index, format, sort, out info);
	Console.WriteLine("mm_read: info = {0}", info);
	if (info != 0)
		return;
	Console.WriteLine("m = {0}, n = {1}, nnz = {2}, matcode = {3}", nrow, ncol, nzero, matcode);
	// Compare data
	for (int i = 0; i < m + 1; i++) {
		if (b_ptr[i] != a_ptr[i])
			Console.WriteLine("error: unmatched ptr {0}", i);
	}
	for (int i = 0; i < nnz; i++) {
		if (b_val[i] != a_val[i])
			Console.WriteLine("error: unmatched val {0}", i);
		if (b_ind[i] != a_ind[i])
			Console.WriteLine("error: unmatched ind {0}", i);
	}
}

// Test chexk routines

static void Test_csr_check()
{
	double[] a_val = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92, 0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int[] a_ptr = { 0, 3, 7, 9, 13, 16 };
	int[] a_ind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	int[] result = new int[10];
	int n = 4, info;

	Console.WriteLine("** Test csr_check");
	csr_check(n, n, a_val, a_ptr, a_ind, result, out info);
	Console.WriteLine("csr_check: info = {0}", info);
	Console.WriteLine("  base = {0}", result[0]);
	Console.WriteLine("  nnz = {0}", result[1]);
	Console.WriteLine("  nnz_lower = {0}", result[2]);
	Console.WriteLine("  nnz_upper = {0}", result[3]);
	Console.WriteLine("  nnz_diag = {0}", result[4]);
	Console.WriteLine("  n_zero = {0}", result[5]);
	Console.WriteLine("  n_zero_diag = {0}", result[6]);
	Console.WriteLine("  empty_rows = {0}", result[7]);
	Console.WriteLine("  unsorted_rows = {0}", result[8]);
	Console.WriteLine("  invalid_inds = {0}", result[9]);
}

static void Test_csr_check_sym()
{
	double[] a_val = { 2.2, -0.11, -0.8, -0.11, 2.93, 0.81, -0.92, 0.81, 2.37, -0.8, -0.92, 0.86 };
	double[] a2_val = { 2.2, 0.11, -0.8, -0.11, 2.93, 0.81, -0.92, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int[] a_ptr = { 0, 3, 7, 9, 12 };
	int[] a_ind = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 3 };
	int[] result = new int[10];
	int n = 4, info;

	Console.WriteLine("** Test csr_check_sym");
	csx_check_sym(n, a_val, a_ptr, a_ind, out info);
	Console.WriteLine("csx_check_sym: info = {0}", info);
	csx_check_sym(n, a2_val, a_ptr, a_ind, out info);
	Console.WriteLine("csx_check_sym: info = {0}", info);
}

// Main program

public static int Main(string[] args)
{
	Test_csr_dense();
	Test_dense_csr();
	Test_csr_coo();
	Test_coo_csr();
	Test_csr_ssr();
	Test_ssr_csr();
	Test_csc_csr();

	Test_csr_dusmv();
	Test_csr_dussv();
	Test_csr_dusmm();
	Test_csr_dussm();
	Test_ssr_dusmv();

	Test_cg1();
	Test_bicg1();

	Test_mmio();

	Test_csr_check();
	Test_csr_check_sym();

	return 0;
}

}
