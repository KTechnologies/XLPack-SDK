/****************************************
 *                                      *
 *  C/C++ interface to XLPack           *
 *  Test program                        *
 *  Version 7.0 (January 31, 2023)      *
 *  (C) 2014-2023  K Technologies       *
 *                                      *
 ****************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "spcnumlib.h"

// Test SPUtils

void test_csr_dense()
{
	double a[] = {
		0.2, 0.11, 0, -0.11, -0.8,
		0.11, 0.93, 0.81, 0.93, -0.92,
		0, -0.81, 0.37, 0.81, 0,
		0.8, 0.92, 0, 0.9, 0.86 };
	double a_val[] = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int a_ptr[] = { 0, 3, 7, 9, 13, 16 };
	int a_ind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	double a2[5*4];
	int m = 5, n = 4, lda = 5, base = 0, info;

	printf("** test csr_dense\n");
	csr_dense(m, n, a_val, a_ptr, a_ind, base, lda, a2, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		for (int j = 0; j < n; j++) {
			for (int i = 0; i < m; i++) {
				if (a[j*lda + i] != a2[j*lda + i])
					printf("error: unmatched %d %d\n", i, j);
			}
		}
	}
}

void test_dense_csr()
{
	double a[] = {
		0.2, 0.11, 0, -0.11, -0.8,
		0.11, 0.93, 0.81, 0.93, -0.92,
		0, -0.81, 0.37, 0.81, 0,
		0.8, 0.92, 0, 0.9, 0.86 };
	double a_val[] = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int a_ptr[] = { 0, 3, 7, 9, 13, 16 };
	int a_ind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	double a2_val[16];
	int a2_ptr[6], a2_ind[16];
	int m = 5, n = 4, lda = 5, nnz = 16, base = 0, info;

	printf("** test dense_csr\n");
	dense_csr(m, n, lda, a, nnz, a2_val, a2_ptr, a2_ind, base, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		for (int i = 0; i < m + 1; i++) {
			if (a2_ptr[i] != a_ptr[i])
				printf("error: unmatched ptr %d\n", i);
		}
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != a_val[i])
				printf("error: unmatched val %d\n", i);
			if (a2_ind[i] != a_ind[i])
				printf("error: unmatched ind %d\n", i);
		}
	}
}

void test_csr_coo()
{
	double a_val[] = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int a_rowptr[] = { 0, 3, 7, 9, 13, 16 };
	int a_colind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	int a_rowind[] = { 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4 };
	double a2_val[16];
	int a2_rowind[16], a2_colind[16];
	int m = 5, n = 4, nnz = 16, base = 0, info;

	printf("** test csr_coo\n");
	csr_coo(m, n, a_val, a_rowptr, a_colind, base, a2_val, a2_rowind, a2_colind, base, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != a_val[i])
				printf("error: unmatched val %d\n", i);
			if (a2_rowind[i] != a_rowind[i])
				printf("error: unmatched rowind %d\n", i);
			if (a2_colind[i] != a_colind[i])
				printf("error: unmatched colind %d\n", i);
		}
	}
}

void test_coo_csr()
{
	double a_val[] = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int a_rowptr[] = { 0, 3, 7, 9, 13, 16 };
	int a_colind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	int a_rowind[] = { 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4 };
	double a2_val[16];
	int a2_rowptr[6], a2_colind[16];
	int m = 5, n = 4, nnz = 16, base = 0, info;

	printf("** test coo_csr\n");
	coo_csr(m, n, nnz, a_val, a_rowind, a_colind, base, a2_val, a2_rowptr, a2_colind, base, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		for (int i = 0; i < m + 1; i++) {
			if (a2_rowptr[i] != a_rowptr[i])
				printf("error: unmatched rowptr %d\n", i);
		}
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != a_val[i])
				printf("error: unmatched val %d\n", i);
			if (a2_colind[i] != a_colind[i])
				printf("error: unmatched colind %d\n", i);
		}
	}
}

void test_csr_ssr()
{
	double a_val[] = { 2.2, -0.11, -0.8, -0.11, 2.93, 0.81, -0.92, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int a_ptr[] = { 0, 3, 7, 9, 12 };
	int a_ind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 3 };
	double s_val[] = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int s_ptr[] = { 0, 1, 3, 5, 8 };
	int s_ind[] = { 0, 0, 1, 1, 2, 0, 1, 3 };
	double a2_val[8];
	int a2_ptr[5], a2_ind[8];
	int n = 4, nnz = 8, base = 0, info;

	printf("** test csr_ssr\n");
	csr_ssr('L', n, a_val, a_ptr, a_ind, base, nnz, a2_val, a2_ptr, a2_ind, base, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		for (int i = 0; i < n + 1; i++) {
			if (a2_ptr[i] != s_ptr[i])
				printf("error: unmatched ptr %d\n", i);
		}
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != s_val[i])
				printf("error: unmatched val %d\n", i);
			if (a2_ind[i] != s_ind[i])
				printf("error: unmatched ind %d\n", i);
		}
	}
}

void test_ssr_csr()
{
	double a_val[] = { 2.2, -0.11, -0.8, -0.11, 2.93, 0.81, -0.92, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int a_ptr[] = { 0, 3, 7, 9, 12 };
	int a_ind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 3 };
	double s_val[] = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int s_ptr[] = { 0, 1, 3, 5, 8 };
	int s_ind[] = { 0, 0, 1, 1, 2, 0, 1, 3 };
	double a2_val[12];
	int a2_ptr[5], a2_ind[12];
	int n = 4, nnz = 12, base = 0, info;

	printf("** test ssr_csr\n");
	ssr_csr('L', n, s_val, s_ptr, s_ind, base, nnz, a2_val, a2_ptr, a2_ind, base, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		for (int i = 0; i < n + 1; i++) {
			if (a2_ptr[i] != a_ptr[i])
				printf("error: unmatched ptr %d\n", i);
		}
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != a_val[i])
				printf("error: unmatched val %d\n", i);
			if (a2_ind[i] != a_ind[i])
				printf("error: unmatched ind %d\n", i);
		}
	}
}

void test_csc_csr()
{
	double a_val[] = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int a_ptr[] = { 0, 3, 7, 9, 13, 16 };
	int a_ind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	double at_val[] = { 0.2, 0.11, -0.11, -0.8, 0.11, 0.93, 0.81, 0.93, -0.92, -0.81, 0.37, 0.81, 0.8, 0.92, 0.9, 0.86 };
	int at_ptr[] = { 0, 4, 9, 12, 16 };
	int at_ind[] = { 0, 1, 3, 4, 0, 1, 2, 3, 4, 1, 2, 3, 0, 1, 3, 4 };
	double a2_val[16];
	int a2_ptr[5], a2_ind[16];
	int m = 5, n = 4, nnz = 16, base = 0, info;

	printf("** test csc_csr\n");
	csc_csr(n, m, a_val, a_ptr, a_ind, base, a2_val, a2_ptr, a2_ind, base, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		for (int i = 0; i < n + 1; i++) {
			if (a2_ptr[i] != at_ptr[i])
				printf("error: unmatched ptr %d\n", i);
		}
		for (int i = 0; i < nnz; i++) {
			if (a2_val[i] != at_val[i])
				printf("error: unmatched val %d\n", i);
			if (a2_ind[i] != at_ind[i])
				printf("error: unmatched ind %d\n", i);
		}
	}
}

// Test SPBLAS

void test_csr_dusmv()
{
	double a[] = {
		0.2, 0.11, 0, -0.11, -0.8,
		0.11, 0.93, 0.81, 0.93, -0.92,
		0, -0.81, 0.37, 0.81, 0,
		0.8, 0.92, 0, 0.9, 0.86 };
	double a_val[] = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int a_ptr[] = { 0, 3, 7, 9, 13, 16 };
	int a_ind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	double x[] = { 0.73, 0.86, 0.44, 0.87 };
	double y[] = { 0.28, 0.30, 0.69, 0.50, 0.15 };
	double y_bup[5];
	double alpha = 2.0, beta = 3.0;
	int m = 5, n = 4, lda = 5, base = 0, info;

	printf("** test csr_dusmv\n");
	for (int i = 0; i < m; i++)
		y_bup[i] = y[i];
	csr_dusmv('N', m, n, alpha, a_val, a_ptr, a_ind, base, x, 1, beta, y, 1, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		double err = 0;
		for (int i = 0; i < m; i++) {
			double s = 0;
			for (int j = 0; j < n; j++)
				s = s + a[lda*j + i]*x[j];
			err = max(err, fabs(y[i] - (alpha*s + beta*y_bup[i])));
		}
		printf("err = %g\n", err);
	}
}

void test_csr_dussv()
{
	double a[] = {
		2.2, -0.11, 0, -0.8,
		0, 2.93, 0.81, -0.92,
		0, 0, 2.37, 0,
		0, 0, 0, 2.86 };
	double s_val[] = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86 };
	int s_ptr[] = { 0, 1, 3, 5, 8 };
	int s_ind[] = { 0, 0, 1, 1, 2, 0, 1, 3 };
	double b[] = { 0.73, 0.86, 0.44, 0.87 };
	double x[4];
	int n = 4, lda = 4, base = 0, info;

	printf("** test csr_dussv\n");
	for (int i = 0; i < n; i++)
		x[i] = b[i];
	csr_dussv('L', 'N', 'N', n, s_val, s_ptr, s_ind, base, x, 1, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		double err = 0;
		for (int i = 0; i < n; i++) {
			double s = 0;
			for (int j = 0; j <= i; j++)
				s = s + a[lda*j + i]*x[j];
			err = max(err, fabs(b[i] - s));
		}
		printf("err = %g\n", err);
	}
}

void test_csr_dusmm()
{
	double a[] = {
		0.2, 0.11, 0, -0.11, -0.8,
		0.11, 0.93, 0.81, 0.93, -0.92,
		0, -0.81, 0.37, 0.81, 0,
		0.8, 0.92, 0, 0.9, 0.86 };
	double a_val[] = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int a_ptr[] = { 0, 3, 7, 9, 13, 16 };
	int a_ind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	double b[] = { 0.73, 0.86, 0.44, 0.87,
		0.0057, 0.89, 0.75, 0.18,
		0.71, 0.46, 0.69, 0.64 };
	double c[] = { 0.28, 0.30, 0.69, 0.50, 0.15,
		0.99, 0.24, 0.4, 0.11, 0.48,
		0.0029, 0.67, 0.38, 0.46, 0.014 };
	double c_bup[5*3];
	double alpha = 2, beta = 3;
	int m = 5, n = 4, l = 3, lda = 5, ldb = 4, ldc = 5, base = 0, info;

	printf("** test csr_dusmm\n");
	for (int i = 0; i < m*l; i++)
		c_bup[i] = c[i];
	csr_dusmm('N', 'C', m, n, l, alpha, a_val, a_ptr, a_ind, base, ldb, b, beta, ldc, c, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		double err = 0;
		for (int k = 0; k < l; k++) {
			for (int i = 0; i < m; i++) {
				double s = 0;
				for (int j = 0; j < n; j++)
					s = s + a[lda*j + i]*b[ldb*k + j];
				err = max(err, fabs(c[ldc*k + i] - (alpha*s + beta*c_bup[ldc*k + i])));
			}
		}
		printf("err = %g\n", err);
	}
}

void test_csr_dussm()
{
	double a[] = {
		2.2, -0.11, 0, -0.8,
		0, 2.93, 0.81, -0.92,
		0, 0, 2.37, 0,
		0, 0, 0, 2.86 };
	double s_val[] = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86 };
	int s_ptr[] = { 0, 1, 3, 5, 8 };
	int s_ind[] = { 0, 0, 1, 1, 2, 0, 1, 3 };
	double b[] = { 0.73, 0.86, 0.44, 0.87,
		0.0057, 0.89, 0.75, 0.18,
		0.71, 0.46, 0.69, 0.64 };
	double x[4*3];
	int n = 4, nrhs = 3, lda = 4, ldb = 4, ldx = 4, base = 0, info;

	printf("** test csr_dussm\n");
	for (int i = 0; i < n*nrhs; i++)
		x[i] = b[i];
	csr_dussm('L', 'N', 'N', 'C', n, nrhs, s_val, s_ptr, s_ind, base, ldx, x, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		double err = 0;
		for (int k = 0; k < nrhs; k++) {
			for (int i = 0; i < n; i++) {
				double s = 0;
				for (int j = 0; j <= i; j++)
					s = s + a[lda*j + i]*x[ldx*k + j];
				err = max(err, fabs(b[ldb*k + i] - s));
			}
		}
		printf("err = %g\n", err);
	}
}

void test_ssr_dusmv()
{
	double a[] = {
		2.2, -0.11, 0, -0.8,
		-0.11, 2.93, 0.81, -0.92,
		0, 0.81, 2.37, 0,
		-0.8, -0.92, 0, 2.86 };
	double s_val[] = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86 };
	int s_ptr[] = { 0, 1, 3, 5, 8 };
	int s_ind[] = { 0, 0, 1, 1, 2, 0, 1, 3 };
	double x[] = { 0.73, 0.86, 0.44, 0.87 };
	double y[] = { 0.28, 0.30, 0.69, 0.50 };
	double y_bup[4];
	double alpha = 2.0, beta = 3.0;
	int n = 4, lda = 4, base = 0, info;

	printf("** test ssr_dusmv\n");
	for (int i = 0; i < n; i++)
		y_bup[i] = y[i];
	ssr_dusmv('L', n, alpha, s_val, s_ptr, s_ind, base, x, 1, beta, y, 1, &info);
	printf("info = %d\n", info);
	if (info == 0) {
		double err = 0;
		for (int i = 0; i < n; i++) {
			double s = 0;
			for (int j = 0; j < n; j++)
				s = s + a[lda*j + i]*x[j];
			err = max(err, fabs(y[i] - (alpha*s + beta*y_bup[i])));
		}
		printf("err = %g\n", err);
	}
}

// Test iterative solvers

void test_cg1()
{
	double s_val[] = { 2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int s_ptr[] = { 0, 1, 3, 5, 8 };
	int s_ind[] = { 0, 0, 1, 1, 2, 0, 1, 3 };
	double b[] = { 1.433, 1.3137, 2.3799, -0.7992 };
	double x[] = { 0, 0, 0, 0 };
	double tol = 1.0e-10, res;
	int maxiter = 100, iter;
	double work[5*4]; int lwork = 5*4;
	int n = 4, info;

	printf("** test cg1\n");
	cg1('L', n, s_val, s_ptr, s_ind, b, x, tol, maxiter, &iter, &res, lwork, work, &info);
	printf("x =\n");
	printf("%g  %g  %g  %g\n", x[0], x[1], x[2], x[3]);
	printf("iter = %d, res = %g, info = %d\n", iter, res, info);
}

void test_bicg1()
{
	double a_val[] = { 2.2, -0.11, -0.8, 0.11, 2.93, 0.81, -0.92, 0.81, 2.37, 0.8, 0.92, 0.86 };
	int a_ptr[] = { 0, 3, 7, 9, 12 };
	int a_ind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 3 };
	double b[] = { 1.433, 1.4325, 2.3799, 0.2488 };
	double x[] = { 0, 0, 0, 0 };
	double tol = 1.0e-10, res;
	int maxiter = 100, iter;
	double work[8*4]; int lwork = 8*4;
	int n = 4, info;

	printf("** test bicg1\n");
	bicg1(n, a_val, a_ptr, a_ind, b, x, tol, maxiter, &iter, &res, lwork, work, &info);
	printf("x =\n");
	printf("%g  %g  %g  %g\n", x[0], x[1], x[2], x[3]);
	printf("iter = %d, res = %g, info = %d\n", iter, res, info);
}

// Test I/O routines

void test_mmio()
{
	char fname[] = "test.mtx", matcode[5];
	double a_val[] = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int a_ptr[] = { 0, 3, 7, 9, 13, 16 };
	int a_ind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	double b_val[16];
	int b_ptr[6], b_ind[16];
	int m = 5, n = 4, nnz = 16, info;
	int skip = 0, base = 0, fchk = 0, format = 0, sort = 0;

	printf("** test mmio\n");
	// Write a (CSR) to MM file
	mm_write(fname, "MCRS", m, n, nnz, a_val, nnz, a_ptr, m + 1, a_ind, nnz, base, fchk, format, &info);
	printf("mm_write: info = %d\n", info);
	if (info != 0)
		return;
	// Read MM file info
	mm_read(fname, matcode, &m, &n, &nnz, b_val, nnz, b_ptr, m + 1, b_ind, nnz, skip, base, format, sort, &info);
	printf("mm_read: info = %d\n", info);
	if (info != 0)
		return;
	matcode[4] = '\0';
	printf("m = %d, n = %d, nnz = %d, matcode = %s\n", m, n, nnz, matcode);
	// Compare data
	for (int i = 0; i < m + 1; i++) {
		if (b_ptr[i] != a_ptr[i])
			printf("error: unmatched ptr %d\n", i);
	}
	for (int i = 0; i < nnz; i++) {
		if (b_val[i] != a_val[i])
			printf("error: unmatched val %d\n", i);
		if (b_ind[i] != a_ind[i])
			printf("error: unmatched ind %d\n", i);
	}
}

// Test chexk routines

void test_csr_check()
{
	double a_val[] = { 0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92, 0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86 };
	int a_ptr[] = { 0, 3, 7, 9, 13, 16 };
	int a_ind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3 };
	int m = 5, n = 4, result[10], info;

	printf("** test csr_check\n");
	csr_check(m, n, a_val, a_ptr, a_ind, result, &info);
	printf("csr_check: info = %d\n", info);
	printf("  base = %d\n", result[0]);
	printf("  nnz = %d\n", result[1]);
	printf("  nnz_lower = %d\n", result[2]);
	printf("  nnz_upper = %d\n", result[3]);
	printf("  nnz_diag = %d\n", result[4]);
	printf("  n_zero = %d\n", result[5]);
	printf("  n_zero_diag = %d\n", result[6]);
	printf("  empty_rows = %d\n", result[7]);
	printf("  unsorted_rows = %d\n", result[8]);
	printf("  invalid_inds = %d\n", result[9]);
}

void test_csr_check_sym()
{
	double a_val[] = { 2.2, -0.11, -0.8, -0.11, 2.93, 0.81, -0.92, 0.81, 2.37, -0.8, -0.92, 0.86 };
	double a2_val[] = { 2.2, 0.11, -0.8, -0.11, 2.93, 0.81, -0.92, 0.81, 2.37, -0.8, -0.92, 0.86 };
	int a_ptr[] = { 0, 3, 7, 9, 12 };
	int a_ind[] = { 0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 3 };
	int n = 4, info;

	printf("** test csr_check_sym\n");
	csx_check_sym(n, a_val, a_ptr, a_ind, &info);
	printf("csx_check_sym: info = %d\n", info);
	csx_check_sym(n, a2_val, a_ptr, a_ind, &info);
	printf("csx_check_sym: info = %d\n", info);
}

// Main program

int main()
{
	test_csr_dense();
	test_dense_csr();
	test_csr_coo();
	test_coo_csr();
	test_csr_ssr();
	test_ssr_csr();
	test_csc_csr();

	test_csr_dusmv();
	test_csr_dussv();
	test_csr_dusmm();
	test_csr_dussm();
	test_ssr_dusmv();

	test_cg1();
	test_bicg1();

	test_mmio();

	test_csr_check();
	test_csr_check_sym();

	return 0;
}
