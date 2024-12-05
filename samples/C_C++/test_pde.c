/****************************************
 *                                      *
 *  C/C++ interface to XLPack           *
 *  Test program                        *
 *  Version 7.0 (January 31, 2023)      *
 *  (C) 2014-2023  K Technologies       *
 *                                      *
 ****************************************/

//#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>

#include "spcnumlib.h"
#include "pdecnumlib.h"

// Test fem2p

// Solve Laplace equation
//   d2u/dx2 + d2u/dy2 = 0
//   u(x, 0) = 0, u(0, y) = 0, u(x, 1) = x, u(1, y) = y (BC 1)

#define NX	8
#define NY	8
#define N	((NX + 1)*(NY + 1))
#define NE	(2*NX*NY)
#define NB	(2*(NX + NY))

void test_fem2p()
{
	char fname[] = "test.vtk";
	double x[N], y[N], z[N], p[N], q[N], f[N];
	int knc[4*NE], ks[3*NB], lb[NB];
	double scale = 1;
	int n, ne, nb, ldknc = 4, ldks = 3;
	// BC 1 (nb1 = nb)
	int nb1, ib1[NB], id[N]; double bv1[NB];
	// BC 2 (not used (nb2 = 0))
//	int nb2 = 0, ldks2 = 3, ks2[3*NB], ldalpha = 2, ldbeta = 2;
//	double alpha[2*NB], beta[2*NB];
	int nb2 = 0, ldks2 = 0, *ks2 = NULL, ldalpha = 0, ldbeta = 0;
	double *alpha = NULL, *beta = NULL;
	// FEM matrix
	double a[N*N], b[N], u[N];
	int ia[N + 1], ja[N*N], nnz = N*N, base = 0;
	//
	int iwork[2*N], info;

	printf("** fem2p\n");
	// Generate mesh data
	mesh23(NX, NY, &n, scale, scale, x, y, &ne, ldknc, knc, ldks, ks, lb, &nb, &info);
	printf("mesh23: info = %d\n", info);
	if (info != 0)
		return;
	if (n != N || ne != NE || nb != NB) {
		printf("invalid n, ne or nb: %d %d %d\n", n, ne, nb);
		return;
	}
	// Set BC 1
	for (int i = 0; i < n; i++)
		id[i] = 0;
	for (int i = 0; i < nb; i++) {
		for (int j = 0; j < 2; j++)
			id[ks[ldks*i + j + 1] - 1] = lb[i];
	}
	int k = 0;
	for (int i = 0; i < n; i++) {
		if (id[i] != 0) {
			ib1[k] = i + 1;
			if (id[i] == 1 || id[i] == 4)
				bv1[k] = 0;
			else if (id[i] == 2)
				bv1[k] = y[i];
			else if (id[i] == 3)
				bv1[k] = x[i];
			k++;
		}
	}
	nb1 = k;
	// Assemble FEM matrix
	for (int i = 0; i < n; i++) {
		p[i] = 1;
		q[i] = 0;
		f[i] = 0;
	}
	fem2p(n, ne, x, y, ldknc, knc, p, q, f, nb1, ib1, bv1, nb2, ldks2, ks2, ldalpha, alpha, ldbeta, beta, nnz, a, ia, ja, base, b, iwork, &info);
	printf("fem2p: info = %d\n", info);
	if (info != 0)
		return;
	// Solve equation by cg
	double tol = 1.0e-10, res, work[5*N];
	int maxiter = 100, iter, lwork = 5*N;
	for (int i = 0; i < n; i++)
		u[i] = 0;
	cg1('F', n, a, ia, ja, b, u, tol, maxiter, &iter, &res, lwork, work, &info);
	printf("cg1: info = %d\n", info);
	if (info != 0)
		return;
	printf("iter = %d, res = %g\n", iter, res);
	printf("Computed solution u(i):\n");
	for (int i = 0; i < n; i += NX + 1) {
		for (int j = 0; j < NX + 1; j++)
			printf("%10.5f", u[i + j]);
		printf("\n");
	}
	// Compute an error
	double averr = 0;
	for (int i = 0; i < n; i++)
		averr = averr + fabs(u[i] - x[i]*y[i]);
	averr = averr/n;
	printf("Av. abs. error = %g\n", averr);
	// Write to VTK file
	for (int i = 0; i < n; i++)
		z[i] = 0;
	writevtkug(fname, n, x, y, z, ne, ldknc, knc, u, &info);
	printf("writevtkug: info = %d\n", info);
	if (info != 0)
		return;
}

// Test gmsh22 I/O

void test_gmsh22()
{
	char fname[] = "test.mesh";
	double x[N], y[N], z[N];
	int knc[4*NE], ks[3*NB], lb[NB];
	double scale = 1;
	int n, ne, nb, ldknc = 4, ldks = 3;
	// gmsh data
	int ne2, ldkc2 = 4, kc2[4*(NE + NB)], ldlb2 = 2, lb2[2*(NE + NB)];
	int n3, ne3, ldkc3 = 4, kc3[4*(NE + NB)], ldlb3 = 2, lb3[2*(NE + NB)];
	double x3[N], y3[N], z3[N];
	// etc
	double err;
	int info;

	printf("** gmsh22\n");
	// Generate mesh data
	mesh23(NX, NY, &n, scale, scale, x, y, &ne, ldknc, knc, ldks, ks, lb, &nb, &info);
	printf("mesh23: info = %d\n", info);
	if (info != 0)
		return;
	printf("n = %d, ne = %d, nb = %d\n", N, NE, NB);
	if (n != N || ne != NE || nb != NB) {
		printf("invalid n, ne or nb: %d %d %d\n", n, ne, nb);
		return;
	}
	// Set gmsh data
	ne2 = ne + nb;
	for (int j = 0; j < ne; j++) {
		for (int i = 0; i < ldknc; i++)
			kc2[ldkc2*j + i] = knc[ldknc*j + i];
		lb2[ldlb2*j] = 1;
		lb2[ldlb2*j + 1] = 11;
	}
	for (int j = 0; j < nb; j++) {
		for (int i = 0; i < ldks; i++)
			kc2[ldkc2*(j + ne) + i] = ks[ldks*j + i];
		lb2[ldlb2*(j + ne)] = 1;
		lb2[ldlb2*(j + ne) + 1] = lb[j];
	}
	for (int i = 0; i < n; i++)
		z[i] = 0;
	// Write to file
	writegmsh22(fname, n, x, y, z, ne2, ldkc2, kc2, ldlb2, lb2, &info);
	printf("writegmsh22: info = %d\n", info);
	if (info != 0)
		return;
	// Read a file
	readgmsh22(fname, &n3, x3, y3, z3, &ne3, ldkc3, kc3, ldlb3, lb3, &info);
	printf("readgmsh22: info = %d\n", info);
	if (info != 0)
		return;
	// Compare data
	if (n3 != n || ne3 != ne2) {
		printf("unmatched data: n = %d, ne = %d\n", n3, ne3);
		return;
	}
	err = 0;
	for (int i = 0; i < n3; i++) {
		if (x3[i] != x[i])
			printf("unmatched x data: %d\n", i);
		if (y3[i] != y[i])
			printf("unmatched y data: %d\n", i);
		if (z3[i] != z[i])
			printf("unmatched z data: %d\n", i);
	}
	for (int j = 0; j < ne3; j++) {
		int ns = 3;
		if (kc3[ldkc3*j] == 2)
			ns = 4;
		for (int i = 0; i < ns; i++) {
			if (kc3[ldkc3*j + i] != kc2[ldkc2*j + i])
				printf("unmatched connection data: %d %d\n", i, j);
		}
		for (int i = 0; i < lb3[ldlb3*j]; i++) {
			if (lb3[ldlb3*j + i] != lb2[ldlb2*j + i])
				printf("unmatched label data: %d %d\n", i, j);
		}
	}
}

int main()
{
	test_fem2p();
	test_gmsh22();

	return 0;
}
