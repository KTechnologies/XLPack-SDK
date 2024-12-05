/****************************************
 *                                      *
 *  C/C++ interface to XLPack           *
 *  Test program                        *
 *  Version 7.0 (January 31, 2023)      *
 *  (C) 2014-2023  K Technologies       *
 *                                      *
 ****************************************/

using System;
using System.IO;
using static System.Math;
using static XLPack;

public class Test
{

// Test fem2p

// Solve Laplace equation
//   d2u/dx2 + d2u/dy2 = 0
//   u(x, 0) = 0, u(0, y) = 0, u(x, 1) = x, u(1, y) = y (BC 1)

static void Test_fem2p()
{
	int nx = 8, ny = 8;
	int n = (nx + 1)*(ny + 1), ne = 2*nx*ny, nb = 2*(nx + ny);
	string fname = "test.vtk";
	double[] x = new double[n];
	double[] y = new double[n];
	double[] z = new double[n];
	double[] p = new double[n];
	double[] q = new double[n];
	double[] f = new double[n];
	int[] knc = new int[4*ne];
	int[] ks = new int[3*nb];
	int[] lb = new int[nb];
	double scale = 1;
	int ldknc = 4, ldks = 3;
	// BC 1 (nb1 = nb)
	int nb1;
	int[] ib1 = new int[nb];
	int[] id = new int[n];
	double[] bv1 = new double[nb];
	// BC 2 (not used (nb2 = 0))
	int nb2 = 0, ldks2 = 3, ldalpha = 2, ldbeta = 2;
	int[] ks2 = new int[3*nb];
	double[] alpha = new double[2*nb];
	double[] beta = new double[2*nb];
	// FEM matrix
	double[] a = new double[n*n];
	double[] b = new double[n];
	double[] u = new double[n];
	int[] ia = new int[n + 1];
	int[] ja = new int[n*n];
	int nnz = n*n, index = 0;
	int info;

	Console.WriteLine("** fem2p");
	// Generate mesh data
	int n_o, ne_o, nb_o;
	mesh23(nx, ny, out n_o, scale, scale, x, y, out ne_o, ldknc, knc, ldks, ks, lb, out nb_o, out info);
	Console.WriteLine("mesh23: info = {0}", info);
	if (info != 0)
		return;
	if (n_o != n || ne_o != ne || nb_o != nb) {
		Console.WriteLine("invalid n, ne or nb: {0} {1} {2}", n_o, ne_o, nb_o);
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
	fem2p(n, ne, x, y, ldknc, knc, p, q, f, nb1, ib1, bv1, nb2, ldks2, ks2, ldalpha, alpha, ldbeta, beta, nnz, a, ia, ja, index, b, out info);
	Console.WriteLine("fem2p: info = {0}", info);
	if (info != 0)
		return;
	// Solve equation by cg
	double tol = 1.0e-10, res;
	int maxiter = 100, iter;
	for (int i = 0; i < n; i++)
		u[i] = 0;
	cg1('F', n, a, ia, ja, b, u, tol, maxiter, out iter, out res, out info);
	Console.WriteLine("cg1: info = {0}", info);
	if (info != 0)
		return;
	Console.WriteLine("iter = {0}, res = {1}", iter, res);
	Console.WriteLine("Computed solution u(i):");
	for (int i = 0; i < n; i += nx + 1) {
		for (int j = 0; j < nx + 1; j++)
			Console.Write("{0,10:F5}", u[i + j]);
		Console.WriteLine();
	}
	// Compute an error
	double averr = 0;
	for (int i = 0; i < n; i++)
		averr = averr + Abs(u[i] - x[i]*y[i]);
	averr = averr/n;
	Console.WriteLine("Av. abs. error = {0}", averr);
	// Write to VTK file
	for (int i = 0; i < n; i++)
		z[i] = 0;
	writevtkug(fname, n, x, y, z, ne, ldknc, knc, u, out info);
	Console.WriteLine("writevtkug: info = {0}", info);
	if (info != 0)
		return;
}

	// Test gmsh22 I/O

static void Test_gmsh22()
{
	int nx = 8, ny = 8;
	int n = (nx + 1)*(ny + 1), ne = 2*nx*ny, nb = 2*(nx + ny);
	string fname = "test.mesh";
	double[] x = new double[n];
	double[] y = new double[n];
	double[] z = new double[n];
	int[] knc = new int[4*ne];
	int[] ks = new int[3*nb];
	int[] lb = new int[nb];
	double scale = 1;
	int ldknc = 4, ldks = 3;
	// gmsh data
	int ne2, ldkc2 = 4, ldlb2 = 2;
	int[] kc2 = new int[4*(ne + nb)];
	int[] lb2 = new int[2*(ne + nb)];
	int n3, ne3, ldkc3 = 4, ldlb3 = 2;
	int[] kc3 = new int[4*(ne + nb)];
	int[] lb3 = new int[2*(ne + nb)];
	double[] x3 = new double[n];
	double[] y3 = new double[n];
	double[] z3 = new double[n];
	int info;

	Console.WriteLine("** gmsh22");
	// Generate mesh data
	int n_o, ne_o, nb_o;
	mesh23(nx, ny, out n_o, scale, scale, x, y, out ne_o, ldknc, knc, ldks, ks, lb, out nb_o, out info);
	Console.WriteLine("mesh23: info = {0}", info);
	if (info != 0)
		return;
	Console.WriteLine("n = {0}, ne = {1}, nb = {2}", n, ne, nb);
	if (n_o != n || ne_o != ne || nb_o != nb) {
		Console.WriteLine("invalid n, ne or nb: {0} {1} {2}", n_o, ne_o, nb_o);
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
	writegmsh22(fname, n, x, y, z, ne2, ldkc2, kc2, ldlb2, lb2, out info);
	Console.WriteLine("writegmsh22: info = {0}", info);
	if (info != 0)
		return;
	// Read a file
	readgmsh22(fname, out n3, x3, y3, z3, out ne3, ldkc3, kc3, ldlb3, lb3, out info);
	Console.WriteLine("readgmsh22: info = {0}", info);
	if (info != 0)
		return;
	// Compare data
	if (n3 != n || ne3 != ne2) {
		Console.WriteLine("unmatched data: n = {0}, ne = {1}", n3, ne3);
		return;
	}
	for (int i = 0; i < n3; i++) {
		if (x3[i] != x[i])
			Console.WriteLine("unmatched x data: {0}", i);
		if (y3[i] != y[i])
			Console.WriteLine("unmatched y data: {0}", i);
		if (z3[i] != z[i])
			Console.WriteLine("unmatched z data: {0}", i);
	}
	for (int j = 0; j < ne3; j++) {
		int ns = 3;
		if (kc3[ldkc3*j] == 2)
			ns = 4;
		for (int i = 0; i < ns; i++) {
			if (kc3[ldkc3*j + i] != kc2[ldkc2*j + i])
				Console.WriteLine("unmatched connection data: {0} {1}", i, j);
		}
		for (int i = 0; i < lb3[ldlb3*j]; i++) {
			if (lb3[ldlb3*j + i] != lb2[ldlb2*j + i])
				Console.WriteLine("unmatched label data: {0} {1}", i, j);
		}
	}
}

public static int Main(string[] args)
{
	Test_fem2p();
	Test_gmsh22();

	return 0;
}

}
