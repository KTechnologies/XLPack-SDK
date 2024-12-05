/****************************************
 *                                      *
 *  C/C++ interface to XLPack           *
 *  Test program                        *
 *  Version 7.0 (January 31, 2023)      *
 *  (C) 2014-2023  K Technologies       *
 *                                      *
 ****************************************/

#include <cmath>
#include <iostream>

#include "cnumlib.h"
#include "Matrix.h"
#include "MatrixIO.h"

// Refer to the following book for Matrix.h and MatrixIO.h:
//   Bjarne Stroustrup "Programming -- Principles and Practice Using C++, 2nd Ed." Addison-Wesley, 2014.
// Matrix.h and MatrixIO.h can be downloaded from the author's homepage https://www.stroustrup.com

using namespace std;
using namespace Numeric_lib;

// C++ library using Matrix

void Dgesv(int n, Matrix<double,2>& a, Matrix<double>& b, int& info)
{
	int *ipiv = new int[n];
	int lda = a.dim2(), nrhs = 1, ldb = n;
	dgesv(n, nrhs, lda, a.data(), ipiv, ldb, b.data(), &info);
	delete[] ipiv;
}

double Dlange(char norm, int m, int n, Matrix<double,2>& a)
{
	double *work = new double[m];
	int lda = a.dim2();
	double anorm = dlange(norm, m, n, lda, a.data(), work);
	delete[] work;
	return anorm;
}

void Dgecon(char norm, int n, Matrix<double,2>& a, double anorm, double& rcond, int& info)
{
	double *work = new double[4*n];
	int *iwork = new int[n];
	int lda = a.dim2();
	dgecon(norm, n, lda, a.data(), anorm, &rcond, work, iwork, &info);
	delete[] work;
	delete[] iwork;
}

void Dposv(char uplo, int n, Matrix<double,2>& a, Matrix<double>& b, int& info)
{
	int lda = a.dim2(), nrhs = 1, ldb = n;
	dposv(uplo, n, nrhs, lda, a.data(), ldb, b.data(), &info);
}

double Dlansy(char norm, char uplo, int n, Matrix<double,2>& a)
{
	double *work = new double[n];
	int lda = a.dim2();
	double anorm = dlansy(norm, uplo, n, lda, a.data(), work);
	delete[] work;
	return anorm;
}

void Dpocon(char uplo, int n, Matrix<double,2>& a, double anorm, double& rcond, int& info)
{
	double *work = new double[3*n];
	int *iwork = new int[n];
	int lda = a.dim2();
	dpocon(uplo, n, lda, a.data(), anorm, &rcond, work, iwork, &info);
	delete[] work;
	delete[] iwork;
}

void Dsyev(char jobz, char uplo, int n, Matrix<double,2>& a, Matrix<double>& w, int& info)
{
	int lwork = -1;
	double *work = new double[1];
	int lda = a.dim2();
	dsyev(jobz, uplo, n, lda, a.data(), w.data(), work, lwork, &info);
	if (info == 0) {
		lwork = (int)work[0];
		delete[] work;
		work = new double[lwork];
		dsyev(jobz, uplo, n, lda, a.data(), w.data(), work, lwork, &info);
	}
	delete[] work;
}

void Dgels(char trans, int m, int n, Matrix<double,2>& a, Matrix<double>& b, int& info)
{
	int lwork = -1;
	double *work = new double[1];
	int lda = a.dim2(), nrhs = 1, ldb = m;
	dgels(trans, m, n, nrhs, lda, a.data(), ldb, b.data(), work, lwork, &info);
	if (info == 0) {
		lwork = (int)work[0];
		delete[] work;
		work = new double[lwork];
		dgels(trans, m, n, nrhs, lda, a.data(), ldb, b.data(), work, lwork, &info);
	}
	delete[] work;
}

void Dgecov(int job, int n, Matrix<double,2>& a, Matrix<double>& ci, int& info)
{
	int lda = a.dim2();
	dgecov(job, n, lda, a.data(), ci.data(), &info);
}

// Test programs

void test_Dgesv()
{
	const int n = 3;
	double init_a[n][n] = {
		{ 0.2, -0.32, -0.8 },
		{ -0.11, 0.81, -0.92 },
		{ -0.93, 0.37, -0.29 }
	};
	Matrix<double, 2> a(init_a);
	double init_b[n] = { -0.3727, 0.4319, -1.4247 };
	Matrix<double> b(init_b);
	double anorm, rcond = 0;
	int info;

	anorm = Dlange('1', n, n, a);
	Dgesv(n, a, b, info);
	if (info == 0)
		Dgecon('1', n, a, anorm, rcond, info);
	cout << "** Dgesv\n";
	cout << "x = " << b << endl;
	cout << "rcond = " << rcond << ", info = " << info << endl;
}

void test_Dposv()
{
	const int n = 3;
	double init_a[n][n] = {
		{ 2.2, 0.0, 0.0 },
		{ -0.11, 2.93, 0.0 },
		{ -0.32, 0.81, 2.37 }
	};
	Matrix<double, 2> a(init_a);
	double init_b[n] = { -1.566, -2.8425, -1.1765 };
	Matrix<double> b(init_b);
	double anorm, rcond = 0;
	int info;

	anorm = Dlansy('1', 'U', n, a);
	Dposv('U', n, a, b, info);
	if (info == 0)
		Dpocon('U', n, a, anorm, rcond, info);
	cout << "** Dposv\n";
	cout << "x = " << b << endl;
	cout << "rcond = " << rcond << ", info = " << info << endl;
}

void test_Dsyev()
{
	const int n = 3;
	double init_a[n][n] = {
		{ 2.2, 0.0, 0.0 },
		{ -0.11, 2.93, 0.0 },
		{ -0.32, 0.81, 2.37 }
	};
	Matrix<double, 2> a(init_a);
	Matrix<double> w(n);
	int info;

	Dsyev('V', 'U', n, a, w, info);
	cout << "** Dsyev\n";
	cout << "Eigenvalues =\n";
	cout << w << endl;
	cout << "Eigenvectors =\n";
	cout << a << endl;
	cout << "info = " << info << endl;
}

void test_Dgels()
{
	const int m = 4, n = 2;
	Matrix<double,2> a(n, m);
	double init_x[m] = { 0.2, 118.2, 337.4, 884.6 };
	Matrix<double> x(init_x);
	double init_y[m] = { 0.1, 118.1, 338.8, 888.0 };
	Matrix<double> y(init_y);
	Matrix<double> ci(n);
	double s;
	int info;

	for (int i = 0; i < m; i++) {
		a(0, i) = 1;
		a(1, i) = x(i);
	}
	Dgels('N', m, n, a, y, info);
	cout << "** Dgels\n";
	cout << "a0 = " << y(0) << ", a1 = " << y(1) << endl;
	cout << "info = " << info << endl;
	if (info == 0) {
		Dgecov(0, n, a, ci, info);
		s = 0.0;
		for (int i = n; i < m; i++)
			s = s + y(i)*y(i);
		s = s / (m - n);
		cout << "Std. dev. = " << sqrt(s*ci(0)) << ", " << sqrt(s*ci(1)) << endl;
		cout << "info = " << info << endl;
	}
}

int main()
{
	test_Dgesv();
	test_Dposv();
	test_Dsyev();
	test_Dgels();
	return 0;
}
