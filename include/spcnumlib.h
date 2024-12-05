/****************************************
 *                                      *
 *  XLPack Numerical Library            *
 *  Version 7.0 (January 31, 2023)      *
 *  (C) 2014-2023  K Technologies       *
 *                                      *
 ****************************************/
#pragma once

#include "cnumlib_complex.h"

#if !defined(_SPCNUMLIB_NO_MANGLING)
#include "spcnumlib_mangling.h"
#endif

#if defined(__cplusplus)
extern "C" {
#endif

/*
 * Solver
 */
extern void _cg(int n, void (*matvec)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _cg1(char uplo, int n, const double val[], const int rowptr[], const int colind[], const double b[], double x[], double tol, int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _cr(int n, void (*matvec)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int mode, int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _cg_r(int n, const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);
extern void _cr_r(int n, const double b[], double x[], int mode, int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);

extern void _sor(int n, void (*matvec)(int, const double [], double []), void (*matsol)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _sor_r(int n, const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);

extern void _fom(int n, void (*matvec)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _diom(int n, void (*matvec)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _fgmres(int n, void (*matvec)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _dqgmres(int n, void (*matvec)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _gcr(int n, void (*matvec)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _orthomin(int n, void (*matvec)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _bicg(int n, void (*matvec)(int, const double [], double []), void (*matvectrans)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*psolvetrans)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _bicg1(int n, const double val[], const int rowptr[], const int colind[], const double b[], double x[], double tol, int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _cgs(int n, void (*matvec)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _gpbicg(int n, void (*matvec)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int mode, int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _qmr(int n, void (*matvec)(int, const double [], double []), void (*matvectrans)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*psolvetrans)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info);
extern void _tfqmr(int n, void (*matvec)(int, const double [], double []), void (*psolve)(int, const double [], double []), void (*chkconv)(int, const double [], double, int, int *), const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info);

extern void _fom_r(int n, const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);
extern void _diom_r(int n, const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);
extern void _fgmres_r(int n, const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);
extern void _dqgmres_r(int n, const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);
extern void _gcr_r(int n, const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);
extern void _orthomin_r(int n, const double b[], double x[], int m, int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);
extern void _bicg_r(int n, const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);
extern void _cgs_r(int n, const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);
extern void _gpbicg_r(int n, const double b[], double x[], int mode, int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);
extern void _qmr_r(int n, const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);
extern void _tfqmr_r(int n, const double b[], double x[], int maxiter, int *iter, double *res, int lwork, double work[], int *info, double xx[], double yy[], int *irev);

extern void _z_cg(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_cr(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int mode, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_cocg(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_cocr(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int mode, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);

extern void _z_cg_r(int n, const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_cr_r(int n, const doublecomplex b[], doublecomplex x[], int mode, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_cocg_r(int n, const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_cocr_r(int n, const doublecomplex b[], doublecomplex x[], int mode, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);

extern void _z_sor(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*matsol)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_sor_r(int n, const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);

extern void _z_fom(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_diom(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_fgmres(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_dqgmres(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_gcr(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_orthomin(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_bicg(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*matvectrans)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*psolvetrans)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_cgs(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_gpbicg(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int mode, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_qmr(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*matvectrans)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*psolvetrans)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);
extern void _z_tfqmr(int n, void (*matvec)(int, const doublecomplex [], doublecomplex []), void (*psolve)(int, const doublecomplex [], doublecomplex []), void (*chkconv)(int, const doublecomplex [], double, int, int *), const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info);

extern void _z_fom_r(int n, const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_diom_r(int n, const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_fgmres_r(int n, const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_dqgmres_r(int n, const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_gcr_r(int n, const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_orthomin_r(int n, const doublecomplex b[], doublecomplex x[], int m, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_bicg_r(int n, const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_cgs_r(int n, const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_gpbicg_r(int n, const doublecomplex b[], doublecomplex x[], int mode, int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_qmr_r(int n, const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);
extern void _z_tfqmr_r(int n, const doublecomplex b[], doublecomplex x[], int maxiter, int *iter, double *res, int lwork, doublecomplex work[], int *info, doublecomplex xx[], doublecomplex yy[], int *irev);

/*
 * Precon
 */
extern void _csx_ds(int n, const double val[], const int ptr[], const int ind[], int base, double d[], int *info);
extern void _csx_ds_solve(int n, const double d[], const double b[], double x[], int *info);

extern void _csx_ssor(int n, double omega, const double val[], const int ptr[], const int ind[], int base, double d[], int iwork[], int *info);
extern void _csc_ssor_solve(char trans, int n, double omega, const double val[], const int colptr[], const int rowind[], int base, const double d[], const double b[], double x[], int *info);
extern void _csr_ssor_solve(char trans, int n, double omega, const double val[], const int rowptr[], const int rowind[], int base, const double d[], const double b[], double x[], int *info);

extern void _ssc_ssor_solve(char uplo, int n, double omega, const double val[], const int colptr[], const int rowind[], int base, const double d[], const double b[], double x[], int *info);
extern void _ssr_ssor_solve(char uplo, int n, double omega, const double val[], const int rowptr[], const int colind[], int base, const double d[], const double b[], double x[], int *info);

void _ssc_ic0(char uplo, int n, const double val[], const int colptr[], const int rowind[], int base, double val2[], int idiag[], double work[], int *info);
void _ssr_ic0(char uplo, int n, const double val[], const int rowptr[], const int colind[], int base, double val2[], int idiag[], double work[], int *info);
void _ssc_ic_solve(char uplo, int n, const double val[], const int colptr[], const int rowind[], int base, const int idiag[], const double b[], double x[], int *info);
void _ssr_ic_solve(char uplo, int n, const double val[], const int rowptr[], const int colind[], int base, const int idiag[], const double b[], double x[], int *info);

extern void _csc_ilu0(int n, const double val[], const int colptr[], const int rowind[], int base, double val2[], double d[], int *info);
extern void _csr_ilu0(int n, const double val[], const int rowptr[], const int colind[], int base, double val2[], double d[], int *info);
extern void _csc_ilu(int n, const double val[], const int colptr[], const int rowind[], int base, int p, int nnz2, double val2[], int colptr2[], int rowind2[], int base2, double d[], double work[], int iwork[], int *info);
extern void _csr_ilu(int n, const double val[], const int rowptr[], const int colind[], int base, int p, int md, int nnz2, double val2[], int rowptr2[], int colind2[], int base2, double d[], double work[], int iwork[], int *info);
extern void _csc_ilu_solve(char trans, int n, const double val[], const int rowind[], const int colptr[], int base, const double d[], const double b[], double x[], int *info);
extern void _csr_ilu_solve(char trans, int n, const double val[], const int colind[], const int rowptr[], int base, const double d[], const double b[], double x[], int *info);

extern void _z_csx_ds(int n, const doublecomplex val[], const int ptr[], const int ind[], int base, doublecomplex d[], int *info);
extern void _z_csx_ds_solve(int n, const doublecomplex d[], const doublecomplex b[], doublecomplex x[], int *info);

extern void _z_csx_ssor(int n, double omega, const doublecomplex val[], const int ptr[], const int ind[], int base, doublecomplex d[], int iwork[], int *info);
extern void _z_csc_ssor_solve(char trans, int n, double omega, const doublecomplex val[], const int colptr[], const int rowind[], int base, const doublecomplex d[], const doublecomplex b[], doublecomplex x[], int *info);
extern void _z_csr_ssor_solve(char trans, int n, double omega, const doublecomplex val[], const int rowptr[], const int colind[], int base, const doublecomplex d[], const doublecomplex b[], doublecomplex x[], int *info);

extern void _z_ssc_ssor_solve(char uplo, int n, double omega, const doublecomplex val[], const int colptr[], const int rowind[], int base, const doublecomplex d[], const doublecomplex b[], doublecomplex x[], int *info);
extern void _z_ssr_ssor_solve(char uplo, int n, double omega, const doublecomplex val[], const int rowptr[], const int colind[], int base, const doublecomplex d[], const doublecomplex b[], doublecomplex x[], int *info);
extern void _z_hsc_ssor_solve(char uplo, int n, double omega, const doublecomplex val[], const int colptr[], const int rowind[], int base, const doublecomplex d[], const doublecomplex b[], doublecomplex x[], int *info);
extern void _z_hsr_ssor_solve(char uplo, int n, double omega, const doublecomplex val[], const int rowptr[], const int colind[], int base, const doublecomplex d[], const doublecomplex b[], doublecomplex x[], int *info);

void _z_hsc_ic0(char uplo, int n, const doublecomplex val[], const int colptr[], const int rowind[], int base, doublecomplex val2[], int idiag[], doublecomplex work[], int *info);
void _z_hsr_ic0(char uplo, int n, const doublecomplex val[], const int rowptr[], const int colind[], int base, doublecomplex val2[], int idiag[], doublecomplex work[], int *info);
void _z_hsc_ic_solve(char uplo, int n, const doublecomplex val[], const int colptr[], const int rowind[], int base, const int idiag[], const doublecomplex b[], doublecomplex x[], int *info);
void _z_hsr_ic_solve(char uplo, int n, const doublecomplex val[], const int rowptr[], const int colind[], int base, const int idiag[], const doublecomplex b[], doublecomplex x[], int *info);

extern void _z_csc_ilu0(int n, const doublecomplex val[], const int colptr[], const int rowind[], int base, doublecomplex val2[], doublecomplex d[], int *info);
extern void _z_csr_ilu0(int n, const doublecomplex val[], const int rowptr[], const int colind[], int base, doublecomplex val2[], doublecomplex d[], int *info);
extern void _z_csc_ilu(int n, const doublecomplex val[], const int colptr[], const int rowind[], int base, int p, int nnz2, doublecomplex val2[], int colptr2[], int rowind2[], int base2, doublecomplex d[], doublecomplex work[], int iwork[], int *info);
extern void _z_csr_ilu(int n, const doublecomplex val[], const int rowptr[], const int colind[], int base, int p, int md, int nnz2, doublecomplex val2[], int rowptr2[], int colind2[], int base2, doublecomplex d[], doublecomplex work[], int iwork[], int *info);
extern void _z_csc_ilu_solve(char trans, int n, const doublecomplex val[], const int rowind[], const int colptr[], int base, const doublecomplex d[], const doublecomplex b[], doublecomplex x[], int *info);
extern void _z_csr_ilu_solve(char trans, int n, const doublecomplex val[], const int colind[], const int rowptr[], int base, const doublecomplex d[], const doublecomplex b[], doublecomplex x[], int *info);

/*
 * SpBlas
 */
void _csc_dusmv(char trans, int m, int n, double alpha, const double val[], const int colptr[], const int rowind[], int base, const double x[], int incx, double beta, double y[], int incy, int *info);
void _csr_dusmv(char trans, int m, int n, double alpha, const double val[], const int rowptr[], const int colind[], int base, const double x[], int incx, double beta, double y[], int incy, int *info);
void _ssc_dusmv(char uplo, int n, double alpha, const double val[], const int colptr[], const int rowind[], int base, const double x[], int incx, double beta, double y[], int incy, int *info);
void _ssr_dusmv(char uplo, int n, double alpha, const double val[], const int rowptr[], const int colind[], int base, const double x[], int incx, double beta, double y[], int incy, int *info);

void _csc_dussv(char uplo, char trans, char diag, int n, const double val[], const int colptr[], const int rowind[], int base, double x[], int incx, int *info);
void _csr_dussv(char uplo, char trans, char diag, int n, const double val[], const int rowptr[], const int colind[], int base, double x[], int incx, int *info);

void _csc_dusmm(char trans, char order, int m, int n, int nrhs, double alpha, const double val[], const int colptr[], const int rowind[], int base, int ldb, const double b[], double beta, int ldc, double c[], int *info);
void _csr_dusmm(char trans, char order, int m, int n, int nrhs, double alpha, const double val[], const int rowptr[], const int colind[], int base, int ldb, const double b[], double beta, int ldc, double c[], int *info);
void _csc_dussm(char uplo, char trans, char diag, char order, int n, int nrhs, const double val[], const int colptr[], const int rowind[], int base, int ldx, double x[], int *info);
void _csr_dussm(char uplo, char trans, char diag, char order, int n, int nrhs, const double val[], const int rowptr[], const int colind[], int base, int ldx, double x[], int *info);

void _csc_zusmv(char trans, int m, int n, doublecomplex alpha, const doublecomplex val[], const int colptr[], const int rowind[], int base, const doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy, int *info);
void _csr_zusmv(char trans, int m, int n, doublecomplex alpha, const doublecomplex val[], const int rowptr[], const int colind[], int base, const doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy, int *info);
void _hsc_zusmv(char uplo, int n, doublecomplex alpha, const doublecomplex val[], const int colptr[], const int rowind[], int base, const doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy, int *info);
void _hsr_zusmv(char uplo, int n, doublecomplex alpha, const doublecomplex val[], const int rowptr[], const int colind[], int base, const doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy, int *info);
void _ssc_zusmv(char uplo, int n, doublecomplex alpha, const doublecomplex val[], const int colptr[], const int rowind[], int base, const doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy, int *info);
void _ssr_zusmv(char uplo, int n, doublecomplex alpha, const doublecomplex val[], const int rowptr[], const int colind[], int base, const doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy, int *info);

void _csc_zussv(char uplo, char trans, char diag, int n, const doublecomplex val[], const int colptr[], const int rowind[], int base, doublecomplex x[], int incx, int *info);
void _csr_zussv(char uplo, char trans, char diag, int n, const doublecomplex val[], const int rowptr[], const int colind[], int base, doublecomplex x[], int incx, int *info);

void _csc_zusmm(char trans, char order, int m, int n, int nrhs, doublecomplex alpha, const doublecomplex val[], const int colptr[], const int rowind[], int base, int ldb, const doublecomplex b[], doublecomplex beta, int ldc, doublecomplex c[], int *info);
void _csr_zusmm(char trans, char order, int m, int n, int nrhs, doublecomplex alpha, const doublecomplex val[], const int rowptr[], const int colind[], int base, int ldb, const doublecomplex b[], doublecomplex beta, int ldc, doublecomplex c[], int *info);
void _csc_zussm(char uplo, char trans, char diag, char order, int n, int nrhs, const doublecomplex val[], const int colptr[], const int rowind[], int base, int ldx, doublecomplex x[], int *info);
void _csr_zussm(char uplo, char trans, char diag, char order, int n, int nrhs, const doublecomplex val[], const int rowptr[], const int colind[], int base, int ldx, doublecomplex x[], int *info);

/*
 * SpUtils
 */
void _coo_csc(int m, int n, int nnz, const double val[], const int rowind[], const int colind[], int base, double val2[], int colptr2[], int rowind2[], int base2, int *info);
void _coo_csr(int m, int n, int nnz, const double val[], const int rowind[], const int colind[], int base, double val2[], int rowptr2[], int colind2[], int base2, int *info);
void _csc_coo(int m, int n, const double val[], const int colptr[], const int rowind[], int base, double val2[], int rowind2[], int colind2[], int base2, int *info);
void _csr_coo(int m, int n, const double val[], const int rowptr[], const int colind[], int base, double val2[], int rowind2[], int colind2[], int base2, int *info);

void _csc_csr(int m, int n, const double val[], const int ptr[], const int ind[], int base, double val2[], int ptr2[], int ind2[], int base2, int *info);

void _csc_ssc(char uplo, int n, const double val[], const int colptr[], const int rowind[], int base, int maxnnz2, double val2[], int colptr2[], int rowind2[], int base2, int *info);
void _csr_ssr(char uplo, int n, const double val[], const int rowptr[], const int colind[], int base, int maxnnz2, double val2[], int rowptr2[], int colind2[], int base2, int *info);
void _ssc_csc(char uplo, int n, const double val[], const int colptr[], const int rowind[], int base, int maxnnz2, double val2[], int colptr2[], int rowind2[], int base2, int *info);
void _ssr_csr(char uplo, int n, const double val[], const int rowptr[], const int colind[], int base, int maxnnz2, double val2[], int rowptr2[], int colind2[], int base2, int *info);

void _csc_dense(int m, int n, const double val[], const int colptr[], const int rowind[], int base, int lda, double a[], int *info);
void _csr_dense(int m, int n, const double val[], const int rowptr[], const int colind[], int base, int lda, double a[], int *info);
void _coo_dense(int m, int n, int nnz, const double val[], const int rowind[], const int colind[], int base, int lda, double a[], int *info);
void _dense_csc(int m, int n, int lda, const double a[], int maxnnz, double val[], int colptr[], int rowind[], int base, int *info);
void _dense_csr(int m, int n, int lda, const double a[], int maxnnz, double val[], int rowptr[], int colind[], int base, int *info);
void _dense_coo(int m, int n, int lda, const double a[], int *nnz, double val[], int rowind[], int colind[], int base, int *info);

void _csx_diag(int m, int n, const double val[], const int ptr[], const int ind[], int base, double alpha, int inv, double d[], int incd, int *info);
void _csx_sort(int n, double val[], const int ptr[], int ind[], int base, int *info);
void _dense_nnz(int m, int n, int lda, const double a[], int *nnz, int *info);

void _z_coo_csc(int m, int n, int nnz, const doublecomplex val[], const int rowind[], const int colind[], int base, doublecomplex val2[], int colptr2[], int rowind2[], int base2, int *info);
void _z_coo_csr(int m, int n, int nnz, const doublecomplex val[], const int rowind[], const int colind[], int base, doublecomplex val2[], int rowptr2[], int colind2[], int base2, int *info);
void _z_csc_coo(int m, int n, const doublecomplex val[], const int colptr[], const int rowind[], int base, doublecomplex val2[], int rowind2[], int colind2[], int base2, int *info);
void _z_csr_coo(int m, int n, const doublecomplex val[], const int rowptr[], const int colind[], int base, doublecomplex val2[], int rowind2[], int colind2[], int base2, int *info);

void _z_csc_csr(int m, int n, const doublecomplex val[], const int ptr[], const int ind[], int base, doublecomplex val2[], int ptr2[], int ind2[], int base2, int conjg, int *info);

void _z_csc_ssc(char uplo, int n, const doublecomplex val[], const int colptr[], const int rowind[], int base, int maxnnz2, doublecomplex val2[], int colptr2[], int rowind2[], int base2, int *info);
void _z_csr_ssr(char uplo, int n, const doublecomplex val[], const int rowptr[], const int colind[], int base, int maxnnz2, doublecomplex val2[], int rowptr2[], int colind2[], int base2, int *info);
void _z_ssc_csc(char uplo, int n, const doublecomplex val[], const int colptr[], const int rowind[], int base, int maxnnz2, doublecomplex val2[], int colptr2[], int rowind2[], int base2, int *info);
void _z_ssr_csr(char uplo, int n, const doublecomplex val[], const int rowptr[], const int colind[], int base, int maxnnz2, doublecomplex val2[], int rowptr2[], int colind2[], int base2, int *info);
void _z_hsc_csc(char uplo, int n, const doublecomplex val[], const int colptr[], const int rowind[], int base, int maxnnz2, doublecomplex val2[], int colptr2[], int rowind2[], int base2, int *info);
void _z_hsr_csr(char uplo, int n, const doublecomplex val[], const int rowptr[], const int colind[], int base, int maxnnz2, doublecomplex val2[], int rowptr2[], int colind2[], int base2, int *info);

void _z_csc_dense(int m, int n, const doublecomplex val[], const int colptr[], const int rowind[], int base, int lda, doublecomplex a[], int *info);
void _z_csr_dense(int m, int n, const doublecomplex val[], const int rowptr[], const int colind[], int base, int lda, doublecomplex a[], int *info);
void _z_coo_dense(int m, int n, int nnz, const doublecomplex val[], const int rowind[], const int colind[], int base, int lda, doublecomplex a[], int *info);
void _z_dense_csc(int m, int n, int lda, const doublecomplex a[], int maxnnz, doublecomplex val[], int colptr[], int rowind[], int base, int *info);
void _z_dense_csr(int m, int n, int lda, const doublecomplex a[], int maxnnz, doublecomplex val[], int rowptr[], int colind[], int base, int *info);
void _z_dense_coo(int m, int n, int lda, const doublecomplex a[], int *nnz, doublecomplex val[], int rowind[], int colind[], int base, int *info);

void _z_csx_diag(int m, int n, const doublecomplex val[], const int ptr[], const int ind[], int base, doublecomplex alpha, int inv, doublecomplex d[], int incd, int *info);
void _csx_diag_ind(int m, int n, const int ptr[], const int ind[], int base, int id[], int *info);
void _z_csx_sort(int n, doublecomplex val[], const int ptr[], int ind[], int base, int *info);
void _z_dense_nnz(int m, int n, int lda, const doublecomplex a[], int *nnz, int *info);

void _csc_dussv_sor(char uplo, int n, const double val[], const int colptr[], const int rowind[], int base, double omega, double x[], int incx, int *info);
void _csr_dussv_sor(char uplo, int n, const double val[], const int rowptr[], const int colind[], int base, double omega, double x[], int incx, int *info);
void _csc_zussv_sor(char uplo, int n, const doublecomplex val[], const int colptr[], const int rowind[], int base, double omega, doublecomplex x[], int incx, int *info);
void _csr_zussv_sor(char uplo, int n, const doublecomplex val[], const int rowptr[], const int colind[], int base, double omega, doublecomplex x[], int incx, int *info);

void _csc_dusadd(int m, int n, double alpha, const double val_a[], const int colptr_a[], const int rowind_a[], int base_a, double beta, const double val_b[], const int colptr_b[], const int rowind_b[], int base_b, double val_c[], int colptr_c[], int rowind_c[], int base_c, int nnz_c, int *info);
void _csr_dusadd(int m, int n, double alpha, const double val_a[], const int rowptr_a[], const int colind_a[], int base_a, double beta, const double val_b[], const int rowptr_b[], const int colind_b[], int base_b, double val_c[], int rowptr_c[], int colind_c[], int base_c, int nnz_c, int *info);
void _csc_dzusadd(int m, int n, doublecomplex alpha, const double val_a[], const int colptr_a[], const int rowind_a[], int base_a, double beta, const double val_b[], const int colptr_b[], const int rowind_b[], int base_b, doublecomplex val_c[], int colptr_c[], int rowind_c[], int base_c, int nnz_c, int *info);
void _csr_dzusadd(int m, int n, doublecomplex alpha, const double val_a[], const int rowptr_a[], const int colind_a[], int base_a, double beta, const double val_b[], const int rowptr_b[], const int colind_b[], int base_b, doublecomplex val_c[], int rowptr_c[], int colind_c[], int base_c, int nnz_c, int *info);
void _csc_zusadd(int m, int n, doublecomplex alpha, const doublecomplex val_a[], const int colptr_a[], const int rowind_a[], int base_a, doublecomplex beta, const doublecomplex val_b[], const int colptr_b[], const int rowind_b[], int base_b, doublecomplex val_c[], int colptr_c[], int rowind_c[], int base_c, int nnz_c, int *info);
void _csr_zusadd(int m, int n, doublecomplex alpha, const doublecomplex val_a[], const int rowptr_a[], const int colind_a[], int base_a, doublecomplex beta, const doublecomplex val_b[], const int rowptr_b[], const int colind_b[], int base_b, doublecomplex val_c[], int rowptr_c[], int colind_c[], int base_c, int nnz_c, int *info);

/*
 * Check matrix routines
 */
void _csc_check(int m, int n, const double val[], const int ptr[], const int ind[], int result[], int *info);
void _csr_check(int m, int n, const double val[], const int ptr[], const int ind[], int result[], int *info);
void _csx_check_sym(int n, const double val[], const int ptr[], const int ind[], int *info);
void _z_csc_check(int m, int n, const doublecomplex val[], const int ptr[], const int ind[], int result[], int *info);
void _z_csr_check(int m, int n, const doublecomplex val[], const int ptr[], const int ind[], int result[], int *info);
void _z_csx_check_sym(int n, const doublecomplex val[], const int ptr[], const int ind[], int *info);

/*
 * Input/output
 */
void _hb_read(char *fname, char *title, char *key, char *mxtype, int *nrow, int *ncol, int *nnz, int *neltvl, double val[], int lval, int ptr[], int lptr, int ind[], int lind, char *rhstyp, int *nrhs, int *nrhsix, double rhsval[], int lrhsval, int rhsptr[], int lrhsptr, int rhsind[], int lrhsind, double sguess[], int lsguess, double xexact[], int lxexact, int skip, int base, int format, int sort, int *info);
void _hb_write(char *fname, char *title, char *key, char *mxtype, int nrow, int ncol, int nnz, int neltvl, double val[], int lval, int ptr[], int lptr, int ind[], int lind, char *rhstyp, int nrhs, int nrhsix, double rhsval[], int lrhsval, int rhsptr[], int lrhsptr, int rhsind[], int lrhsind, double sguess[], int lsguess, double xexact[], int lxexact, int base, int fchk, int format, int int_n, int int_w, int dbl_n, int dbl_w, int dbl_d, int *info);
void _mm_read(char *fname, char *matcode, int *nrow, int *ncol, int *nnz, double val[], int lval, int ptr[], int lptr, int ind[], int lind, int skip, int base, int format, int sort, int *info);
void _mm_write(char *fname, char *matcode, int nrow, int ncol, int nnz, double val[], int lval, int ptr[], int lptr, int ind[], int lind, int base, int fchk, int format, int *info);

/*
 * SuerLU
 */
void _dgssv(char col_perm, double thresh, char sym_mode, int n, const double val[], const int ptr[], const int ind[], int format, int base, int perm_c[], int perm_r[], int nrhs, int ldb, double b[], double *rcond, int *info);
void _dgssvx(char equil, char col_perm, char trans, char refine, double thresh, char sym_mode, int n, double val[], const int ptr[], const int ind[], int format, int base, int perm_c[], int perm_r[], int etree[], char *equed, double r[], double c[], int nrhs, int ldb, double b[], int ldx, double x[], double *rpg, double *rcond, double ferr[], double berr[], int *info);
void _zgssv(char col_perm, double thresh, char sym_mode, int n, const doublecomplex val[], const int ptr[], const int ind[], int format, int base, int perm_c[], int perm_r[], int nrhs, int ldb, doublecomplex b[], double *rcond, int *info);
void _zgssvx(char equil, char col_perm, char trans, char refine, double thresh, char sym_mode, int n, doublecomplex val[], const int ptr[], const int ind[], int format, int base, int perm_c[], int perm_r[], int etree[], char *equed, double r[], double c[], int nrhs, int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rpg, double *rcond, double ferr[], double berr[], int *info);

/*
 * Arpack
 */
void _dssev(char jobz, char uplo, int n, const double val[], const int ptr[], const int ind[], int base, int format, double d[], int ldz, double z[], const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int maxiter, double workd[], double workl[], int lworkl, int iwork[], int *info);
void _dssevs(char jobz, char uplo, int n, const double val[], const int ptr[], const int ind[], int base, int format, double sigma, double d[], int ldz, double z[], const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int maxiter, double workd[], double workl[], int lworkl, int iwork[], int *info);
void _dgsev(char jobz, int n, const double val[], const int ptr[], const int ind[], int base, int format, double dr[], double di[], int ldz, double z[], const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int maxiter, double workd[], double workl[], int lworkl, double workev[], int iwork[], int *info);
void _dgsevs(char jobz, int n, const double val[], const int ptr[], const int ind[], int base, int format, double sigmar, double sigmai, double dr[], double di[], int ldz, double z[], const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int maxiter, double workd[], double workl[], int lworkl, double workev[], int iwork[], int *info);
void _dssgvs(char jobz, char uplo, char uplo2, int n, const double val[], const int ptr[], const int ind[], int base, const double val2[], const int ptr2[], const int ind2[], int base2, int format, double sigma, double d[], int ldz, double z[], const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int maxiter, double workd[], double workl[], int lworkl, int iwork[], int *info);
void _dssgv(char jobz, char uplo, char uplo2, int n, const double val[], const int ptr[], const int ind[], int base, const double val2[], const int ptr2[], const int ind2[], int base2, int format, double d[], int ldz, double z[], const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int maxiter, double workd[], double workl[], int lworkl, int iwork[], int *info);
void _dgsgv(char jobz, char uplo2, int n, const double val[], const int ptr[], const int ind[], int base, const double val2[], const int ptr2[], const int ind2[], int base2, int format, double dr[], double di[], int ldz, double z[], const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int maxiter, double workd[], double workl[], int lworkl, double workev[], int iwork[], int *info);
void _dgsgvs(char jobz, char uplo2, int n, const double val[], const int ptr[], const int ind[], int base, const double val2[], const int ptr2[], const int ind2[], int base2, int format, double sigmar, double sigmai, double dr[], double di[], int ldz, double z[], const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int maxiter, int mode, double workd[], double workl[], int lworkl, double workev[], int iwork[], int *info);
void _zgsev(char jobz, int n, const doublecomplex val[], const int ptr[], const int ind[], int base, int format, doublecomplex d[], int ldz, doublecomplex z[], const char *which, int nev, double tol, doublecomplex resid[], int ncv, int ldv, doublecomplex v[], int maxiter, doublecomplex workd[], doublecomplex workl[], int lworkl, doublecomplex workev[], double rwork[], int iwork[], int *info);
void _zgsevs(char jobz, int n, const doublecomplex val[], const int ptr[], const int ind[], int base, int format, doublecomplex sigma, doublecomplex d[], int ldz, doublecomplex z[], const char *which, int nev, double tol, doublecomplex resid[], int ncv, int ldv, doublecomplex v[], int maxiter, doublecomplex workd[], doublecomplex workl[], int lworkl, doublecomplex workev[], double rwork[], int iwork[], int *info);
void _zgsgv(char jobz, char uplo2, int n, const doublecomplex val[], const int ptr[], const int ind[], int base, const doublecomplex val2[], const int ptr2[], const int ind2[], int base2, int format, doublecomplex d[], int ldz, doublecomplex z[], const char *which, int nev, double tol, doublecomplex resid[], int ncv, int ldv, doublecomplex v[], int maxiter, doublecomplex workd[], doublecomplex workl[], int lworkl, doublecomplex workev[], double rwork[], int iwork[], int *info);
void _zgsgvs(char jobz, char uplo2, int n, const doublecomplex val[], const int ptr[], const int ind[], int base, const doublecomplex val2[], const int ptr2[], const int ind2[], int base2, int format, doublecomplex sigma, doublecomplex d[], int ldz, doublecomplex z[], const char *which, int nev, double tol, doublecomplex resid[], int ncv, int ldv, doublecomplex v[], int maxiter, doublecomplex workd[], doublecomplex workl[], int lworkl, doublecomplex workev[], double rwork[], int iwork[], int *info);

void _dgssvd(char jobv, char jobu, int m, int n, const double val[], const int ptr[], const int ind[], int base, int format, double s[], const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int ldu, double u[], int maxiter, double workd[], double workl[], int lworkl, double workv[], int iwork[], int *info);
void _zgssvd(char jobv, char jobu, int m, int n, const doublecomplex val[], const int ptr[], const int ind[], int base, int format, doublecomplex s[], const char *which, int nev, double tol, doublecomplex resid[], int ncv, int ldv, doublecomplex v[], int ldu, doublecomplex u[], int maxiter, doublecomplex workd[], doublecomplex workl[], int lworkl, doublecomplex workev[], doublecomplex workv[], double rwork[], int iwork[], int *info);

void _dnaupd(int *ido, char bmat, int n, const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int iparam[], int ipntr[], double workd[], double workl[], int lworkl, int *info);
void _dneupd(int rvec, char howmny, int select[], double dr[], double di[], int ldz, double z[], double sigmar, double sigmai, double workev[], char bmat, int n, const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int iparam[], int ipntr[], double workd[], double workl[], int lworkl, int *info);
void _dsaupd(int *ido, char bmat, int n, const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int iparam[], int ipntr[], double workd[], double workl[], int lworkl, int *info);
void _dseupd(int rvec, char howmny, int select[], double d[], int ldz, double z[], double sigma, char bmat, int n, const char *which, int nev, double tol, double resid[], int ncv, int ldv, double v[], int iparam[], int ipntr[], double workd[], double workl[], int lworkl, int *info);
void _znaupd(int* ido, char bmat, int n, char const* which, int nev, double tol, doublecomplex resid[], int ncv, int ldv, doublecomplex v[], int iparam[], int ipntr[], doublecomplex workd[], doublecomplex workl[], int lworkl, double rwork[], int* info);
void _zneupd(int rvec, char howmny, int select[], doublecomplex d[], int ldz, doublecomplex z[], doublecomplex sigma, doublecomplex workev[], char bmat, int n, const char *which, int nev, double tol, doublecomplex resid[], int ncv, int ldv, doublecomplex v[], int iparam[], int ipntr[], doublecomplex workd[], doublecomplex workl[], int lworkl, double rwork[], int* info);

#if defined(__cplusplus)
}
#endif
