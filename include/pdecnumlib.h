/****************************************
 *                                      *
 *  XLPack Numerical Library            *
 *  Version 7.0 (January 31, 2023)      *
 *  (C) 2014-2023  K Technologies       *
 *                                      *
 ****************************************/
#pragma once

#if !defined(_PDECNUMLIB_NO_MANGLING)
#include "pdecnumlib_mangling.h"
#endif

#if defined(__cplusplus)
extern "C" {
#endif

void fem2p(int n, int ne, double x[], double y[], int ldknc, int knc[], double p[], double q[], double f[], int nb1, int ib[], double bv[], int nb2, int ldks2, int ks2[], int ldalpha, double alpha[], int ldbeta, double beta[], int asize, double val[], int rowptr[], int colind[], int base, double b[], int iwork[], int *info);
void _fem2pd(int n, int ne, double x[], double y[], int ldknc, int knc[], double p[], double q[], double f[], int nb1, int ib[], double bv[], int nb2, int ldks2, int ks2[], int ldalpha, double alpha[], int ldbeta, double beta[], int lda, double a[], double b[], int *info);
void fem3p(int n, int ne, double x[], double y[], double z[], int ldknc, int knc[], double p[], double q[], double f[], int nb1, int ib[], double bv[], int nb2, int ldks2, int ks2[], int ldalpha, double alpha[], int ldbeta, double beta[], int asize, double val[], int rowptr[], int colind[], int base, double b[], int iwork[], int *info);
void _fem3pd(int n, int ne, double x[], double y[], double z[], int ldknc, int knc[], double p[], double q[], double f[], int nb1, int ib[], double bv[], int nb2, int ldks2, int ks2[], int ldalpha, double alpha[], int ldbeta, double beta[], int lda, double a[], double b[], int *info);

void _mesh23(int nx, int ny, int *n, double sclx, double scly, double x[], double y[], int *ne, int ldknc, int knc[], int ldks, int ks[], int lb[], int *nb, int *info);
void _mesh24(int nx, int ny, int *n, double sclx, double scly, double x[], double y[], int *ne, int ldknc, int knc[], int ldks, int ks[], int lb[], int *nb, int *info);
void _mesh29(int nx, int ny, int *n, double sclx, double scly, double x[], double y[], int *ne, int ldknc, int knc[], int ldks, int ks[], int lb[], int *nb, int *info);

void _mesh34(int nx, int ny, int nz, int *n, double sclx, double scly, double sclz, double x[], double y[], double z[], int *ne, int ldknc, int knc[], int ldks, int ks[], int lb[], int *nb, int *info);
void _mesh35(int nx, int ny, int nz, int *n, double sclx, double scly, double sclz, double x[], double y[], double z[], int *ne, int ldknc, int knc[], int ldks, int ks[], int lb[], int *nb, int *info);
void _mesh36(int nx, int ny, int nz, int *n, double sclx, double scly, double sclz, double x[], double y[], double z[], int *ne, int ldknc, int knc[], int ldks, int ks[], int lb[], int *nb, int *info);
void _mesh38(int nx, int ny, int nz, int *n, double sclx, double scly, double sclz, double x[], double y[], double z[], int *ne, int ldknc, int knc[], int ldks, int ks[], int lb[], int *nb, int *info);

void _readgmsh22(char *fname, int *n, double x[], double y[], double z[], int *ne, int ldkc, int kc[], int ldlb, int lb[], int *info);
void _readmsh2(char *fname, int *n, double x[], double y[], int ln[], int *ne, int ldkc, int kc[], int le[], int *nb, int ldks, int ks[], int lb[], int *info);
void _readmsh3(char *fname, int *n, double x[], double y[], double z[], int ln[], int *ne, int ldkc, int kc[], int le[], int *nb, int ldks, int ks[], int lb[], int *info);
void _readmesh(char *fname, int *ndim, int *n, double x[], double y[], double z[], int ln[], int *ne, int ldkc, int kc[], int le[], int *nb, int ldks, int ks[], int lb[], int *info);
void _writecsv2(char *fname, int n, double x[], double y[], double u[], int *info);
void _writecsv3(char *fname, int n, double x[], double y[], double z[], double u[], int *info);
void _writegmsh22(char *fname, int n, double x[], double y[], double z[], int ne, int ldkc, int kc[], int ldlb, int lb[], int *info);
void _writevtkug(char *fname, int n, double x[], double y[], double z[], int ne, int ldkc, int kc[], double u[], int *info);

#if defined(__cplusplus)
}
#endif
