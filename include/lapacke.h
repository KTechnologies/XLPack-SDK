/*****************************************************************************
  Copyright (c) 2014, Intel Corp.
  All rights reserved.

  Redistribution and use in source and binary forms, with or without
  modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
  THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************
* Contents: Native C interface to LAPACK
* Author: Intel Corporation
* Generated August, 2015
*****************************************************************************/

/* Modified for XLPack 6.0, Feb. 12, 2022, K Technologies */

#ifndef _LAPACKE_H_
#define _LAPACKE_H_

/*
*  Turn on HAVE_LAPACK_CONFIG_H to redefine C-LAPACK datatypes
*/
#ifdef HAVE_LAPACK_CONFIG_H
#include "lapacke_config.h"
#endif

#include <stdlib.h>

#ifndef lapack_int
#define lapack_int     int
#endif

#ifndef lapack_logical
#define lapack_logical lapack_int
#endif

/* Complex types are structures equivalent to the
* Fortran complex types COMPLEX(4) and COMPLEX(8).
*
* One can also redefine the types with his own types
* for example by including in the code definitions like
*
* #define lapack_complex_float std::complex<float>
* #define lapack_complex_double std::complex<double>
*
* or define these types in the command line:
*
* -Dlapack_complex_float="std::complex<float>"
* -Dlapack_complex_double="std::complex<double>"
*/

#ifndef LAPACK_COMPLEX_CUSTOM

/* Complex type (single precision) */
#ifndef lapack_complex_float

#if defined(__cplusplus)
#include <cstddef>
#include <complex>
#define lapack_complex_float	std::complex<float>
#else
#include <stddef.h>
#include <complex.h>
#if defined(__clang__) || defined(__GNUC__) || defined(__APPLE__)
#define lapack_complex_float	float _Complex
#if !defined(_NO_VLARRAY)
#define _VLARRAY
#endif
#else	// Visual C/C++
#define lapack_complex_float	_Fcomplex
#endif
#endif

#endif

#ifndef lapack_complex_float_real
#define lapack_complex_float_real(z)	(creal(z))
#endif

#ifndef lapack_complex_float_imag
#define lapack_complex_float_imag(z)	(cimag(z))
#endif

lapack_complex_float lapack_make_complex_float( float re, float im );

/* Complex type (double precision) */
#ifndef lapack_complex_double

#if defined(__cplusplus)
#include <cstddef>
#include <complex>
#define lapack_complex_double	std::complex<double>
#else
#include <stddef.h>
#include <complex.h>
#if defined(__clang__) || defined(__GNUC__) || defined(__APPLE__)
#define lapack_complex_double	double _Complex
#if !defined(_NO_VLARRAY)
#define _VLARRAY
#endif
#else	// Visual C/C++
#define lapack_complex_double	_Dcomplex
#endif
#endif

#endif

#ifndef lapack_complex_double_real
#define lapack_complex_double_real(z)	(creal(z))
#endif

#ifndef lapack_complex_double_imag
#define lapack_complex_double_imag(z)	(cimag(z))
#endif

lapack_complex_double lapack_make_complex_double( double re, double im );

#endif	// LAPACK_COMPLEX_CUSTOM

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef LAPACKE_malloc
#define LAPACKE_malloc( size ) malloc( size )
#endif
#ifndef LAPACKE_free
#define LAPACKE_free( p )      free( p )
#endif

#define LAPACK_C2INT( x ) (lapack_int)(*((float*)&x ))
#define LAPACK_Z2INT( x ) (lapack_int)(*((double*)&x ))

#define LAPACK_ROW_MAJOR               101
#define LAPACK_COL_MAJOR               102

#define LAPACK_WORK_MEMORY_ERROR       -1010
#define LAPACK_TRANSPOSE_MEMORY_ERROR  -1011

/* Callback logical functions of one, two, or three arguments are used
*  to select eigenvalues to sort to the top left of the Schur form.
*  The value is selected if function returns TRUE (non-zero). */

typedef lapack_logical (*LAPACK_D_SELECT2) ( const double*, const double* );
typedef lapack_logical (*LAPACK_D_SELECT3)
    ( const double*, const double*, const double* );

typedef lapack_logical (*LAPACK_Z_SELECT1) ( const lapack_complex_double* );
typedef lapack_logical (*LAPACK_Z_SELECT2)
    ( const lapack_complex_double*, const lapack_complex_double* );

/* LAPACKE function prototypes */

lapack_int LAPACKE_ddisna( char job, lapack_int m, lapack_int n,
                           const double* d, double* sep );


lapack_int LAPACKE_dgbcon( int matrix_layout, char norm, lapack_int n,
                           lapack_int kl, lapack_int ku, const double* ab,
                           lapack_int ldab, const lapack_int* ipiv,
                           double anorm, double* rcond );
lapack_int LAPACKE_zgbcon( int matrix_layout, char norm, lapack_int n,
                           lapack_int kl, lapack_int ku,
                           const lapack_complex_double* ab, lapack_int ldab,
                           const lapack_int* ipiv, double anorm,
                           double* rcond );

lapack_int LAPACKE_dgbsv( int matrix_layout, lapack_int n, lapack_int kl,
                          lapack_int ku, lapack_int nrhs, double* ab,
                          lapack_int ldab, lapack_int* ipiv, double* b,
                          lapack_int ldb );
lapack_int LAPACKE_zgbsv( int matrix_layout, lapack_int n, lapack_int kl,
                          lapack_int ku, lapack_int nrhs,
                          lapack_complex_double* ab, lapack_int ldab,
                          lapack_int* ipiv, lapack_complex_double* b,
                          lapack_int ldb );

lapack_int LAPACKE_dgbsvx( int matrix_layout, char fact, char trans,
                           lapack_int n, lapack_int kl, lapack_int ku,
                           lapack_int nrhs, double* ab, lapack_int ldab,
                           double* afb, lapack_int ldafb, lapack_int* ipiv,
                           char* equed, double* r, double* c, double* b,
                           lapack_int ldb, double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr,
                           double* rpivot );
lapack_int LAPACKE_zgbsvx( int matrix_layout, char fact, char trans,
                           lapack_int n, lapack_int kl, lapack_int ku,
                           lapack_int nrhs, lapack_complex_double* ab,
                           lapack_int ldab, lapack_complex_double* afb,
                           lapack_int ldafb, lapack_int* ipiv, char* equed,
                           double* r, double* c, lapack_complex_double* b,
                           lapack_int ldb, lapack_complex_double* x,
                           lapack_int ldx, double* rcond, double* ferr,
                           double* berr, double* rpivot );

lapack_int LAPACKE_dgbtrf( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int kl, lapack_int ku, double* ab,
                           lapack_int ldab, lapack_int* ipiv );
lapack_int LAPACKE_zgbtrf( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int kl, lapack_int ku,
                           lapack_complex_double* ab, lapack_int ldab,
                           lapack_int* ipiv );

lapack_int LAPACKE_dgbtrs( int matrix_layout, char trans, lapack_int n,
                           lapack_int kl, lapack_int ku, lapack_int nrhs,
                           const double* ab, lapack_int ldab,
                           const lapack_int* ipiv, double* b, lapack_int ldb );
lapack_int LAPACKE_zgbtrs( int matrix_layout, char trans, lapack_int n,
                           lapack_int kl, lapack_int ku, lapack_int nrhs,
                           const lapack_complex_double* ab, lapack_int ldab,
                           const lapack_int* ipiv, lapack_complex_double* b,
                           lapack_int ldb );

lapack_int LAPACKE_dgebak( int matrix_layout, char job, char side, lapack_int n,
                           lapack_int ilo, lapack_int ihi, const double* scale,
                           lapack_int m, double* v, lapack_int ldv );
lapack_int LAPACKE_zgebak( int matrix_layout, char job, char side, lapack_int n,
                           lapack_int ilo, lapack_int ihi, const double* scale,
                           lapack_int m, lapack_complex_double* v,
                           lapack_int ldv );

lapack_int LAPACKE_dgebal( int matrix_layout, char job, lapack_int n, double* a,
                           lapack_int lda, lapack_int* ilo, lapack_int* ihi,
                           double* scale );
lapack_int LAPACKE_zgebal( int matrix_layout, char job, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_int* ilo, lapack_int* ihi, double* scale );

lapack_int LAPACKE_dgecon( int matrix_layout, char norm, lapack_int n,
                           const double* a, lapack_int lda, double anorm,
                           double* rcond );
lapack_int LAPACKE_zgecon( int matrix_layout, char norm, lapack_int n,
                           const lapack_complex_double* a, lapack_int lda,
                           double anorm, double* rcond );

lapack_int LAPACKE_dgees( int matrix_layout, char jobvs, char sort,
                          LAPACK_D_SELECT2 select, lapack_int n, double* a,
                          lapack_int lda, lapack_int* sdim, double* wr,
                          double* wi, double* vs, lapack_int ldvs );
lapack_int LAPACKE_zgees( int matrix_layout, char jobvs, char sort,
                          LAPACK_Z_SELECT1 select, lapack_int n,
                          lapack_complex_double* a, lapack_int lda,
                          lapack_int* sdim, lapack_complex_double* w,
                          lapack_complex_double* vs, lapack_int ldvs );

lapack_int LAPACKE_dgeesx( int matrix_layout, char jobvs, char sort,
                           LAPACK_D_SELECT2 select, char sense, lapack_int n,
                           double* a, lapack_int lda, lapack_int* sdim,
                           double* wr, double* wi, double* vs, lapack_int ldvs,
                           double* rconde, double* rcondv );
lapack_int LAPACKE_zgeesx( int matrix_layout, char jobvs, char sort,
                           LAPACK_Z_SELECT1 select, char sense, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_int* sdim, lapack_complex_double* w,
                           lapack_complex_double* vs, lapack_int ldvs,
                           double* rconde, double* rcondv );

lapack_int LAPACKE_dgeev( int matrix_layout, char jobvl, char jobvr,
                          lapack_int n, double* a, lapack_int lda, double* wr,
                          double* wi, double* vl, lapack_int ldvl, double* vr,
                          lapack_int ldvr );
lapack_int LAPACKE_zgeev( int matrix_layout, char jobvl, char jobvr,
                          lapack_int n, lapack_complex_double* a,
                          lapack_int lda, lapack_complex_double* w,
                          lapack_complex_double* vl, lapack_int ldvl,
                          lapack_complex_double* vr, lapack_int ldvr );

lapack_int LAPACKE_dgeevx( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, lapack_int n, double* a,
                           lapack_int lda, double* wr, double* wi, double* vl,
                           lapack_int ldvl, double* vr, lapack_int ldvr,
                           lapack_int* ilo, lapack_int* ihi, double* scale,
                           double* abnrm, double* rconde, double* rcondv );
lapack_int LAPACKE_zgeevx( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_complex_double* w, lapack_complex_double* vl,
                           lapack_int ldvl, lapack_complex_double* vr,
                           lapack_int ldvr, lapack_int* ilo, lapack_int* ihi,
                           double* scale, double* abnrm, double* rconde,
                           double* rcondv );

lapack_int LAPACKE_dgehrd( int matrix_layout, lapack_int n, lapack_int ilo,
                           lapack_int ihi, double* a, lapack_int lda,
                           double* tau );
lapack_int LAPACKE_zgehrd( int matrix_layout, lapack_int n, lapack_int ilo,
                           lapack_int ihi, lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* tau );

lapack_int LAPACKE_dgejsv( int matrix_layout, char joba, char jobu, char jobv,
                           char jobr, char jobt, char jobp, lapack_int m,
                           lapack_int n, double* a, lapack_int lda, double* sva,
                           double* u, lapack_int ldu, double* v, lapack_int ldv,
                           double* stat, lapack_int* istat );
lapack_int LAPACKE_zgejsv( int matrix_layout, char joba, char jobu, char jobv,
                           char jobr, char jobt, char jobp, lapack_int m,
                           lapack_int n, lapack_complex_double* a, lapack_int lda, double* sva,
                           lapack_complex_double* u, lapack_int ldu, lapack_complex_double* v, lapack_int ldv,
                           double* stat, lapack_int* istat );

lapack_int LAPACKE_dgelqf( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, double* tau );
lapack_int LAPACKE_zgelqf( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_complex_double* tau );

lapack_int LAPACKE_dgels( int matrix_layout, char trans, lapack_int m,
                          lapack_int n, lapack_int nrhs, double* a,
                          lapack_int lda, double* b, lapack_int ldb );
lapack_int LAPACKE_zgels( int matrix_layout, char trans, lapack_int m,
                          lapack_int n, lapack_int nrhs,
                          lapack_complex_double* a, lapack_int lda,
                          lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dgelsd( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int nrhs, double* a, lapack_int lda,
                           double* b, lapack_int ldb, double* s, double rcond,
                           lapack_int* rank );
lapack_int LAPACKE_zgelsd( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int nrhs, lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* b,
                           lapack_int ldb, double* s, double rcond,
                           lapack_int* rank );

lapack_int LAPACKE_dgelss( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int nrhs, double* a, lapack_int lda,
                           double* b, lapack_int ldb, double* s, double rcond,
                           lapack_int* rank );
lapack_int LAPACKE_zgelss( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int nrhs, lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* b,
                           lapack_int ldb, double* s, double rcond,
                           lapack_int* rank );

lapack_int LAPACKE_dgelsy( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int nrhs, double* a, lapack_int lda,
                           double* b, lapack_int ldb, lapack_int* jpvt,
                           double rcond, lapack_int* rank );
lapack_int LAPACKE_zgelsy( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int nrhs, lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* b,
                           lapack_int ldb, lapack_int* jpvt, double rcond,
                           lapack_int* rank );

lapack_int LAPACKE_dgeqp3( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, lapack_int* jpvt,
                           double* tau );
lapack_int LAPACKE_zgeqp3( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_int* jpvt, lapack_complex_double* tau );

lapack_int LAPACKE_dgeqrf( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, double* tau );
lapack_int LAPACKE_zgeqrf( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_complex_double* tau );

lapack_int LAPACKE_dgesdd( int matrix_layout, char jobz, lapack_int m,
                           lapack_int n, double* a, lapack_int lda, double* s,
                           double* u, lapack_int ldu, double* vt,
                           lapack_int ldvt );
lapack_int LAPACKE_zgesdd( int matrix_layout, char jobz, lapack_int m,
                           lapack_int n, lapack_complex_double* a,
                           lapack_int lda, double* s, lapack_complex_double* u,
                           lapack_int ldu, lapack_complex_double* vt,
                           lapack_int ldvt );

lapack_int LAPACKE_dgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* a, lapack_int lda, lapack_int* ipiv,
                          double* b, lapack_int ldb );
lapack_int LAPACKE_zgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          lapack_complex_double* a, lapack_int lda,
                          lapack_int* ipiv, lapack_complex_double* b,
                          lapack_int ldb );

lapack_int LAPACKE_dsgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                           double* a, lapack_int lda, lapack_int* ipiv,
                           double* b, lapack_int ldb, double* x, lapack_int ldx,
                           lapack_int* iter );
lapack_int LAPACKE_zcgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_int* ipiv, lapack_complex_double* b,
                           lapack_int ldb, lapack_complex_double* x,
                           lapack_int ldx, lapack_int* iter );

lapack_int LAPACKE_dgesvd( int matrix_layout, char jobu, char jobvt,
                           lapack_int m, lapack_int n, double* a,
                           lapack_int lda, double* s, double* u, lapack_int ldu,
                           double* vt, lapack_int ldvt, double* superb );
lapack_int LAPACKE_zgesvd( int matrix_layout, char jobu, char jobvt,
                           lapack_int m, lapack_int n, lapack_complex_double* a,
                           lapack_int lda, double* s, lapack_complex_double* u,
                           lapack_int ldu, lapack_complex_double* vt,
                           lapack_int ldvt, double* superb );

lapack_int LAPACKE_dgesvdq( int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,
                           lapack_int m, lapack_int n, double* a,
                           lapack_int lda, double* s, double* u, lapack_int ldu,
                           double* v, lapack_int ldv, lapack_int* numrank);
lapack_int LAPACKE_zgesvdq( int matrix_layout, char joba, char jobp, char jobr, char jobu, char jobv,
                           lapack_int m, lapack_int n, lapack_complex_double* a,
                           lapack_int lda, double* s, lapack_complex_double* u,
                           lapack_int ldu, lapack_complex_double* v,
                           lapack_int ldv, lapack_int* numrank );

lapack_int LAPACKE_dgesvdx( int matrix_layout, char jobu, char jobvt, char range,
                           lapack_int m, lapack_int n, double* a,
                           lapack_int lda, double vl, double vu,
                           lapack_int il, lapack_int iu, lapack_int* ns,
                           double* s, double* u, lapack_int ldu,
                           double* vt, lapack_int ldvt,
                           lapack_int* superb );
lapack_int LAPACKE_zgesvdx( int matrix_layout, char jobu, char jobvt, char range,
                           lapack_int m, lapack_int n, lapack_complex_double* a,
                           lapack_int lda, double vl, double vu,
                           lapack_int il, lapack_int iu, lapack_int* ns,
                           double* s, lapack_complex_double* u, lapack_int ldu,
                           lapack_complex_double* vt, lapack_int ldvt,
                           lapack_int* superb );

lapack_int LAPACKE_dgesvx( int matrix_layout, char fact, char trans,
                           lapack_int n, lapack_int nrhs, double* a,
                           lapack_int lda, double* af, lapack_int ldaf,
                           lapack_int* ipiv, char* equed, double* r, double* c,
                           double* b, lapack_int ldb, double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr,
                           double* rpivot );
lapack_int LAPACKE_zgesvx( int matrix_layout, char fact, char trans,
                           lapack_int n, lapack_int nrhs,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_complex_double* af, lapack_int ldaf,
                           lapack_int* ipiv, char* equed, double* r, double* c,
                           lapack_complex_double* b, lapack_int ldb,
                           lapack_complex_double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr,
                           double* rpivot );

lapack_int LAPACKE_dgetrf( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, lapack_int* ipiv );
lapack_int LAPACKE_zgetrf( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_int* ipiv );

lapack_int LAPACKE_dgetri( int matrix_layout, lapack_int n, double* a,
                           lapack_int lda, const lapack_int* ipiv );
lapack_int LAPACKE_zgetri( int matrix_layout, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           const lapack_int* ipiv );

lapack_int LAPACKE_dgetrs( int matrix_layout, char trans, lapack_int n,
                           lapack_int nrhs, const double* a, lapack_int lda,
                           const lapack_int* ipiv, double* b, lapack_int ldb );
lapack_int LAPACKE_zgetrs( int matrix_layout, char trans, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* a,
                           lapack_int lda, const lapack_int* ipiv,
                           lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dgetsls( int matrix_layout, char trans, lapack_int m,
                            lapack_int n, lapack_int nrhs, double* a,
                            lapack_int lda, double* b, lapack_int ldb );
lapack_int LAPACKE_zgetsls( int matrix_layout, char trans, lapack_int m,
                            lapack_int n, lapack_int nrhs,
                            lapack_complex_double* a, lapack_int lda,
                            lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dgges( int matrix_layout, char jobvsl, char jobvsr, char sort,
                          LAPACK_D_SELECT3 selctg, lapack_int n, double* a,
                          lapack_int lda, double* b, lapack_int ldb,
                          lapack_int* sdim, double* alphar, double* alphai,
                          double* beta, double* vsl, lapack_int ldvsl,
                          double* vsr, lapack_int ldvsr );
lapack_int LAPACKE_zgges( int matrix_layout, char jobvsl, char jobvsr, char sort,
                          LAPACK_Z_SELECT2 selctg, lapack_int n,
                          lapack_complex_double* a, lapack_int lda,
                          lapack_complex_double* b, lapack_int ldb,
                          lapack_int* sdim, lapack_complex_double* alpha,
                          lapack_complex_double* beta,
                          lapack_complex_double* vsl, lapack_int ldvsl,
                          lapack_complex_double* vsr, lapack_int ldvsr );

lapack_int LAPACKE_dggesx( int matrix_layout, char jobvsl, char jobvsr,
                           char sort, LAPACK_D_SELECT3 selctg, char sense,
                           lapack_int n, double* a, lapack_int lda, double* b,
                           lapack_int ldb, lapack_int* sdim, double* alphar,
                           double* alphai, double* beta, double* vsl,
                           lapack_int ldvsl, double* vsr, lapack_int ldvsr,
                           double* rconde, double* rcondv );
lapack_int LAPACKE_zggesx( int matrix_layout, char jobvsl, char jobvsr,
                           char sort, LAPACK_Z_SELECT2 selctg, char sense,
                           lapack_int n, lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* b,
                           lapack_int ldb, lapack_int* sdim,
                           lapack_complex_double* alpha,
                           lapack_complex_double* beta,
                           lapack_complex_double* vsl, lapack_int ldvsl,
                           lapack_complex_double* vsr, lapack_int ldvsr,
                           double* rconde, double* rcondv );

lapack_int LAPACKE_dggev( int matrix_layout, char jobvl, char jobvr,
                          lapack_int n, double* a, lapack_int lda, double* b,
                          lapack_int ldb, double* alphar, double* alphai,
                          double* beta, double* vl, lapack_int ldvl, double* vr,
                          lapack_int ldvr );
lapack_int LAPACKE_zggev( int matrix_layout, char jobvl, char jobvr,
                          lapack_int n, lapack_complex_double* a,
                          lapack_int lda, lapack_complex_double* b,
                          lapack_int ldb, lapack_complex_double* alpha,
                          lapack_complex_double* beta,
                          lapack_complex_double* vl, lapack_int ldvl,
                          lapack_complex_double* vr, lapack_int ldvr );

lapack_int LAPACKE_dggevx( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, lapack_int n, double* a,
                           lapack_int lda, double* b, lapack_int ldb,
                           double* alphar, double* alphai, double* beta,
                           double* vl, lapack_int ldvl, double* vr,
                           lapack_int ldvr, lapack_int* ilo, lapack_int* ihi,
                           double* lscale, double* rscale, double* abnrm,
                           double* bbnrm, double* rconde, double* rcondv );
lapack_int LAPACKE_zggevx( int matrix_layout, char balanc, char jobvl,
                           char jobvr, char sense, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_complex_double* b, lapack_int ldb,
                           lapack_complex_double* alpha,
                           lapack_complex_double* beta,
                           lapack_complex_double* vl, lapack_int ldvl,
                           lapack_complex_double* vr, lapack_int ldvr,
                           lapack_int* ilo, lapack_int* ihi, double* lscale,
                           double* rscale, double* abnrm, double* bbnrm,
                           double* rconde, double* rcondv );

lapack_int LAPACKE_dggglm( int matrix_layout, lapack_int n, lapack_int m,
                           lapack_int p, double* a, lapack_int lda, double* b,
                           lapack_int ldb, double* d, double* x, double* y );
lapack_int LAPACKE_zggglm( int matrix_layout, lapack_int n, lapack_int m,
                           lapack_int p, lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* b,
                           lapack_int ldb, lapack_complex_double* d,
                           lapack_complex_double* x, lapack_complex_double* y );

lapack_int LAPACKE_dgglse( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int p, double* a, lapack_int lda, double* b,
                           lapack_int ldb, double* c, double* d, double* x );
lapack_int LAPACKE_zgglse( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int p, lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* b,
                           lapack_int ldb, lapack_complex_double* c,
                           lapack_complex_double* d, lapack_complex_double* x );

lapack_int LAPACKE_dggsvd3( int matrix_layout, char jobu, char jobv, char jobq,
                            lapack_int m, lapack_int n, lapack_int p,
                            lapack_int* k, lapack_int* l, double* a,
                            lapack_int lda, double* b, lapack_int ldb,
                            double* alpha, double* beta, double* u,
                            lapack_int ldu, double* v, lapack_int ldv, double* q,
                            lapack_int ldq, lapack_int* iwork );
lapack_int LAPACKE_zggsvd3( int matrix_layout, char jobu, char jobv, char jobq,
                            lapack_int m, lapack_int n, lapack_int p,
                            lapack_int* k, lapack_int* l,
                            lapack_complex_double* a, lapack_int lda,
                            lapack_complex_double* b, lapack_int ldb,
                            double* alpha, double* beta,
                            lapack_complex_double* u, lapack_int ldu,
                            lapack_complex_double* v, lapack_int ldv,
                            lapack_complex_double* q, lapack_int ldq,
                            lapack_int* iwork );

lapack_int LAPACKE_dgtcon( char norm, lapack_int n, const double* dl,
                           const double* d, const double* du, const double* du2,
                           const lapack_int* ipiv, double anorm,
                           double* rcond );
lapack_int LAPACKE_zgtcon( char norm, lapack_int n,
                           const lapack_complex_double* dl,
                           const lapack_complex_double* d,
                           const lapack_complex_double* du,
                           const lapack_complex_double* du2,
                           const lapack_int* ipiv, double anorm,
                           double* rcond );

lapack_int LAPACKE_dgtsv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* dl, double* d, double* du, double* b,
                          lapack_int ldb );
lapack_int LAPACKE_zgtsv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          lapack_complex_double* dl, lapack_complex_double* d,
                          lapack_complex_double* du, lapack_complex_double* b,
                          lapack_int ldb );

lapack_int LAPACKE_dgtsvx( int matrix_layout, char fact, char trans,
                           lapack_int n, lapack_int nrhs, const double* dl,
                           const double* d, const double* du, double* dlf,
                           double* df, double* duf, double* du2,
                           lapack_int* ipiv, const double* b, lapack_int ldb,
                           double* x, lapack_int ldx, double* rcond,
                           double* ferr, double* berr );
lapack_int LAPACKE_zgtsvx( int matrix_layout, char fact, char trans,
                           lapack_int n, lapack_int nrhs,
                           const lapack_complex_double* dl,
                           const lapack_complex_double* d,
                           const lapack_complex_double* du,
                           lapack_complex_double* dlf,
                           lapack_complex_double* df,
                           lapack_complex_double* duf,
                           lapack_complex_double* du2, lapack_int* ipiv,
                           const lapack_complex_double* b, lapack_int ldb,
                           lapack_complex_double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr );

lapack_int LAPACKE_dgttrf( lapack_int n, double* dl, double* d, double* du,
                           double* du2, lapack_int* ipiv );
lapack_int LAPACKE_zgttrf( lapack_int n, lapack_complex_double* dl,
                           lapack_complex_double* d, lapack_complex_double* du,
                           lapack_complex_double* du2, lapack_int* ipiv );

lapack_int LAPACKE_dgttrs( int matrix_layout, char trans, lapack_int n,
                           lapack_int nrhs, const double* dl, const double* d,
                           const double* du, const double* du2,
                           const lapack_int* ipiv, double* b, lapack_int ldb );
lapack_int LAPACKE_zgttrs( int matrix_layout, char trans, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* dl,
                           const lapack_complex_double* d,
                           const lapack_complex_double* du,
                           const lapack_complex_double* du2,
                           const lapack_int* ipiv, lapack_complex_double* b,
                           lapack_int ldb );

lapack_int LAPACKE_zhbev( int matrix_layout, char jobz, char uplo, lapack_int n,
                          lapack_int kd, lapack_complex_double* ab,
                          lapack_int ldab, double* w, lapack_complex_double* z,
                          lapack_int ldz );

lapack_int LAPACKE_zhbevd( int matrix_layout, char jobz, char uplo, lapack_int n,
                           lapack_int kd, lapack_complex_double* ab,
                           lapack_int ldab, double* w, lapack_complex_double* z,
                           lapack_int ldz );

lapack_int LAPACKE_zhbevx( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, lapack_int kd,
                           lapack_complex_double* ab, lapack_int ldab,
                           lapack_complex_double* q, lapack_int ldq, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w,
                           lapack_complex_double* z, lapack_int ldz,
                           lapack_int* ifail );

lapack_int LAPACKE_zhbgv( int matrix_layout, char jobz, char uplo, lapack_int n,
                          lapack_int ka, lapack_int kb,
                          lapack_complex_double* ab, lapack_int ldab,
                          lapack_complex_double* bb, lapack_int ldbb, double* w,
                          lapack_complex_double* z, lapack_int ldz );

lapack_int LAPACKE_zhbgvd( int matrix_layout, char jobz, char uplo, lapack_int n,
                           lapack_int ka, lapack_int kb,
                           lapack_complex_double* ab, lapack_int ldab,
                           lapack_complex_double* bb, lapack_int ldbb,
                           double* w, lapack_complex_double* z,
                           lapack_int ldz );

lapack_int LAPACKE_zhbgvx( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, lapack_int ka, lapack_int kb,
                           lapack_complex_double* ab, lapack_int ldab,
                           lapack_complex_double* bb, lapack_int ldbb,
                           lapack_complex_double* q, lapack_int ldq, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w,
                           lapack_complex_double* z, lapack_int ldz,
                           lapack_int* ifail );

lapack_int LAPACKE_zhbtrd( int matrix_layout, char vect, char uplo, lapack_int n,
                           lapack_int kd, lapack_complex_double* ab,
                           lapack_int ldab, double* d, double* e,
                           lapack_complex_double* q, lapack_int ldq );

lapack_int LAPACKE_zhecon( int matrix_layout, char uplo, lapack_int n,
                           const lapack_complex_double* a, lapack_int lda,
                           const lapack_int* ipiv, double anorm,
                           double* rcond );

lapack_int LAPACKE_zheev( int matrix_layout, char jobz, char uplo, lapack_int n,
                          lapack_complex_double* a, lapack_int lda, double* w );

lapack_int LAPACKE_zheevd( int matrix_layout, char jobz, char uplo, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           double* w );

lapack_int LAPACKE_zheevr( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, lapack_complex_double* a,
                           lapack_int lda, double vl, double vu, lapack_int il,
                           lapack_int iu, double abstol, lapack_int* m,
                           double* w, lapack_complex_double* z, lapack_int ldz,
                           lapack_int* isuppz );

lapack_int LAPACKE_zheevx( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, lapack_complex_double* a,
                           lapack_int lda, double vl, double vu, lapack_int il,
                           lapack_int iu, double abstol, lapack_int* m,
                           double* w, lapack_complex_double* z, lapack_int ldz,
                           lapack_int* ifail );

lapack_int LAPACKE_zhegv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, lapack_complex_double* a,
                          lapack_int lda, lapack_complex_double* b,
                          lapack_int ldb, double* w );

lapack_int LAPACKE_zhegvd( int matrix_layout, lapack_int itype, char jobz,
                           char uplo, lapack_int n, lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* b,
                           lapack_int ldb, double* w );

lapack_int LAPACKE_zhegvx( int matrix_layout, lapack_int itype, char jobz,
                           char range, char uplo, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_complex_double* b, lapack_int ldb, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w,
                           lapack_complex_double* z, lapack_int ldz,
                           lapack_int* ifail );

lapack_int LAPACKE_zhesv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, lapack_complex_double* a,
                          lapack_int lda, lapack_int* ipiv,
                          lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_zhesvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* af,
                           lapack_int ldaf, lapack_int* ipiv,
                           const lapack_complex_double* b, lapack_int ldb,
                           lapack_complex_double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr );

lapack_int LAPACKE_zhetrd( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* a, lapack_int lda, double* d,
                           double* e, lapack_complex_double* tau );

lapack_int LAPACKE_zhetrf( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_int* ipiv );

lapack_int LAPACKE_zhetri( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           const lapack_int* ipiv );

lapack_int LAPACKE_zhetrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* a,
                           lapack_int lda, const lapack_int* ipiv,
                           lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_zhpcon( int matrix_layout, char uplo, lapack_int n,
                           const lapack_complex_double* ap,
                           const lapack_int* ipiv, double anorm,
                           double* rcond );

lapack_int LAPACKE_zhpev( int matrix_layout, char jobz, char uplo, lapack_int n,
                          lapack_complex_double* ap, double* w,
                          lapack_complex_double* z, lapack_int ldz );

lapack_int LAPACKE_zhpevd( int matrix_layout, char jobz, char uplo, lapack_int n,
                           lapack_complex_double* ap, double* w,
                           lapack_complex_double* z, lapack_int ldz );

lapack_int LAPACKE_zhpevx( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, lapack_complex_double* ap, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w,
                           lapack_complex_double* z, lapack_int ldz,
                           lapack_int* ifail );

lapack_int LAPACKE_zhpgv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, lapack_complex_double* ap,
                          lapack_complex_double* bp, double* w,
                          lapack_complex_double* z, lapack_int ldz );

lapack_int LAPACKE_zhpgvd( int matrix_layout, lapack_int itype, char jobz,
                           char uplo, lapack_int n, lapack_complex_double* ap,
                           lapack_complex_double* bp, double* w,
                           lapack_complex_double* z, lapack_int ldz );

lapack_int LAPACKE_zhpgvx( int matrix_layout, lapack_int itype, char jobz,
                           char range, char uplo, lapack_int n,
                           lapack_complex_double* ap, lapack_complex_double* bp,
                           double vl, double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w,
                           lapack_complex_double* z, lapack_int ldz,
                           lapack_int* ifail );

lapack_int LAPACKE_zhpsv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, lapack_complex_double* ap,
                          lapack_int* ipiv, lapack_complex_double* b,
                          lapack_int ldb );

lapack_int LAPACKE_zhpsvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* ap,
                           lapack_complex_double* afp, lapack_int* ipiv,
                           const lapack_complex_double* b, lapack_int ldb,
                           lapack_complex_double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr );

lapack_int LAPACKE_zhptrd( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* ap, double* d, double* e,
                           lapack_complex_double* tau );

lapack_int LAPACKE_zhptrf( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* ap, lapack_int* ipiv );

lapack_int LAPACKE_zhptri( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* ap, const lapack_int* ipiv );

lapack_int LAPACKE_zhptrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* ap,
                           const lapack_int* ipiv, lapack_complex_double* b,
                           lapack_int ldb );

lapack_int LAPACKE_dhsein( int matrix_layout, char job, char eigsrc, char initv,
                           lapack_logical* select, lapack_int n,
                           const double* h, lapack_int ldh, double* wr,
                           const double* wi, double* vl, lapack_int ldvl,
                           double* vr, lapack_int ldvr, lapack_int mm,
                           lapack_int* m, lapack_int* ifaill,
                           lapack_int* ifailr );
lapack_int LAPACKE_zhsein( int matrix_layout, char job, char eigsrc, char initv,
                           const lapack_logical* select, lapack_int n,
                           const lapack_complex_double* h, lapack_int ldh,
                           lapack_complex_double* w, lapack_complex_double* vl,
                           lapack_int ldvl, lapack_complex_double* vr,
                           lapack_int ldvr, lapack_int mm, lapack_int* m,
                           lapack_int* ifaill, lapack_int* ifailr );

lapack_int LAPACKE_dhseqr( int matrix_layout, char job, char compz, lapack_int n,
                           lapack_int ilo, lapack_int ihi, double* h,
                           lapack_int ldh, double* wr, double* wi, double* z,
                           lapack_int ldz );
lapack_int LAPACKE_zhseqr( int matrix_layout, char job, char compz, lapack_int n,
                           lapack_int ilo, lapack_int ihi,
                           lapack_complex_double* h, lapack_int ldh,
                           lapack_complex_double* w, lapack_complex_double* z,
                           lapack_int ldz );

float LAPACKE_slamch( char cmach );
double LAPACKE_dlamch( char cmach );

double LAPACKE_dlangb( int matrix_layout, char norm, lapack_int n,
                           lapack_int kl, lapack_int ku, const double* ab,
                           lapack_int ldab );
double LAPACKE_zlangb( int matrix_layout, char norm, lapack_int n,
                           lapack_int kl, lapack_int ku,
                           const lapack_complex_double* ab, lapack_int ldab );

double LAPACKE_dlange( int matrix_layout, char norm, lapack_int m,
                           lapack_int n, const double* a, lapack_int lda );
double LAPACKE_zlange( int matrix_layout, char norm, lapack_int m,
                           lapack_int n, const lapack_complex_double* a,
                           lapack_int lda );

double LAPACKE_dlangt( char norm, lapack_int n, const double* dl,
                           const double* d, const double* du );
double LAPACKE_zlangt( char norm, lapack_int n, const lapack_complex_double* dl,
                           const lapack_complex_double* d,
                           const lapack_complex_double* du );

double LAPACKE_zlanhb( int matrix_layout, char norm, char uplo, lapack_int n,
                           lapack_int k, const lapack_complex_double* ab,
                           lapack_int ldab );

double LAPACKE_zlanhe( int matrix_layout, char norm, char uplo, lapack_int n,
                           const lapack_complex_double* a, lapack_int lda );

double LAPACKE_zlanhp( int matrix_layout, char norm, char uplo, lapack_int n,
                           const lapack_complex_double* ap );

double LAPACKE_zlanht( char norm, lapack_int n, const double* d,
                           const lapack_complex_double* e );

double LAPACKE_dlansy( int matrix_layout, char norm, char uplo, lapack_int n,
                           const double* a, lapack_int lda );
double LAPACKE_zlansy( int matrix_layout, char norm, char uplo, lapack_int n,
                           const lapack_complex_double* a, lapack_int lda );

double LAPACKE_dlansb( int matrix_layout, char norm, char uplo, lapack_int n,
                           lapack_int k, const double* ab, lapack_int ldab );

double LAPACKE_dlansp( int matrix_layout, char norm, char uplo, lapack_int n,
                           const double* ap );
double LAPACKE_zlansp( int matrix_layout, char norm, char uplo, lapack_int n,
                           const lapack_complex_double* ap );

double LAPACKE_dlanst( char norm, lapack_int n, const double* d, const double* e );

double LAPACKE_dlantr( int matrix_layout, char norm, char uplo, char diag,
                           lapack_int m, lapack_int n, const double* a,
                           lapack_int lda );
double LAPACKE_zlantr( int matrix_layout, char norm, char uplo, char diag,
                           lapack_int m, lapack_int n, const lapack_complex_double* a,
                           lapack_int lda );

lapack_int LAPACKE_dlatms( int matrix_layout, lapack_int m, lapack_int n,
                           char dist, lapack_int* iseed, char sym, double* d,
                           lapack_int mode, double cond, double dmax,
                           lapack_int kl, lapack_int ku, char pack, double* a,
                           lapack_int lda );
lapack_int LAPACKE_zlatms( int matrix_layout, lapack_int m, lapack_int n,
                           char dist, lapack_int* iseed, char sym, double* d,
                           lapack_int mode, double cond, double dmax,
                           lapack_int kl, lapack_int ku, char pack,
                           lapack_complex_double* a, lapack_int lda );

lapack_int LAPACKE_dopgtr( int matrix_layout, char uplo, lapack_int n,
                           const double* ap, const double* tau, double* q,
                           lapack_int ldq );

lapack_int LAPACKE_dopmtr( int matrix_layout, char side, char uplo, char trans,
                           lapack_int m, lapack_int n, const double* ap,
                           const double* tau, double* c, lapack_int ldc );

lapack_int LAPACKE_dorghr( int matrix_layout, lapack_int n, lapack_int ilo,
                           lapack_int ihi, double* a, lapack_int lda,
                           const double* tau );

lapack_int LAPACKE_dorglq( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int k, double* a, lapack_int lda,
                           const double* tau );

lapack_int LAPACKE_dorgqr( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int k, double* a, lapack_int lda,
                           const double* tau );

lapack_int LAPACKE_dorgtr( int matrix_layout, char uplo, lapack_int n, double* a,
                           lapack_int lda, const double* tau );

lapack_int LAPACKE_dormhr( int matrix_layout, char side, char trans,
                           lapack_int m, lapack_int n, lapack_int ilo,
                           lapack_int ihi, const double* a, lapack_int lda,
                           const double* tau, double* c, lapack_int ldc );

lapack_int LAPACKE_dormlq( int matrix_layout, char side, char trans,
                           lapack_int m, lapack_int n, lapack_int k,
                           const double* a, lapack_int lda, const double* tau,
                           double* c, lapack_int ldc );

lapack_int LAPACKE_dormqr( int matrix_layout, char side, char trans,
                           lapack_int m, lapack_int n, lapack_int k,
                           const double* a, lapack_int lda, const double* tau,
                           double* c, lapack_int ldc );

lapack_int LAPACKE_dormtr( int matrix_layout, char side, char uplo, char trans,
                           lapack_int m, lapack_int n, const double* a,
                           lapack_int lda, const double* tau, double* c,
                           lapack_int ldc );

lapack_int LAPACKE_dpbcon( int matrix_layout, char uplo, lapack_int n,
                           lapack_int kd, const double* ab, lapack_int ldab,
                           double anorm, double* rcond );
lapack_int LAPACKE_zpbcon( int matrix_layout, char uplo, lapack_int n,
                           lapack_int kd, const lapack_complex_double* ab,
                           lapack_int ldab, double anorm, double* rcond );

lapack_int LAPACKE_dpbsv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int kd, lapack_int nrhs, double* ab,
                          lapack_int ldab, double* b, lapack_int ldb );
lapack_int LAPACKE_zpbsv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int kd, lapack_int nrhs,
                          lapack_complex_double* ab, lapack_int ldab,
                          lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dpbsvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int kd, lapack_int nrhs, double* ab,
                           lapack_int ldab, double* afb, lapack_int ldafb,
                           char* equed, double* s, double* b, lapack_int ldb,
                           double* x, lapack_int ldx, double* rcond,
                           double* ferr, double* berr );
lapack_int LAPACKE_zpbsvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int kd, lapack_int nrhs,
                           lapack_complex_double* ab, lapack_int ldab,
                           lapack_complex_double* afb, lapack_int ldafb,
                           char* equed, double* s, lapack_complex_double* b,
                           lapack_int ldb, lapack_complex_double* x,
                           lapack_int ldx, double* rcond, double* ferr,
                           double* berr );

lapack_int LAPACKE_dpbtrf( int matrix_layout, char uplo, lapack_int n,
                           lapack_int kd, double* ab, lapack_int ldab );
lapack_int LAPACKE_zpbtrf( int matrix_layout, char uplo, lapack_int n,
                           lapack_int kd, lapack_complex_double* ab,
                           lapack_int ldab );

lapack_int LAPACKE_dpbtrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int kd, lapack_int nrhs, const double* ab,
                           lapack_int ldab, double* b, lapack_int ldb );
lapack_int LAPACKE_zpbtrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int kd, lapack_int nrhs,
                           const lapack_complex_double* ab, lapack_int ldab,
                           lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dpocon( int matrix_layout, char uplo, lapack_int n,
                           const double* a, lapack_int lda, double anorm,
                           double* rcond );
lapack_int LAPACKE_zpocon( int matrix_layout, char uplo, lapack_int n,
                           const lapack_complex_double* a, lapack_int lda,
                           double anorm, double* rcond );

lapack_int LAPACKE_dposv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, double* a, lapack_int lda, double* b,
                          lapack_int ldb );
lapack_int LAPACKE_zposv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, lapack_complex_double* a,
                          lapack_int lda, lapack_complex_double* b,
                          lapack_int ldb );
lapack_int LAPACKE_dsposv( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, double* a, lapack_int lda,
                           double* b, lapack_int ldb, double* x, lapack_int ldx,
                           lapack_int* iter );
lapack_int LAPACKE_zcposv( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* b,
                           lapack_int ldb, lapack_complex_double* x,
                           lapack_int ldx, lapack_int* iter );

lapack_int LAPACKE_dposvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int nrhs, double* a, lapack_int lda,
                           double* af, lapack_int ldaf, char* equed, double* s,
                           double* b, lapack_int ldb, double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr );
lapack_int LAPACKE_zposvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int nrhs, lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* af,
                           lapack_int ldaf, char* equed, double* s,
                           lapack_complex_double* b, lapack_int ldb,
                           lapack_complex_double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr );

lapack_int LAPACKE_dpotrf( int matrix_layout, char uplo, lapack_int n, double* a,
                           lapack_int lda );
lapack_int LAPACKE_zpotrf( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* a, lapack_int lda );

lapack_int LAPACKE_dpotri( int matrix_layout, char uplo, lapack_int n, double* a,
                           lapack_int lda );
lapack_int LAPACKE_zpotri( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* a, lapack_int lda );

lapack_int LAPACKE_dpotrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, const double* a, lapack_int lda,
                           double* b, lapack_int ldb );
lapack_int LAPACKE_zpotrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* b,
                           lapack_int ldb );

lapack_int LAPACKE_dppcon( int matrix_layout, char uplo, lapack_int n,
                           const double* ap, double anorm, double* rcond );
lapack_int LAPACKE_zppcon( int matrix_layout, char uplo, lapack_int n,
                           const lapack_complex_double* ap, double anorm,
                           double* rcond );

lapack_int LAPACKE_dppsv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, double* ap, double* b,
                          lapack_int ldb );
lapack_int LAPACKE_zppsv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, lapack_complex_double* ap,
                          lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dppsvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int nrhs, double* ap, double* afp,
                           char* equed, double* s, double* b, lapack_int ldb,
                           double* x, lapack_int ldx, double* rcond,
                           double* ferr, double* berr );
lapack_int LAPACKE_zppsvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int nrhs, lapack_complex_double* ap,
                           lapack_complex_double* afp, char* equed, double* s,
                           lapack_complex_double* b, lapack_int ldb,
                           lapack_complex_double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr );

lapack_int LAPACKE_dpptrf( int matrix_layout, char uplo, lapack_int n,
                           double* ap );
lapack_int LAPACKE_zpptrf( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* ap );

lapack_int LAPACKE_dpptri( int matrix_layout, char uplo, lapack_int n,
                           double* ap );
lapack_int LAPACKE_zpptri( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* ap );

lapack_int LAPACKE_dpptrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, const double* ap, double* b,
                           lapack_int ldb );
lapack_int LAPACKE_zpptrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* ap,
                           lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dptcon( lapack_int n, const double* d, const double* e,
                           double anorm, double* rcond );
lapack_int LAPACKE_zptcon( lapack_int n, const double* d,
                           const lapack_complex_double* e, double anorm,
                           double* rcond );

lapack_int LAPACKE_dpteqr( int matrix_layout, char compz, lapack_int n,
                           double* d, double* e, double* z, lapack_int ldz );
lapack_int LAPACKE_zpteqr( int matrix_layout, char compz, lapack_int n,
                           double* d, double* e, lapack_complex_double* z,
                           lapack_int ldz );

lapack_int LAPACKE_dptsv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* d, double* e, double* b, lapack_int ldb );
lapack_int LAPACKE_zptsv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* d, lapack_complex_double* e,
                          lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dptsvx( int matrix_layout, char fact, lapack_int n,
                           lapack_int nrhs, const double* d, const double* e,
                           double* df, double* ef, const double* b,
                           lapack_int ldb, double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr );
lapack_int LAPACKE_zptsvx( int matrix_layout, char fact, lapack_int n,
                           lapack_int nrhs, const double* d,
                           const lapack_complex_double* e, double* df,
                           lapack_complex_double* ef,
                           const lapack_complex_double* b, lapack_int ldb,
                           lapack_complex_double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr );

lapack_int LAPACKE_dpttrf( lapack_int n, double* d, double* e );
lapack_int LAPACKE_zpttrf( lapack_int n, double* d, lapack_complex_double* e );

lapack_int LAPACKE_dpttrs( int matrix_layout, lapack_int n, lapack_int nrhs,
                           const double* d, const double* e, double* b,
                           lapack_int ldb );
lapack_int LAPACKE_zpttrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, const double* d,
                           const lapack_complex_double* e,
                           lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dsbev( int matrix_layout, char jobz, char uplo, lapack_int n,
                          lapack_int kd, double* ab, lapack_int ldab, double* w,
                          double* z, lapack_int ldz );

lapack_int LAPACKE_dsbevd( int matrix_layout, char jobz, char uplo, lapack_int n,
                           lapack_int kd, double* ab, lapack_int ldab,
                           double* w, double* z, lapack_int ldz );

lapack_int LAPACKE_dsbevx( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, lapack_int kd, double* ab,
                           lapack_int ldab, double* q, lapack_int ldq,
                           double vl, double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w, double* z,
                           lapack_int ldz, lapack_int* ifail );

lapack_int LAPACKE_dsbgv( int matrix_layout, char jobz, char uplo, lapack_int n,
                          lapack_int ka, lapack_int kb, double* ab,
                          lapack_int ldab, double* bb, lapack_int ldbb,
                          double* w, double* z, lapack_int ldz );

lapack_int LAPACKE_dsbgvd( int matrix_layout, char jobz, char uplo, lapack_int n,
                           lapack_int ka, lapack_int kb, double* ab,
                           lapack_int ldab, double* bb, lapack_int ldbb,
                           double* w, double* z, lapack_int ldz );

lapack_int LAPACKE_dsbgvx( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, lapack_int ka, lapack_int kb,
                           double* ab, lapack_int ldab, double* bb,
                           lapack_int ldbb, double* q, lapack_int ldq,
                           double vl, double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w, double* z,
                           lapack_int ldz, lapack_int* ifail );

lapack_int LAPACKE_dsbtrd( int matrix_layout, char vect, char uplo, lapack_int n,
                           lapack_int kd, double* ab, lapack_int ldab,
                           double* d, double* e, double* q, lapack_int ldq );

lapack_int LAPACKE_dspev( int matrix_layout, char jobz, char uplo, lapack_int n,
                          double* ap, double* w, double* z, lapack_int ldz );

lapack_int LAPACKE_dspevd( int matrix_layout, char jobz, char uplo, lapack_int n,
                           double* ap, double* w, double* z, lapack_int ldz );

lapack_int LAPACKE_dspevx( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, double* ap, double vl, double vu,
                           lapack_int il, lapack_int iu, double abstol,
                           lapack_int* m, double* w, double* z, lapack_int ldz,
                           lapack_int* ifail );

lapack_int LAPACKE_dspgv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, double* ap, double* bp,
                          double* w, double* z, lapack_int ldz );

lapack_int LAPACKE_dspgvd( int matrix_layout, lapack_int itype, char jobz,
                           char uplo, lapack_int n, double* ap, double* bp,
                           double* w, double* z, lapack_int ldz );

lapack_int LAPACKE_dspgvx( int matrix_layout, lapack_int itype, char jobz,
                           char range, char uplo, lapack_int n, double* ap,
                           double* bp, double vl, double vu, lapack_int il,
                           lapack_int iu, double abstol, lapack_int* m,
                           double* w, double* z, lapack_int ldz,
                           lapack_int* ifail );

lapack_int LAPACKE_dspcon( int matrix_layout, char uplo, lapack_int n,
                           const double* ap, const lapack_int* ipiv,
                           double anorm, double* rcond );
lapack_int LAPACKE_zspcon( int matrix_layout, char uplo, lapack_int n,
                           const lapack_complex_double* ap,
                           const lapack_int* ipiv, double anorm,
                           double* rcond );

lapack_int LAPACKE_dspsv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, double* ap, lapack_int* ipiv,
                          double* b, lapack_int ldb );
lapack_int LAPACKE_zspsv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, lapack_complex_double* ap,
                          lapack_int* ipiv, lapack_complex_double* b,
                          lapack_int ldb );

lapack_int LAPACKE_dspsvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int nrhs, const double* ap, double* afp,
                           lapack_int* ipiv, const double* b, lapack_int ldb,
                           double* x, lapack_int ldx, double* rcond,
                           double* ferr, double* berr );
lapack_int LAPACKE_zspsvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* ap,
                           lapack_complex_double* afp, lapack_int* ipiv,
                           const lapack_complex_double* b, lapack_int ldb,
                           lapack_complex_double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr );

lapack_int LAPACKE_dsptrd( int matrix_layout, char uplo, lapack_int n,
                           double* ap, double* d, double* e, double* tau );

lapack_int LAPACKE_dsptrf( int matrix_layout, char uplo, lapack_int n,
                           double* ap, lapack_int* ipiv );
lapack_int LAPACKE_zsptrf( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* ap, lapack_int* ipiv );

lapack_int LAPACKE_dsptri( int matrix_layout, char uplo, lapack_int n,
                           double* ap, const lapack_int* ipiv );
lapack_int LAPACKE_zsptri( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* ap, const lapack_int* ipiv );

lapack_int LAPACKE_dsptrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, const double* ap,
                           const lapack_int* ipiv, double* b, lapack_int ldb );
lapack_int LAPACKE_zsptrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* ap,
                           const lapack_int* ipiv, lapack_complex_double* b,
                           lapack_int ldb );

lapack_int LAPACKE_dstebz( char range, char order, lapack_int n, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, const double* d, const double* e,
                           lapack_int* m, lapack_int* nsplit, double* w,
                           lapack_int* iblock, lapack_int* isplit );

lapack_int LAPACKE_dstedc( int matrix_layout, char compz, lapack_int n,
                           double* d, double* e, double* z, lapack_int ldz );
lapack_int LAPACKE_zstedc( int matrix_layout, char compz, lapack_int n,
                           double* d, double* e, lapack_complex_double* z,
                           lapack_int ldz );

lapack_int LAPACKE_dstein( int matrix_layout, lapack_int n, const double* d,
                           const double* e, lapack_int m, const double* w,
                           const lapack_int* iblock, const lapack_int* isplit,
                           double* z, lapack_int ldz, lapack_int* ifailv );
lapack_int LAPACKE_zstein( int matrix_layout, lapack_int n, const double* d,
                           const double* e, lapack_int m, const double* w,
                           const lapack_int* iblock, const lapack_int* isplit,
                           lapack_complex_double* z, lapack_int ldz,
                           lapack_int* ifailv );

lapack_int LAPACKE_dstemr( int matrix_layout, char jobz, char range,
                           lapack_int n, double* d, double* e, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           lapack_int* m, double* w, double* z, lapack_int ldz,
                           lapack_int nzc, lapack_int* isuppz,
                           lapack_logical* tryrac );
lapack_int LAPACKE_zstemr( int matrix_layout, char jobz, char range,
                           lapack_int n, double* d, double* e, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           lapack_int* m, double* w, lapack_complex_double* z,
                           lapack_int ldz, lapack_int nzc, lapack_int* isuppz,
                           lapack_logical* tryrac );

lapack_int LAPACKE_dsteqr( int matrix_layout, char compz, lapack_int n,
                           double* d, double* e, double* z, lapack_int ldz );
lapack_int LAPACKE_zsteqr( int matrix_layout, char compz, lapack_int n,
                           double* d, double* e, lapack_complex_double* z,
                           lapack_int ldz );

lapack_int LAPACKE_dsterf( lapack_int n, double* d, double* e );

lapack_int LAPACKE_dstev( int matrix_layout, char jobz, lapack_int n, double* d,
                          double* e, double* z, lapack_int ldz );

lapack_int LAPACKE_dstevd( int matrix_layout, char jobz, lapack_int n, double* d,
                           double* e, double* z, lapack_int ldz );

lapack_int LAPACKE_dstevr( int matrix_layout, char jobz, char range,
                           lapack_int n, double* d, double* e, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w, double* z,
                           lapack_int ldz, lapack_int* isuppz );

lapack_int LAPACKE_dstevx( int matrix_layout, char jobz, char range,
                           lapack_int n, double* d, double* e, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w, double* z,
                           lapack_int ldz, lapack_int* ifail );

lapack_int LAPACKE_dsycon( int matrix_layout, char uplo, lapack_int n,
                           const double* a, lapack_int lda,
                           const lapack_int* ipiv, double anorm,
                           double* rcond );
lapack_int LAPACKE_zsycon( int matrix_layout, char uplo, lapack_int n,
                           const lapack_complex_double* a, lapack_int lda,
                           const lapack_int* ipiv, double anorm,
                           double* rcond );

lapack_int LAPACKE_dsyev( int matrix_layout, char jobz, char uplo, lapack_int n,
                          double* a, lapack_int lda, double* w );

lapack_int LAPACKE_dsyevd( int matrix_layout, char jobz, char uplo, lapack_int n,
                           double* a, lapack_int lda, double* w );

lapack_int LAPACKE_dsyevr( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, double* a, lapack_int lda, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w, double* z,
                           lapack_int ldz, lapack_int* isuppz );

lapack_int LAPACKE_dsyevx( int matrix_layout, char jobz, char range, char uplo,
                           lapack_int n, double* a, lapack_int lda, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w, double* z,
                           lapack_int ldz, lapack_int* ifail );

lapack_int LAPACKE_dsygv( int matrix_layout, lapack_int itype, char jobz,
                          char uplo, lapack_int n, double* a, lapack_int lda,
                          double* b, lapack_int ldb, double* w );

lapack_int LAPACKE_dsygvd( int matrix_layout, lapack_int itype, char jobz,
                           char uplo, lapack_int n, double* a, lapack_int lda,
                           double* b, lapack_int ldb, double* w );

lapack_int LAPACKE_dsygvx( int matrix_layout, lapack_int itype, char jobz,
                           char range, char uplo, lapack_int n, double* a,
                           lapack_int lda, double* b, lapack_int ldb, double vl,
                           double vu, lapack_int il, lapack_int iu,
                           double abstol, lapack_int* m, double* w, double* z,
                           lapack_int ldz, lapack_int* ifail );

lapack_int LAPACKE_zsyr( int matrix_layout, char uplo, lapack_int n,
                             lapack_complex_double alpha,
                             const lapack_complex_double* x, lapack_int incx,
                             lapack_complex_double* a, lapack_int lda );

lapack_int LAPACKE_dsysv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, double* a, lapack_int lda,
                          lapack_int* ipiv, double* b, lapack_int ldb );
lapack_int LAPACKE_zsysv( int matrix_layout, char uplo, lapack_int n,
                          lapack_int nrhs, lapack_complex_double* a,
                          lapack_int lda, lapack_int* ipiv,
                          lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dsysvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int nrhs, const double* a, lapack_int lda,
                           double* af, lapack_int ldaf, lapack_int* ipiv,
                           const double* b, lapack_int ldb, double* x,
                           lapack_int ldx, double* rcond, double* ferr,
                           double* berr );
lapack_int LAPACKE_zsysvx( int matrix_layout, char fact, char uplo, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* a,
                           lapack_int lda, lapack_complex_double* af,
                           lapack_int ldaf, lapack_int* ipiv,
                           const lapack_complex_double* b, lapack_int ldb,
                           lapack_complex_double* x, lapack_int ldx,
                           double* rcond, double* ferr, double* berr );

lapack_int LAPACKE_dsytrd( int matrix_layout, char uplo, lapack_int n, double* a,
                           lapack_int lda, double* d, double* e, double* tau );

lapack_int LAPACKE_dsytrf( int matrix_layout, char uplo, lapack_int n, double* a,
                           lapack_int lda, lapack_int* ipiv );
lapack_int LAPACKE_zsytrf( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           lapack_int* ipiv );

lapack_int LAPACKE_dsytri( int matrix_layout, char uplo, lapack_int n, double* a,
                           lapack_int lda, const lapack_int* ipiv );
lapack_int LAPACKE_zsytri( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           const lapack_int* ipiv );

lapack_int LAPACKE_dsytrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, const double* a, lapack_int lda,
                           const lapack_int* ipiv, double* b, lapack_int ldb );
lapack_int LAPACKE_zsytrs( int matrix_layout, char uplo, lapack_int n,
                           lapack_int nrhs, const lapack_complex_double* a,
                           lapack_int lda, const lapack_int* ipiv,
                           lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dtbcon( int matrix_layout, char norm, char uplo, char diag,
                           lapack_int n, lapack_int kd, const double* ab,
                           lapack_int ldab, double* rcond );
lapack_int LAPACKE_ztbcon( int matrix_layout, char norm, char uplo, char diag,
                           lapack_int n, lapack_int kd,
                           const lapack_complex_double* ab, lapack_int ldab,
                           double* rcond );

lapack_int LAPACKE_dtbtrs( int matrix_layout, char uplo, char trans, char diag,
                           lapack_int n, lapack_int kd, lapack_int nrhs,
                           const double* ab, lapack_int ldab, double* b,
                           lapack_int ldb );
lapack_int LAPACKE_ztbtrs( int matrix_layout, char uplo, char trans, char diag,
                           lapack_int n, lapack_int kd, lapack_int nrhs,
                           const lapack_complex_double* ab, lapack_int ldab,
                           lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dtpcon( int matrix_layout, char norm, char uplo, char diag,
                           lapack_int n, const double* ap, double* rcond );
lapack_int LAPACKE_ztpcon( int matrix_layout, char norm, char uplo, char diag,
                           lapack_int n, const lapack_complex_double* ap,
                           double* rcond );

lapack_int LAPACKE_dtptri( int matrix_layout, char uplo, char diag, lapack_int n,
                           double* ap );
lapack_int LAPACKE_ztptri( int matrix_layout, char uplo, char diag, lapack_int n,
                           lapack_complex_double* ap );

lapack_int LAPACKE_dtptrs( int matrix_layout, char uplo, char trans, char diag,
                           lapack_int n, lapack_int nrhs, const double* ap,
                           double* b, lapack_int ldb );
lapack_int LAPACKE_ztptrs( int matrix_layout, char uplo, char trans, char diag,
                           lapack_int n, lapack_int nrhs,
                           const lapack_complex_double* ap,
                           lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dtrcon( int matrix_layout, char norm, char uplo, char diag,
                           lapack_int n, const double* a, lapack_int lda,
                           double* rcond );
lapack_int LAPACKE_ztrcon( int matrix_layout, char norm, char uplo, char diag,
                           lapack_int n, const lapack_complex_double* a,
                           lapack_int lda, double* rcond );

lapack_int LAPACKE_dtrevc3(int matrix_layout, char side, char howmny,
                           lapack_logical* select, lapack_int n, const double* t,
                           lapack_int ldt, double* vl, lapack_int ldvl, double* vr,
                           lapack_int ldvr, lapack_int mm, lapack_int* m);
lapack_int LAPACKE_ztrevc3(int matrix_layout, char side, char howmny,
                           const lapack_logical* select, lapack_int n, lapack_complex_double* t,
                           lapack_int ldt, lapack_complex_double* vl, lapack_int ldvl,
                           lapack_complex_double* vr, lapack_int ldvr, lapack_int mm, lapack_int* m);

lapack_int LAPACKE_dtrexc( int matrix_layout, char compq, lapack_int n,
                           double* t, lapack_int ldt, double* q, lapack_int ldq,
                           lapack_int* ifst, lapack_int* ilst );
lapack_int LAPACKE_ztrexc( int matrix_layout, char compq, lapack_int n,
                           lapack_complex_double* t, lapack_int ldt,
                           lapack_complex_double* q, lapack_int ldq,
                           lapack_int ifst, lapack_int ilst );

lapack_int LAPACKE_dtrsen( int matrix_layout, char job, char compq,
                           const lapack_logical* select, lapack_int n,
                           double* t, lapack_int ldt, double* q, lapack_int ldq,
                           double* wr, double* wi, lapack_int* m, double* s,
                           double* sep );
lapack_int LAPACKE_ztrsen( int matrix_layout, char job, char compq,
                           const lapack_logical* select, lapack_int n,
                           lapack_complex_double* t, lapack_int ldt,
                           lapack_complex_double* q, lapack_int ldq,
                           lapack_complex_double* w, lapack_int* m, double* s,
                           double* sep );

lapack_int LAPACKE_dtrsna( int matrix_layout, char job, char howmny,
                           const lapack_logical* select, lapack_int n,
                           const double* t, lapack_int ldt, const double* vl,
                           lapack_int ldvl, const double* vr, lapack_int ldvr,
                           double* s, double* sep, lapack_int mm,
                           lapack_int* m );
lapack_int LAPACKE_ztrsna( int matrix_layout, char job, char howmny,
                           const lapack_logical* select, lapack_int n,
                           const lapack_complex_double* t, lapack_int ldt,
                           const lapack_complex_double* vl, lapack_int ldvl,
                           const lapack_complex_double* vr, lapack_int ldvr,
                           double* s, double* sep, lapack_int mm,
                           lapack_int* m );

lapack_int LAPACKE_dtrsyl( int matrix_layout, char trana, char tranb,
                           lapack_int isgn, lapack_int m, lapack_int n,
                           const double* a, lapack_int lda, const double* b,
                           lapack_int ldb, double* c, lapack_int ldc,
                           double* scale );
lapack_int LAPACKE_ztrsyl( int matrix_layout, char trana, char tranb,
                           lapack_int isgn, lapack_int m, lapack_int n,
                           const lapack_complex_double* a, lapack_int lda,
                           const lapack_complex_double* b, lapack_int ldb,
                           lapack_complex_double* c, lapack_int ldc,
                           double* scale );

lapack_int LAPACKE_dtrtri( int matrix_layout, char uplo, char diag, lapack_int n,
                           double* a, lapack_int lda );
lapack_int LAPACKE_ztrtri( int matrix_layout, char uplo, char diag, lapack_int n,
                           lapack_complex_double* a, lapack_int lda );

lapack_int LAPACKE_dtrtrs( int matrix_layout, char uplo, char trans, char diag,
                           lapack_int n, lapack_int nrhs, const double* a,
                           lapack_int lda, double* b, lapack_int ldb );
lapack_int LAPACKE_ztrtrs( int matrix_layout, char uplo, char trans, char diag,
                           lapack_int n, lapack_int nrhs,
                           const lapack_complex_double* a, lapack_int lda,
                           lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_zunghr( int matrix_layout, lapack_int n, lapack_int ilo,
                           lapack_int ihi, lapack_complex_double* a,
                           lapack_int lda, const lapack_complex_double* tau );

lapack_int LAPACKE_zunglq( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int k, lapack_complex_double* a,
                           lapack_int lda, const lapack_complex_double* tau );

lapack_int LAPACKE_zungqr( int matrix_layout, lapack_int m, lapack_int n,
                           lapack_int k, lapack_complex_double* a,
                           lapack_int lda, const lapack_complex_double* tau );

lapack_int LAPACKE_zungtr( int matrix_layout, char uplo, lapack_int n,
                           lapack_complex_double* a, lapack_int lda,
                           const lapack_complex_double* tau );

lapack_int LAPACKE_zunmhr( int matrix_layout, char side, char trans,
                           lapack_int m, lapack_int n, lapack_int ilo,
                           lapack_int ihi, const lapack_complex_double* a,
                           lapack_int lda, const lapack_complex_double* tau,
                           lapack_complex_double* c, lapack_int ldc );

lapack_int LAPACKE_zunmlq( int matrix_layout, char side, char trans,
                           lapack_int m, lapack_int n, lapack_int k,
                           const lapack_complex_double* a, lapack_int lda,
                           const lapack_complex_double* tau,
                           lapack_complex_double* c, lapack_int ldc );

lapack_int LAPACKE_zunmqr( int matrix_layout, char side, char trans,
                           lapack_int m, lapack_int n, lapack_int k,
                           const lapack_complex_double* a, lapack_int lda,
                           const lapack_complex_double* tau,
                           lapack_complex_double* c, lapack_int ldc );

lapack_int LAPACKE_zunmtr( int matrix_layout, char side, char uplo, char trans,
                           lapack_int m, lapack_int n,
                           const lapack_complex_double* a, lapack_int lda,
                           const lapack_complex_double* tau,
                           lapack_complex_double* c, lapack_int ldc );

lapack_int LAPACKE_zupgtr( int matrix_layout, char uplo, lapack_int n,
                           const lapack_complex_double* ap,
                           const lapack_complex_double* tau,
                           lapack_complex_double* q, lapack_int ldq );

lapack_int LAPACKE_zupmtr( int matrix_layout, char side, char uplo, char trans,
                           lapack_int m, lapack_int n,
                           const lapack_complex_double* ap,
                           const lapack_complex_double* tau,
                           lapack_complex_double* c, lapack_int ldc );

void LAPACKE_ilaver( lapack_int* vers_major,
                     lapack_int* vers_minor,
                     lapack_int* vers_patch );

/*----------*/

lapack_int LAPACKE_ddisna_work( char job, lapack_int m, lapack_int n,
                                const double* d, double* sep );

lapack_int LAPACKE_dgbcon_work( int matrix_layout, char norm, lapack_int n,
                                lapack_int kl, lapack_int ku, const double* ab,
                                lapack_int ldab, const lapack_int* ipiv,
                                double anorm, double* rcond, double* work,
                                lapack_int* iwork );
lapack_int LAPACKE_zgbcon_work( int matrix_layout, char norm, lapack_int n,
                                lapack_int kl, lapack_int ku,
                                const lapack_complex_double* ab,
                                lapack_int ldab, const lapack_int* ipiv,
                                double anorm, double* rcond,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_dgbsv_work( int matrix_layout, lapack_int n, lapack_int kl,
                               lapack_int ku, lapack_int nrhs, double* ab,
                               lapack_int ldab, lapack_int* ipiv, double* b,
                               lapack_int ldb );
lapack_int LAPACKE_zgbsv_work( int matrix_layout, lapack_int n, lapack_int kl,
                               lapack_int ku, lapack_int nrhs,
                               lapack_complex_double* ab, lapack_int ldab,
                               lapack_int* ipiv, lapack_complex_double* b,
                               lapack_int ldb );

lapack_int LAPACKE_dgbsvx_work( int matrix_layout, char fact, char trans,
                                lapack_int n, lapack_int kl, lapack_int ku,
                                lapack_int nrhs, double* ab, lapack_int ldab,
                                double* afb, lapack_int ldafb, lapack_int* ipiv,
                                char* equed, double* r, double* c, double* b,
                                lapack_int ldb, double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, lapack_int* iwork );
lapack_int LAPACKE_zgbsvx_work( int matrix_layout, char fact, char trans,
                                lapack_int n, lapack_int kl, lapack_int ku,
                                lapack_int nrhs, lapack_complex_double* ab,
                                lapack_int ldab, lapack_complex_double* afb,
                                lapack_int ldafb, lapack_int* ipiv, char* equed,
                                double* r, double* c, lapack_complex_double* b,
                                lapack_int ldb, lapack_complex_double* x,
                                lapack_int ldx, double* rcond, double* ferr,
                                double* berr, lapack_complex_double* work,
                                double* rwork );

lapack_int LAPACKE_dgbtrf_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int kl, lapack_int ku, double* ab,
                                lapack_int ldab, lapack_int* ipiv );
lapack_int LAPACKE_zgbtrf_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int kl, lapack_int ku,
                                lapack_complex_double* ab, lapack_int ldab,
                                lapack_int* ipiv );

lapack_int LAPACKE_dgbtrs_work( int matrix_layout, char trans, lapack_int n,
                                lapack_int kl, lapack_int ku, lapack_int nrhs,
                                const double* ab, lapack_int ldab,
                                const lapack_int* ipiv, double* b,
                                lapack_int ldb );
lapack_int LAPACKE_zgbtrs_work( int matrix_layout, char trans, lapack_int n,
                                lapack_int kl, lapack_int ku, lapack_int nrhs,
                                const lapack_complex_double* ab,
                                lapack_int ldab, const lapack_int* ipiv,
                                lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dgebak_work( int matrix_layout, char job, char side,
                                lapack_int n, lapack_int ilo, lapack_int ihi,
                                const double* scale, lapack_int m, double* v,
                                lapack_int ldv );
lapack_int LAPACKE_zgebak_work( int matrix_layout, char job, char side,
                                lapack_int n, lapack_int ilo, lapack_int ihi,
                                const double* scale, lapack_int m,
                                lapack_complex_double* v, lapack_int ldv );

lapack_int LAPACKE_dgebal_work( int matrix_layout, char job, lapack_int n,
                                double* a, lapack_int lda, lapack_int* ilo,
                                lapack_int* ihi, double* scale );
lapack_int LAPACKE_zgebal_work( int matrix_layout, char job, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_int* ilo, lapack_int* ihi,
                                double* scale );

lapack_int LAPACKE_dgecon_work( int matrix_layout, char norm, lapack_int n,
                                const double* a, lapack_int lda, double anorm,
                                double* rcond, double* work,
                                lapack_int* iwork );
lapack_int LAPACKE_zgecon_work( int matrix_layout, char norm, lapack_int n,
                                const lapack_complex_double* a, lapack_int lda,
                                double anorm, double* rcond,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_dgees_work( int matrix_layout, char jobvs, char sort,
                               LAPACK_D_SELECT2 select, lapack_int n, double* a,
                               lapack_int lda, lapack_int* sdim, double* wr,
                               double* wi, double* vs, lapack_int ldvs,
                               double* work, lapack_int lwork,
                               lapack_logical* bwork );
lapack_int LAPACKE_zgees_work( int matrix_layout, char jobvs, char sort,
                               LAPACK_Z_SELECT1 select, lapack_int n,
                               lapack_complex_double* a, lapack_int lda,
                               lapack_int* sdim, lapack_complex_double* w,
                               lapack_complex_double* vs, lapack_int ldvs,
                               lapack_complex_double* work, lapack_int lwork,
                               double* rwork, lapack_logical* bwork );

lapack_int LAPACKE_dgeesx_work( int matrix_layout, char jobvs, char sort,
                                LAPACK_D_SELECT2 select, char sense,
                                lapack_int n, double* a, lapack_int lda,
                                lapack_int* sdim, double* wr, double* wi,
                                double* vs, lapack_int ldvs, double* rconde,
                                double* rcondv, double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork,
                                lapack_logical* bwork );
lapack_int LAPACKE_zgeesx_work( int matrix_layout, char jobvs, char sort,
                                LAPACK_Z_SELECT1 select, char sense,
                                lapack_int n, lapack_complex_double* a,
                                lapack_int lda, lapack_int* sdim,
                                lapack_complex_double* w,
                                lapack_complex_double* vs, lapack_int ldvs,
                                double* rconde, double* rcondv,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork, lapack_logical* bwork );

lapack_int LAPACKE_dgeev_work( int matrix_layout, char jobvl, char jobvr,
                               lapack_int n, double* a, lapack_int lda,
                               double* wr, double* wi, double* vl,
                               lapack_int ldvl, double* vr, lapack_int ldvr,
                               double* work, lapack_int lwork );
lapack_int LAPACKE_zgeev_work( int matrix_layout, char jobvl, char jobvr,
                               lapack_int n, lapack_complex_double* a,
                               lapack_int lda, lapack_complex_double* w,
                               lapack_complex_double* vl, lapack_int ldvl,
                               lapack_complex_double* vr, lapack_int ldvr,
                               lapack_complex_double* work, lapack_int lwork,
                               double* rwork );

lapack_int LAPACKE_dgeevx_work( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, lapack_int n, double* a,
                                lapack_int lda, double* wr, double* wi,
                                double* vl, lapack_int ldvl, double* vr,
                                lapack_int ldvr, lapack_int* ilo,
                                lapack_int* ihi, double* scale, double* abnrm,
                                double* rconde, double* rcondv, double* work,
                                lapack_int lwork, lapack_int* iwork );
lapack_int LAPACKE_zgeevx_work( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_complex_double* w,
                                lapack_complex_double* vl, lapack_int ldvl,
                                lapack_complex_double* vr, lapack_int ldvr,
                                lapack_int* ilo, lapack_int* ihi, double* scale,
                                double* abnrm, double* rconde, double* rcondv,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork );

lapack_int LAPACKE_dgehrd_work( int matrix_layout, lapack_int n, lapack_int ilo,
                                lapack_int ihi, double* a, lapack_int lda,
                                double* tau, double* work, lapack_int lwork );
lapack_int LAPACKE_zgehrd_work( int matrix_layout, lapack_int n, lapack_int ilo,
                                lapack_int ihi, lapack_complex_double* a,
                                lapack_int lda, lapack_complex_double* tau,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_dgejsv_work( int matrix_layout, char joba, char jobu,
                                char jobv, char jobr, char jobt, char jobp,
                                lapack_int m, lapack_int n, double* a,
                                lapack_int lda, double* sva, double* u,
                                lapack_int ldu, double* v, lapack_int ldv,
                                double* work, lapack_int lwork,
                                lapack_int* iwork );
lapack_int LAPACKE_zgejsv_work( int matrix_layout, char joba, char jobu,
                                char jobv, char jobr, char jobt, char jobp,
                                lapack_int m, lapack_int n, lapack_complex_double* a,
                                lapack_int lda, double* sva, lapack_complex_double* u,
                                lapack_int ldu, lapack_complex_double* v, lapack_int ldv,
                                lapack_complex_double* cwork, lapack_int lwork,
                                double* work, lapack_int lrwork,
                                lapack_int* iwork );

lapack_int LAPACKE_dgelqf_work( int matrix_layout, lapack_int m, lapack_int n,
                                double* a, lapack_int lda, double* tau,
                                double* work, lapack_int lwork );
lapack_int LAPACKE_zgelqf_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_complex_double* tau,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_dgels_work( int matrix_layout, char trans, lapack_int m,
                               lapack_int n, lapack_int nrhs, double* a,
                               lapack_int lda, double* b, lapack_int ldb,
                               double* work, lapack_int lwork );
lapack_int LAPACKE_zgels_work( int matrix_layout, char trans, lapack_int m,
                               lapack_int n, lapack_int nrhs,
                               lapack_complex_double* a, lapack_int lda,
                               lapack_complex_double* b, lapack_int ldb,
                               lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_dgelsd_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int nrhs, double* a, lapack_int lda,
                                double* b, lapack_int ldb, double* s,
                                double rcond, lapack_int* rank, double* work,
                                lapack_int lwork, lapack_int* iwork );
lapack_int LAPACKE_zgelsd_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int nrhs, lapack_complex_double* a,
                                lapack_int lda, lapack_complex_double* b,
                                lapack_int ldb, double* s, double rcond,
                                lapack_int* rank, lapack_complex_double* work,
                                lapack_int lwork, double* rwork,
                                lapack_int* iwork );

lapack_int LAPACKE_dgelss_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int nrhs, double* a, lapack_int lda,
                                double* b, lapack_int ldb, double* s,
                                double rcond, lapack_int* rank, double* work,
                                lapack_int lwork );
lapack_int LAPACKE_zgelss_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int nrhs, lapack_complex_double* a,
                                lapack_int lda, lapack_complex_double* b,
                                lapack_int ldb, double* s, double rcond,
                                lapack_int* rank, lapack_complex_double* work,
                                lapack_int lwork, double* rwork );

lapack_int LAPACKE_dgelsy_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int nrhs, double* a, lapack_int lda,
                                double* b, lapack_int ldb, lapack_int* jpvt,
                                double rcond, lapack_int* rank, double* work,
                                lapack_int lwork );
lapack_int LAPACKE_zgelsy_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int nrhs, lapack_complex_double* a,
                                lapack_int lda, lapack_complex_double* b,
                                lapack_int ldb, lapack_int* jpvt, double rcond,
                                lapack_int* rank, lapack_complex_double* work,
                                lapack_int lwork, double* rwork );

lapack_int LAPACKE_dgeqp3_work( int matrix_layout, lapack_int m, lapack_int n,
                                double* a, lapack_int lda, lapack_int* jpvt,
                                double* tau, double* work, lapack_int lwork );
lapack_int LAPACKE_zgeqp3_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_int* jpvt, lapack_complex_double* tau,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork );

lapack_int LAPACKE_dgeqrf_work( int matrix_layout, lapack_int m, lapack_int n,
                                double* a, lapack_int lda, double* tau,
                                double* work, lapack_int lwork );
lapack_int LAPACKE_zgeqrf_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_complex_double* tau,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_dgesdd_work( int matrix_layout, char jobz, lapack_int m,
                                lapack_int n, double* a, lapack_int lda,
                                double* s, double* u, lapack_int ldu,
                                double* vt, lapack_int ldvt, double* work,
                                lapack_int lwork, lapack_int* iwork );
lapack_int LAPACKE_zgesdd_work( int matrix_layout, char jobz, lapack_int m,
                                lapack_int n, lapack_complex_double* a,
                                lapack_int lda, double* s,
                                lapack_complex_double* u, lapack_int ldu,
                                lapack_complex_double* vt, lapack_int ldvt,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork, lapack_int* iwork );

lapack_int LAPACKE_dgesv_work( int matrix_layout, lapack_int n, lapack_int nrhs,
                               double* a, lapack_int lda, lapack_int* ipiv,
                               double* b, lapack_int ldb );
lapack_int LAPACKE_zgesv_work( int matrix_layout, lapack_int n, lapack_int nrhs,
                               lapack_complex_double* a, lapack_int lda,
                               lapack_int* ipiv, lapack_complex_double* b,
                               lapack_int ldb );

lapack_int LAPACKE_dsgesv_work( int matrix_layout, lapack_int n, lapack_int nrhs,
                                double* a, lapack_int lda, lapack_int* ipiv,
                                double* b, lapack_int ldb, double* x,
                                lapack_int ldx, double* work, float* swork,
                                lapack_int* iter );
lapack_int LAPACKE_zcgesv_work( int matrix_layout, lapack_int n, lapack_int nrhs,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_int* ipiv, lapack_complex_double* b,
                                lapack_int ldb, lapack_complex_double* x,
                                lapack_int ldx, lapack_complex_double* work,
                                lapack_complex_float* swork, double* rwork,
                                lapack_int* iter );

lapack_int LAPACKE_dgesvd_work( int matrix_layout, char jobu, char jobvt,
                                lapack_int m, lapack_int n, double* a,
                                lapack_int lda, double* s, double* u,
                                lapack_int ldu, double* vt, lapack_int ldvt,
                                double* work, lapack_int lwork );
lapack_int LAPACKE_zgesvd_work( int matrix_layout, char jobu, char jobvt,
                                lapack_int m, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                double* s, lapack_complex_double* u,
                                lapack_int ldu, lapack_complex_double* vt,
                                lapack_int ldvt, lapack_complex_double* work,
                                lapack_int lwork, double* rwork );

lapack_int LAPACKE_dgesvdq_work( int matrix_layout, char joba, char jobp,
                                char jobr, char jobu, char jobv,
                                lapack_int m, lapack_int n, double* a,
                                lapack_int lda, double* s, double* u,
                                lapack_int ldu, double* v, lapack_int ldv,
                                lapack_int* numrank,
                                lapack_int* iwork, lapack_int liwork,
                                double* work, lapack_int lwork,
                                double* rwork, lapack_int lrwork);
lapack_int LAPACKE_zgesvdq_work( int matrix_layout, char joba, char jobp,
                                char jobr, char jobu, char jobv,
                                lapack_int m, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                double* s, lapack_complex_double* u,
                                lapack_int ldu, lapack_complex_double* v,
                                lapack_int ldv, lapack_int* numrank,
                                lapack_int* iwork, lapack_int liwork,
                                lapack_complex_double* cwork, lapack_int lcwork,
                                double* rwork, lapack_int lrwork);

lapack_int LAPACKE_dgesvdx_work( int matrix_layout, char jobu, char jobvt, char range,
                                 lapack_int m, lapack_int n, double* a,
                                 lapack_int lda, double vl, double vu,
                                 lapack_int il, lapack_int iu, lapack_int* ns,
                                 double* s, double* u, lapack_int ldu,
                                 double* vt, lapack_int ldvt,
                                 double* work, lapack_int lwork, lapack_int* iwork );
lapack_int LAPACKE_zgesvdx_work( int matrix_layout, char jobu, char jobvt, char range,
                                 lapack_int m, lapack_int n, lapack_complex_double* a,
                                 lapack_int lda, double vl, double vu,
                                 lapack_int il, lapack_int iu, lapack_int* ns,
                                 double* s, lapack_complex_double* u, lapack_int ldu,
                                 lapack_complex_double* vt, lapack_int ldvt,
                                 lapack_complex_double* work, lapack_int lwork,
                                 double* rwork, lapack_int* iwork );

lapack_int LAPACKE_dgesvx_work( int matrix_layout, char fact, char trans,
                                lapack_int n, lapack_int nrhs, double* a,
                                lapack_int lda, double* af, lapack_int ldaf,
                                lapack_int* ipiv, char* equed, double* r,
                                double* c, double* b, lapack_int ldb, double* x,
                                lapack_int ldx, double* rcond, double* ferr,
                                double* berr, double* work, lapack_int* iwork );
lapack_int LAPACKE_zgesvx_work( int matrix_layout, char fact, char trans,
                                lapack_int n, lapack_int nrhs,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_complex_double* af, lapack_int ldaf,
                                lapack_int* ipiv, char* equed, double* r,
                                double* c, lapack_complex_double* b,
                                lapack_int ldb, lapack_complex_double* x,
                                lapack_int ldx, double* rcond, double* ferr,
                                double* berr, lapack_complex_double* work,
                                double* rwork );

lapack_int LAPACKE_dgetrf_work( int matrix_layout, lapack_int m, lapack_int n,
                                double* a, lapack_int lda, lapack_int* ipiv );
lapack_int LAPACKE_zgetrf_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_int* ipiv );

lapack_int LAPACKE_dgetri_work( int matrix_layout, lapack_int n, double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                double* work, lapack_int lwork );
lapack_int LAPACKE_zgetri_work( int matrix_layout, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                const lapack_int* ipiv,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_dgetrs_work( int matrix_layout, char trans, lapack_int n,
                                lapack_int nrhs, const double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                double* b, lapack_int ldb );
lapack_int LAPACKE_zgetrs_work( int matrix_layout, char trans, lapack_int n,
                                lapack_int nrhs, const lapack_complex_double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dgetsls_work( int matrix_layout, char trans, lapack_int m,
                                 lapack_int n, lapack_int nrhs, double* a,
                                 lapack_int lda, double* b, lapack_int ldb,
                                 double* work, lapack_int lwork );
lapack_int LAPACKE_zgetsls_work( int matrix_layout, char trans, lapack_int m,
                                 lapack_int n, lapack_int nrhs,
                                 lapack_complex_double* a, lapack_int lda,
                                 lapack_complex_double* b, lapack_int ldb,
                                 lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_dgges_work( int matrix_layout, char jobvsl, char jobvsr,
                               char sort, LAPACK_D_SELECT3 selctg, lapack_int n,
                               double* a, lapack_int lda, double* b,
                               lapack_int ldb, lapack_int* sdim, double* alphar,
                               double* alphai, double* beta, double* vsl,
                               lapack_int ldvsl, double* vsr, lapack_int ldvsr,
                               double* work, lapack_int lwork,
                               lapack_logical* bwork );
lapack_int LAPACKE_zgges_work( int matrix_layout, char jobvsl, char jobvsr,
                               char sort, LAPACK_Z_SELECT2 selctg, lapack_int n,
                               lapack_complex_double* a, lapack_int lda,
                               lapack_complex_double* b, lapack_int ldb,
                               lapack_int* sdim, lapack_complex_double* alpha,
                               lapack_complex_double* beta,
                               lapack_complex_double* vsl, lapack_int ldvsl,
                               lapack_complex_double* vsr, lapack_int ldvsr,
                               lapack_complex_double* work, lapack_int lwork,
                               double* rwork, lapack_logical* bwork );

lapack_int LAPACKE_dggesx_work( int matrix_layout, char jobvsl, char jobvsr,
                                char sort, LAPACK_D_SELECT3 selctg, char sense,
                                lapack_int n, double* a, lapack_int lda,
                                double* b, lapack_int ldb, lapack_int* sdim,
                                double* alphar, double* alphai, double* beta,
                                double* vsl, lapack_int ldvsl, double* vsr,
                                lapack_int ldvsr, double* rconde,
                                double* rcondv, double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork,
                                lapack_logical* bwork );
lapack_int LAPACKE_zggesx_work( int matrix_layout, char jobvsl, char jobvsr,
                                char sort, LAPACK_Z_SELECT2 selctg, char sense,
                                lapack_int n, lapack_complex_double* a,
                                lapack_int lda, lapack_complex_double* b,
                                lapack_int ldb, lapack_int* sdim,
                                lapack_complex_double* alpha,
                                lapack_complex_double* beta,
                                lapack_complex_double* vsl, lapack_int ldvsl,
                                lapack_complex_double* vsr, lapack_int ldvsr,
                                double* rconde, double* rcondv,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork, lapack_int* iwork,
                                lapack_int liwork, lapack_logical* bwork );

lapack_int LAPACKE_dggev_work( int matrix_layout, char jobvl, char jobvr,
                               lapack_int n, double* a, lapack_int lda,
                               double* b, lapack_int ldb, double* alphar,
                               double* alphai, double* beta, double* vl,
                               lapack_int ldvl, double* vr, lapack_int ldvr,
                               double* work, lapack_int lwork );
lapack_int LAPACKE_zggev_work( int matrix_layout, char jobvl, char jobvr,
                               lapack_int n, lapack_complex_double* a,
                               lapack_int lda, lapack_complex_double* b,
                               lapack_int ldb, lapack_complex_double* alpha,
                               lapack_complex_double* beta,
                               lapack_complex_double* vl, lapack_int ldvl,
                               lapack_complex_double* vr, lapack_int ldvr,
                               lapack_complex_double* work, lapack_int lwork,
                               double* rwork );

lapack_int LAPACKE_dggevx_work( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, lapack_int n, double* a,
                                lapack_int lda, double* b, lapack_int ldb,
                                double* alphar, double* alphai, double* beta,
                                double* vl, lapack_int ldvl, double* vr,
                                lapack_int ldvr, lapack_int* ilo,
                                lapack_int* ihi, double* lscale, double* rscale,
                                double* abnrm, double* bbnrm, double* rconde,
                                double* rcondv, double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_logical* bwork );
lapack_int LAPACKE_zggevx_work( int matrix_layout, char balanc, char jobvl,
                                char jobvr, char sense, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_complex_double* b, lapack_int ldb,
                                lapack_complex_double* alpha,
                                lapack_complex_double* beta,
                                lapack_complex_double* vl, lapack_int ldvl,
                                lapack_complex_double* vr, lapack_int ldvr,
                                lapack_int* ilo, lapack_int* ihi,
                                double* lscale, double* rscale, double* abnrm,
                                double* bbnrm, double* rconde, double* rcondv,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork, lapack_int* iwork,
                                lapack_logical* bwork );

lapack_int LAPACKE_dggglm_work( int matrix_layout, lapack_int n, lapack_int m,
                                lapack_int p, double* a, lapack_int lda,
                                double* b, lapack_int ldb, double* d, double* x,
                                double* y, double* work, lapack_int lwork );
lapack_int LAPACKE_zggglm_work( int matrix_layout, lapack_int n, lapack_int m,
                                lapack_int p, lapack_complex_double* a,
                                lapack_int lda, lapack_complex_double* b,
                                lapack_int ldb, lapack_complex_double* d,
                                lapack_complex_double* x,
                                lapack_complex_double* y,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_dgglse_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int p, double* a, lapack_int lda,
                                double* b, lapack_int ldb, double* c, double* d,
                                double* x, double* work, lapack_int lwork );
lapack_int LAPACKE_zgglse_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int p, lapack_complex_double* a,
                                lapack_int lda, lapack_complex_double* b,
                                lapack_int ldb, lapack_complex_double* c,
                                lapack_complex_double* d,
                                lapack_complex_double* x,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_dggsvd3_work( int matrix_layout, char jobu, char jobv,
                                 char jobq, lapack_int m, lapack_int n,
                                 lapack_int p, lapack_int* k, lapack_int* l,
                                 double* a, lapack_int lda, double* b,
                                 lapack_int ldb, double* alpha, double* beta,
                                 double* u, lapack_int ldu, double* v,
                                 lapack_int ldv, double* q, lapack_int ldq,
                                 double* work, lapack_int lwork,
                                 lapack_int* iwork );
lapack_int LAPACKE_zggsvd3_work( int matrix_layout, char jobu, char jobv,
                                 char jobq, lapack_int m, lapack_int n,
                                 lapack_int p, lapack_int* k, lapack_int* l,
                                 lapack_complex_double* a, lapack_int lda,
                                 lapack_complex_double* b, lapack_int ldb,
                                 double* alpha, double* beta,
                                 lapack_complex_double* u, lapack_int ldu,
                                 lapack_complex_double* v, lapack_int ldv,
                                 lapack_complex_double* q, lapack_int ldq,
                                 lapack_complex_double* work, lapack_int lwork,
                                 double* rwork, lapack_int* iwork );

lapack_int LAPACKE_dgtcon_work( char norm, lapack_int n, const double* dl,
                                const double* d, const double* du,
                                const double* du2, const lapack_int* ipiv,
                                double anorm, double* rcond, double* work,
                                lapack_int* iwork );
lapack_int LAPACKE_zgtcon_work( char norm, lapack_int n,
                                const lapack_complex_double* dl,
                                const lapack_complex_double* d,
                                const lapack_complex_double* du,
                                const lapack_complex_double* du2,
                                const lapack_int* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );

lapack_int LAPACKE_dgtsv_work( int matrix_layout, lapack_int n, lapack_int nrhs,
                               double* dl, double* d, double* du, double* b,
                               lapack_int ldb );
lapack_int LAPACKE_zgtsv_work( int matrix_layout, lapack_int n, lapack_int nrhs,
                               lapack_complex_double* dl,
                               lapack_complex_double* d,
                               lapack_complex_double* du,
                               lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dgtsvx_work( int matrix_layout, char fact, char trans,
                                lapack_int n, lapack_int nrhs, const double* dl,
                                const double* d, const double* du, double* dlf,
                                double* df, double* duf, double* du2,
                                lapack_int* ipiv, const double* b,
                                lapack_int ldb, double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, lapack_int* iwork );
lapack_int LAPACKE_zgtsvx_work( int matrix_layout, char fact, char trans,
                                lapack_int n, lapack_int nrhs,
                                const lapack_complex_double* dl,
                                const lapack_complex_double* d,
                                const lapack_complex_double* du,
                                lapack_complex_double* dlf,
                                lapack_complex_double* df,
                                lapack_complex_double* duf,
                                lapack_complex_double* du2, lapack_int* ipiv,
                                const lapack_complex_double* b, lapack_int ldb,
                                lapack_complex_double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_dgttrf_work( lapack_int n, double* dl, double* d, double* du,
                                double* du2, lapack_int* ipiv );
lapack_int LAPACKE_zgttrf_work( lapack_int n, lapack_complex_double* dl,
                                lapack_complex_double* d,
                                lapack_complex_double* du,
                                lapack_complex_double* du2, lapack_int* ipiv );

lapack_int LAPACKE_dgttrs_work( int matrix_layout, char trans, lapack_int n,
                                lapack_int nrhs, const double* dl,
                                const double* d, const double* du,
                                const double* du2, const lapack_int* ipiv,
                                double* b, lapack_int ldb );
lapack_int LAPACKE_zgttrs_work( int matrix_layout, char trans, lapack_int n,
                                lapack_int nrhs,
                                const lapack_complex_double* dl,
                                const lapack_complex_double* d,
                                const lapack_complex_double* du,
                                const lapack_complex_double* du2,
                                const lapack_int* ipiv,
                                lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_zhbev_work( int matrix_layout, char jobz, char uplo,
                               lapack_int n, lapack_int kd,
                               lapack_complex_double* ab, lapack_int ldab,
                               double* w, lapack_complex_double* z,
                               lapack_int ldz, lapack_complex_double* work,
                               double* rwork );

lapack_int LAPACKE_zhbevd_work( int matrix_layout, char jobz, char uplo,
                                lapack_int n, lapack_int kd,
                                lapack_complex_double* ab, lapack_int ldab,
                                double* w, lapack_complex_double* z,
                                lapack_int ldz, lapack_complex_double* work,
                                lapack_int lwork, double* rwork,
                                lapack_int lrwork, lapack_int* iwork,
                                lapack_int liwork );

lapack_int LAPACKE_zhbevx_work( int matrix_layout, char jobz, char range,
                                char uplo, lapack_int n, lapack_int kd,
                                lapack_complex_double* ab, lapack_int ldab,
                                lapack_complex_double* q, lapack_int ldq,
                                double vl, double vu, lapack_int il,
                                lapack_int iu, double abstol, lapack_int* m,
                                double* w, lapack_complex_double* z,
                                lapack_int ldz, lapack_complex_double* work,
                                double* rwork, lapack_int* iwork,
                                lapack_int* ifail );

lapack_int LAPACKE_zhbgv_work( int matrix_layout, char jobz, char uplo,
                               lapack_int n, lapack_int ka, lapack_int kb,
                               lapack_complex_double* ab, lapack_int ldab,
                               lapack_complex_double* bb, lapack_int ldbb,
                               double* w, lapack_complex_double* z,
                               lapack_int ldz, lapack_complex_double* work,
                               double* rwork );

lapack_int LAPACKE_zhbgvd_work( int matrix_layout, char jobz, char uplo,
                                lapack_int n, lapack_int ka, lapack_int kb,
                                lapack_complex_double* ab, lapack_int ldab,
                                lapack_complex_double* bb, lapack_int ldbb,
                                double* w, lapack_complex_double* z,
                                lapack_int ldz, lapack_complex_double* work,
                                lapack_int lwork, double* rwork,
                                lapack_int lrwork, lapack_int* iwork,
                                lapack_int liwork );

lapack_int LAPACKE_zhbgvx_work( int matrix_layout, char jobz, char range,
                                char uplo, lapack_int n, lapack_int ka,
                                lapack_int kb, lapack_complex_double* ab,
                                lapack_int ldab, lapack_complex_double* bb,
                                lapack_int ldbb, lapack_complex_double* q,
                                lapack_int ldq, double vl, double vu,
                                lapack_int il, lapack_int iu, double abstol,
                                lapack_int* m, double* w,
                                lapack_complex_double* z, lapack_int ldz,
                                lapack_complex_double* work, double* rwork,
                                lapack_int* iwork, lapack_int* ifail );

lapack_int LAPACKE_zhbtrd_work( int matrix_layout, char vect, char uplo,
                                lapack_int n, lapack_int kd,
                                lapack_complex_double* ab, lapack_int ldab,
                                double* d, double* e, lapack_complex_double* q,
                                lapack_int ldq, lapack_complex_double* work );

lapack_int LAPACKE_zhecon_work( int matrix_layout, char uplo, lapack_int n,
                                const lapack_complex_double* a, lapack_int lda,
                                const lapack_int* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );

lapack_int LAPACKE_zheev_work( int matrix_layout, char jobz, char uplo,
                               lapack_int n, lapack_complex_double* a,
                               lapack_int lda, double* w,
                               lapack_complex_double* work, lapack_int lwork,
                               double* rwork );

lapack_int LAPACKE_zheevd_work( int matrix_layout, char jobz, char uplo,
                                lapack_int n, lapack_complex_double* a,
                                lapack_int lda, double* w,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork, lapack_int lrwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_zheevr_work( int matrix_layout, char jobz, char range,
                                char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                double vl, double vu, lapack_int il,
                                lapack_int iu, double abstol, lapack_int* m,
                                double* w, lapack_complex_double* z,
                                lapack_int ldz, lapack_int* isuppz,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork, lapack_int lrwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_zheevx_work( int matrix_layout, char jobz, char range,
                                char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                double vl, double vu, lapack_int il,
                                lapack_int iu, double abstol, lapack_int* m,
                                double* w, lapack_complex_double* z,
                                lapack_int ldz, lapack_complex_double* work,
                                lapack_int lwork, double* rwork,
                                lapack_int* iwork, lapack_int* ifail );

lapack_int LAPACKE_zhegv_work( int matrix_layout, lapack_int itype, char jobz,
                               char uplo, lapack_int n,
                               lapack_complex_double* a, lapack_int lda,
                               lapack_complex_double* b, lapack_int ldb,
                               double* w, lapack_complex_double* work,
                               lapack_int lwork, double* rwork );

lapack_int LAPACKE_zhegvd_work( int matrix_layout, lapack_int itype, char jobz,
                                char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_complex_double* b, lapack_int ldb,
                                double* w, lapack_complex_double* work,
                                lapack_int lwork, double* rwork,
                                lapack_int lrwork, lapack_int* iwork,
                                lapack_int liwork );

lapack_int LAPACKE_zhegvx_work( int matrix_layout, lapack_int itype, char jobz,
                                char range, char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_complex_double* b, lapack_int ldb,
                                double vl, double vu, lapack_int il,
                                lapack_int iu, double abstol, lapack_int* m,
                                double* w, lapack_complex_double* z,
                                lapack_int ldz, lapack_complex_double* work,
                                lapack_int lwork, double* rwork,
                                lapack_int* iwork, lapack_int* ifail );

lapack_int LAPACKE_zhesv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int nrhs, lapack_complex_double* a,
                               lapack_int lda, lapack_int* ipiv,
                               lapack_complex_double* b, lapack_int ldb,
                               lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_zhesvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int nrhs,
                                const lapack_complex_double* a, lapack_int lda,
                                lapack_complex_double* af, lapack_int ldaf,
                                lapack_int* ipiv,
                                const lapack_complex_double* b, lapack_int ldb,
                                lapack_complex_double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork );

lapack_int LAPACKE_zhetrd_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                double* d, double* e,
                                lapack_complex_double* tau,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_zhetrf_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_int* ipiv, lapack_complex_double* work,
                                lapack_int lwork );

lapack_int LAPACKE_zhetri_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                const lapack_int* ipiv,
                                lapack_complex_double* work );

lapack_int LAPACKE_zhetrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs, const lapack_complex_double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_zhpcon_work( int matrix_layout, char uplo, lapack_int n,
                                const lapack_complex_double* ap,
                                const lapack_int* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );

lapack_int LAPACKE_zhpev_work( int matrix_layout, char jobz, char uplo,
                               lapack_int n, lapack_complex_double* ap,
                               double* w, lapack_complex_double* z,
                               lapack_int ldz, lapack_complex_double* work,
                               double* rwork );

lapack_int LAPACKE_zhpevd_work( int matrix_layout, char jobz, char uplo,
                                lapack_int n, lapack_complex_double* ap,
                                double* w, lapack_complex_double* z,
                                lapack_int ldz, lapack_complex_double* work,
                                lapack_int lwork, double* rwork,
                                lapack_int lrwork, lapack_int* iwork,
                                lapack_int liwork );

lapack_int LAPACKE_zhpevx_work( int matrix_layout, char jobz, char range,
                                char uplo, lapack_int n,
                                lapack_complex_double* ap, double vl, double vu,
                                lapack_int il, lapack_int iu, double abstol,
                                lapack_int* m, double* w,
                                lapack_complex_double* z, lapack_int ldz,
                                lapack_complex_double* work, double* rwork,
                                lapack_int* iwork, lapack_int* ifail );

lapack_int LAPACKE_zhpgv_work( int matrix_layout, lapack_int itype, char jobz,
                               char uplo, lapack_int n,
                               lapack_complex_double* ap,
                               lapack_complex_double* bp, double* w,
                               lapack_complex_double* z, lapack_int ldz,
                               lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_zhpgvd_work( int matrix_layout, lapack_int itype, char jobz,
                                char uplo, lapack_int n,
                                lapack_complex_double* ap,
                                lapack_complex_double* bp, double* w,
                                lapack_complex_double* z, lapack_int ldz,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork, lapack_int lrwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_zhpgvx_work( int matrix_layout, lapack_int itype, char jobz,
                                char range, char uplo, lapack_int n,
                                lapack_complex_double* ap,
                                lapack_complex_double* bp, double vl, double vu,
                                lapack_int il, lapack_int iu, double abstol,
                                lapack_int* m, double* w,
                                lapack_complex_double* z, lapack_int ldz,
                                lapack_complex_double* work, double* rwork,
                                lapack_int* iwork, lapack_int* ifail );

lapack_int LAPACKE_zhpsv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int nrhs, lapack_complex_double* ap,
                               lapack_int* ipiv, lapack_complex_double* b,
                               lapack_int ldb );

lapack_int LAPACKE_zhpsvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int nrhs,
                                const lapack_complex_double* ap,
                                lapack_complex_double* afp, lapack_int* ipiv,
                                const lapack_complex_double* b, lapack_int ldb,
                                lapack_complex_double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_zhptrd_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* ap, double* d, double* e,
                                lapack_complex_double* tau );

lapack_int LAPACKE_zhptrf_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* ap, lapack_int* ipiv );

lapack_int LAPACKE_zhptri_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* ap,
                                const lapack_int* ipiv,
                                lapack_complex_double* work );

lapack_int LAPACKE_zhptrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs,
                                const lapack_complex_double* ap,
                                const lapack_int* ipiv,
                                lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dhsein_work( int matrix_layout, char job, char eigsrc,
                                char initv, lapack_logical* select,
                                lapack_int n, const double* h, lapack_int ldh,
                                double* wr, const double* wi, double* vl,
                                lapack_int ldvl, double* vr, lapack_int ldvr,
                                lapack_int mm, lapack_int* m, double* work,
                                lapack_int* ifaill, lapack_int* ifailr );
lapack_int LAPACKE_zhsein_work( int matrix_layout, char job, char eigsrc,
                                char initv, const lapack_logical* select,
                                lapack_int n, const lapack_complex_double* h,
                                lapack_int ldh, lapack_complex_double* w,
                                lapack_complex_double* vl, lapack_int ldvl,
                                lapack_complex_double* vr, lapack_int ldvr,
                                lapack_int mm, lapack_int* m,
                                lapack_complex_double* work, double* rwork,
                                lapack_int* ifaill, lapack_int* ifailr );

lapack_int LAPACKE_dhseqr_work( int matrix_layout, char job, char compz,
                                lapack_int n, lapack_int ilo, lapack_int ihi,
                                double* h, lapack_int ldh, double* wr,
                                double* wi, double* z, lapack_int ldz,
                                double* work, lapack_int lwork );
lapack_int LAPACKE_zhseqr_work( int matrix_layout, char job, char compz,
                                lapack_int n, lapack_int ilo, lapack_int ihi,
                                lapack_complex_double* h, lapack_int ldh,
                                lapack_complex_double* w,
                                lapack_complex_double* z, lapack_int ldz,
                                lapack_complex_double* work, lapack_int lwork );

float LAPACKE_slamch_work( char cmach );
double LAPACKE_dlamch_work( char cmach );

double LAPACKE_dlange_work( int matrix_layout, char norm, lapack_int m,
                                lapack_int n, const double* a, lapack_int lda,
                                double* work );
double LAPACKE_zlange_work( int matrix_layout, char norm, lapack_int m,
                                lapack_int n, const lapack_complex_double* a,
                                lapack_int lda, double* work );

double LAPACKE_zlanhe_work( int matrix_layout, char norm, char uplo,
                                lapack_int n, const lapack_complex_double* a,
                                lapack_int lda, double* work );

double LAPACKE_dlansy_work( int matrix_layout, char norm, char uplo,
                                lapack_int n, const double* a, lapack_int lda,
                                double* work );
double LAPACKE_zlansy_work( int matrix_layout, char norm, char uplo,
                                lapack_int n, const lapack_complex_double* a,
                                lapack_int lda, double* work );

double LAPACKE_dlantr_work( int matrix_layout, char norm, char uplo,
                                char diag, lapack_int m, lapack_int n,
                                const double* a, lapack_int lda, double* work );
double LAPACKE_zlantr_work( int matrix_layout, char norm, char uplo,
                                char diag, lapack_int m, lapack_int n,
                                const lapack_complex_double* a, lapack_int lda,
                                double* work );

lapack_int LAPACKE_dlatms_work( int matrix_layout, lapack_int m, lapack_int n,
                                char dist, lapack_int* iseed, char sym,
                                double* d, lapack_int mode, double cond,
                                double dmax, lapack_int kl, lapack_int ku,
                                char pack, double* a, lapack_int lda,
                                double* work );
lapack_int LAPACKE_zlatms_work( int matrix_layout, lapack_int m, lapack_int n,
                                char dist, lapack_int* iseed, char sym,
                                double* d, lapack_int mode, double cond,
                                double dmax, lapack_int kl, lapack_int ku,
                                char pack, lapack_complex_double* a,
                                lapack_int lda, lapack_complex_double* work );

lapack_int LAPACKE_dopgtr_work( int matrix_layout, char uplo, lapack_int n,
                                const double* ap, const double* tau, double* q,
                                lapack_int ldq, double* work );

lapack_int LAPACKE_dopmtr_work( int matrix_layout, char side, char uplo,
                                char trans, lapack_int m, lapack_int n,
                                const double* ap, const double* tau, double* c,
                                lapack_int ldc, double* work );

lapack_int LAPACKE_dorghr_work( int matrix_layout, lapack_int n, lapack_int ilo,
                                lapack_int ihi, double* a, lapack_int lda,
                                const double* tau, double* work,
                                lapack_int lwork );

lapack_int LAPACKE_dorglq_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int k, double* a, lapack_int lda,
                                const double* tau, double* work,
                                lapack_int lwork );

lapack_int LAPACKE_dorgqr_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int k, double* a, lapack_int lda,
                                const double* tau, double* work,
                                lapack_int lwork );

lapack_int LAPACKE_dorgtr_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda, const double* tau,
                                double* work, lapack_int lwork );

lapack_int LAPACKE_dormhr_work( int matrix_layout, char side, char trans,
                                lapack_int m, lapack_int n, lapack_int ilo,
                                lapack_int ihi, const double* a, lapack_int lda,
                                const double* tau, double* c, lapack_int ldc,
                                double* work, lapack_int lwork );

lapack_int LAPACKE_dormlq_work( int matrix_layout, char side, char trans,
                                lapack_int m, lapack_int n, lapack_int k,
                                const double* a, lapack_int lda,
                                const double* tau, double* c, lapack_int ldc,
                                double* work, lapack_int lwork );

lapack_int LAPACKE_dormqr_work( int matrix_layout, char side, char trans,
                                lapack_int m, lapack_int n, lapack_int k,
                                const double* a, lapack_int lda,
                                const double* tau, double* c, lapack_int ldc,
                                double* work, lapack_int lwork );

lapack_int LAPACKE_dormtr_work( int matrix_layout, char side, char uplo,
                                char trans, lapack_int m, lapack_int n,
                                const double* a, lapack_int lda,
                                const double* tau, double* c, lapack_int ldc,
                                double* work, lapack_int lwork );

lapack_int LAPACKE_dpbcon_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int kd, const double* ab,
                                lapack_int ldab, double anorm, double* rcond,
                                double* work, lapack_int* iwork );
lapack_int LAPACKE_zpbcon_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int kd, const lapack_complex_double* ab,
                                lapack_int ldab, double anorm, double* rcond,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_dpbsv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int kd, lapack_int nrhs, double* ab,
                               lapack_int ldab, double* b, lapack_int ldb );
lapack_int LAPACKE_zpbsv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int kd, lapack_int nrhs,
                               lapack_complex_double* ab, lapack_int ldab,
                               lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dpbsvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int kd, lapack_int nrhs,
                                double* ab, lapack_int ldab, double* afb,
                                lapack_int ldafb, char* equed, double* s,
                                double* b, lapack_int ldb, double* x,
                                lapack_int ldx, double* rcond, double* ferr,
                                double* berr, double* work, lapack_int* iwork );
lapack_int LAPACKE_zpbsvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int kd, lapack_int nrhs,
                                lapack_complex_double* ab, lapack_int ldab,
                                lapack_complex_double* afb, lapack_int ldafb,
                                char* equed, double* s,
                                lapack_complex_double* b, lapack_int ldb,
                                lapack_complex_double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_dpbtrf_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int kd, double* ab, lapack_int ldab );
lapack_int LAPACKE_zpbtrf_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int kd, lapack_complex_double* ab,
                                lapack_int ldab );

lapack_int LAPACKE_dpbtrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int kd, lapack_int nrhs,
                                const double* ab, lapack_int ldab, double* b,
                                lapack_int ldb );
lapack_int LAPACKE_zpbtrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int kd, lapack_int nrhs,
                                const lapack_complex_double* ab,
                                lapack_int ldab, lapack_complex_double* b,
                                lapack_int ldb );

lapack_int LAPACKE_dpocon_work( int matrix_layout, char uplo, lapack_int n,
                                const double* a, lapack_int lda, double anorm,
                                double* rcond, double* work,
                                lapack_int* iwork );
lapack_int LAPACKE_zpocon_work( int matrix_layout, char uplo, lapack_int n,
                                const lapack_complex_double* a, lapack_int lda,
                                double anorm, double* rcond,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_dposv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int nrhs, double* a, lapack_int lda,
                               double* b, lapack_int ldb );
lapack_int LAPACKE_zposv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int nrhs, lapack_complex_double* a,
                               lapack_int lda, lapack_complex_double* b,
                               lapack_int ldb );
lapack_int LAPACKE_dsposv_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs, double* a, lapack_int lda,
                                double* b, lapack_int ldb, double* x,
                                lapack_int ldx, double* work, float* swork,
                                lapack_int* iter );
lapack_int LAPACKE_zcposv_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs, lapack_complex_double* a,
                                lapack_int lda, lapack_complex_double* b,
                                lapack_int ldb, lapack_complex_double* x,
                                lapack_int ldx, lapack_complex_double* work,
                                lapack_complex_float* swork, double* rwork,
                                lapack_int* iter );

lapack_int LAPACKE_dposvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int nrhs, double* a,
                                lapack_int lda, double* af, lapack_int ldaf,
                                char* equed, double* s, double* b,
                                lapack_int ldb, double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, lapack_int* iwork );
lapack_int LAPACKE_zposvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int nrhs,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_complex_double* af, lapack_int ldaf,
                                char* equed, double* s,
                                lapack_complex_double* b, lapack_int ldb,
                                lapack_complex_double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_dpotrf_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda );
lapack_int LAPACKE_zpotrf_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda );

lapack_int LAPACKE_dpotri_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda );
lapack_int LAPACKE_zpotri_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda );

lapack_int LAPACKE_dpotrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs, const double* a,
                                lapack_int lda, double* b, lapack_int ldb );
lapack_int LAPACKE_zpotrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs, const lapack_complex_double* a,
                                lapack_int lda, lapack_complex_double* b,
                                lapack_int ldb );

lapack_int LAPACKE_dppcon_work( int matrix_layout, char uplo, lapack_int n,
                                const double* ap, double anorm, double* rcond,
                                double* work, lapack_int* iwork );
lapack_int LAPACKE_zppcon_work( int matrix_layout, char uplo, lapack_int n,
                                const lapack_complex_double* ap, double anorm,
                                double* rcond, lapack_complex_double* work,
                                double* rwork );

lapack_int LAPACKE_dppsv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int nrhs, double* ap, double* b,
                               lapack_int ldb );
lapack_int LAPACKE_zppsv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int nrhs, lapack_complex_double* ap,
                               lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dppsvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int nrhs, double* ap,
                                double* afp, char* equed, double* s, double* b,
                                lapack_int ldb, double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, lapack_int* iwork );
lapack_int LAPACKE_zppsvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int nrhs,
                                lapack_complex_double* ap,
                                lapack_complex_double* afp, char* equed,
                                double* s, lapack_complex_double* b,
                                lapack_int ldb, lapack_complex_double* x,
                                lapack_int ldx, double* rcond, double* ferr,
                                double* berr, lapack_complex_double* work,
                                double* rwork );

lapack_int LAPACKE_dpptrf_work( int matrix_layout, char uplo, lapack_int n,
                                double* ap );
lapack_int LAPACKE_zpptrf_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* ap );

lapack_int LAPACKE_dpptri_work( int matrix_layout, char uplo, lapack_int n,
                                double* ap );
lapack_int LAPACKE_zpptri_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* ap );

lapack_int LAPACKE_dpptrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs, const double* ap, double* b,
                                lapack_int ldb );
lapack_int LAPACKE_zpptrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs,
                                const lapack_complex_double* ap,
                                lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dptcon_work( lapack_int n, const double* d, const double* e,
                                double anorm, double* rcond, double* work );
lapack_int LAPACKE_zptcon_work( lapack_int n, const double* d,
                                const lapack_complex_double* e, double anorm,
                                double* rcond, double* work );

lapack_int LAPACKE_dpteqr_work( int matrix_layout, char compz, lapack_int n,
                                double* d, double* e, double* z, lapack_int ldz,
                                double* work );
lapack_int LAPACKE_zpteqr_work( int matrix_layout, char compz, lapack_int n,
                                double* d, double* e, lapack_complex_double* z,
                                lapack_int ldz, double* work );

lapack_int LAPACKE_dptsv_work( int matrix_layout, lapack_int n, lapack_int nrhs,
                               double* d, double* e, double* b,
                               lapack_int ldb );
lapack_int LAPACKE_zptsv_work( int matrix_layout, lapack_int n, lapack_int nrhs,
                               double* d, lapack_complex_double* e,
                               lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dptsvx_work( int matrix_layout, char fact, lapack_int n,
                                lapack_int nrhs, const double* d,
                                const double* e, double* df, double* ef,
                                const double* b, lapack_int ldb, double* x,
                                lapack_int ldx, double* rcond, double* ferr,
                                double* berr, double* work );
lapack_int LAPACKE_zptsvx_work( int matrix_layout, char fact, lapack_int n,
                                lapack_int nrhs, const double* d,
                                const lapack_complex_double* e, double* df,
                                lapack_complex_double* ef,
                                const lapack_complex_double* b, lapack_int ldb,
                                lapack_complex_double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_dpttrf_work( lapack_int n, double* d, double* e );
lapack_int LAPACKE_zpttrf_work( lapack_int n, double* d,
                                lapack_complex_double* e );

lapack_int LAPACKE_dpttrs_work( int matrix_layout, lapack_int n, lapack_int nrhs,
                                const double* d, const double* e, double* b,
                                lapack_int ldb );
lapack_int LAPACKE_zpttrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs, const double* d,
                                const lapack_complex_double* e,
                                lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dsbev_work( int matrix_layout, char jobz, char uplo,
                               lapack_int n, lapack_int kd, double* ab,
                               lapack_int ldab, double* w, double* z,
                               lapack_int ldz, double* work );

lapack_int LAPACKE_dsbevd_work( int matrix_layout, char jobz, char uplo,
                                lapack_int n, lapack_int kd, double* ab,
                                lapack_int ldab, double* w, double* z,
                                lapack_int ldz, double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_dsbevx_work( int matrix_layout, char jobz, char range,
                                char uplo, lapack_int n, lapack_int kd,
                                double* ab, lapack_int ldab, double* q,
                                lapack_int ldq, double vl, double vu,
                                lapack_int il, lapack_int iu, double abstol,
                                lapack_int* m, double* w, double* z,
                                lapack_int ldz, double* work,
                                lapack_int* iwork, lapack_int* ifail );

lapack_int LAPACKE_dsbgv_work( int matrix_layout, char jobz, char uplo,
                               lapack_int n, lapack_int ka, lapack_int kb,
                               double* ab, lapack_int ldab, double* bb,
                               lapack_int ldbb, double* w, double* z,
                               lapack_int ldz, double* work );

lapack_int LAPACKE_dsbgvd_work( int matrix_layout, char jobz, char uplo,
                                lapack_int n, lapack_int ka, lapack_int kb,
                                double* ab, lapack_int ldab, double* bb,
                                lapack_int ldbb, double* w, double* z,
                                lapack_int ldz, double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_dsbgvx_work( int matrix_layout, char jobz, char range,
                                char uplo, lapack_int n, lapack_int ka,
                                lapack_int kb, double* ab, lapack_int ldab,
                                double* bb, lapack_int ldbb, double* q,
                                lapack_int ldq, double vl, double vu,
                                lapack_int il, lapack_int iu, double abstol,
                                lapack_int* m, double* w, double* z,
                                lapack_int ldz, double* work, lapack_int* iwork,
                                lapack_int* ifail );

lapack_int LAPACKE_dsbtrd_work( int matrix_layout, char vect, char uplo,
                                lapack_int n, lapack_int kd, double* ab,
                                lapack_int ldab, double* d, double* e,
                                double* q, lapack_int ldq, double* work );

lapack_int LAPACKE_dspcon_work( int matrix_layout, char uplo, lapack_int n,
                                const double* ap, const lapack_int* ipiv,
                                double anorm, double* rcond, double* work,
                                lapack_int* iwork );
lapack_int LAPACKE_zspcon_work( int matrix_layout, char uplo, lapack_int n,
                                const lapack_complex_double* ap,
                                const lapack_int* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );

lapack_int LAPACKE_dspev_work( int matrix_layout, char jobz, char uplo,
                               lapack_int n, double* ap, double* w, double* z,
                               lapack_int ldz, double* work );

lapack_int LAPACKE_dspevd_work( int matrix_layout, char jobz, char uplo,
                                lapack_int n, double* ap, double* w, double* z,
                                lapack_int ldz, double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_dspevx_work( int matrix_layout, char jobz, char range,
                                char uplo, lapack_int n, double* ap, double vl,
                                double vu, lapack_int il, lapack_int iu,
                                double abstol, lapack_int* m, double* w,
                                double* z, lapack_int ldz, double* work,
                                lapack_int* iwork, lapack_int* ifail );

lapack_int LAPACKE_dspgv_work( int matrix_layout, lapack_int itype, char jobz,
                               char uplo, lapack_int n, double* ap, double* bp,
                               double* w, double* z, lapack_int ldz,
                               double* work );

lapack_int LAPACKE_dspgvd_work( int matrix_layout, lapack_int itype, char jobz,
                                char uplo, lapack_int n, double* ap, double* bp,
                                double* w, double* z, lapack_int ldz,
                                double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_dspgvx_work( int matrix_layout, lapack_int itype, char jobz,
                                char range, char uplo, lapack_int n, double* ap,
                                double* bp, double vl, double vu, lapack_int il,
                                lapack_int iu, double abstol, lapack_int* m,
                                double* w, double* z, lapack_int ldz,
                                double* work, lapack_int* iwork,
                                lapack_int* ifail );

lapack_int LAPACKE_dspsv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int nrhs, double* ap, lapack_int* ipiv,
                               double* b, lapack_int ldb );
lapack_int LAPACKE_zspsv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int nrhs, lapack_complex_double* ap,
                               lapack_int* ipiv, lapack_complex_double* b,
                               lapack_int ldb );

lapack_int LAPACKE_dspsvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int nrhs, const double* ap,
                                double* afp, lapack_int* ipiv, const double* b,
                                lapack_int ldb, double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, lapack_int* iwork );
lapack_int LAPACKE_zspsvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int nrhs,
                                const lapack_complex_double* ap,
                                lapack_complex_double* afp, lapack_int* ipiv,
                                const lapack_complex_double* b, lapack_int ldb,
                                lapack_complex_double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_dsptrd_work( int matrix_layout, char uplo, lapack_int n,
                                double* ap, double* d, double* e, double* tau );

lapack_int LAPACKE_dsptrf_work( int matrix_layout, char uplo, lapack_int n,
                                double* ap, lapack_int* ipiv );
lapack_int LAPACKE_zsptrf_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* ap, lapack_int* ipiv );

lapack_int LAPACKE_dsptri_work( int matrix_layout, char uplo, lapack_int n,
                                double* ap, const lapack_int* ipiv,
                                double* work );
lapack_int LAPACKE_zsptri_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* ap,
                                const lapack_int* ipiv,
                                lapack_complex_double* work );

lapack_int LAPACKE_dsptrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs, const double* ap,
                                const lapack_int* ipiv, double* b,
                                lapack_int ldb );
lapack_int LAPACKE_zsptrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs,
                                const lapack_complex_double* ap,
                                const lapack_int* ipiv,
                                lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dstebz_work( char range, char order, lapack_int n, double vl,
                                double vu, lapack_int il, lapack_int iu,
                                double abstol, const double* d, const double* e,
                                lapack_int* m, lapack_int* nsplit, double* w,
                                lapack_int* iblock, lapack_int* isplit,
                                double* work, lapack_int* iwork );

lapack_int LAPACKE_dstedc_work( int matrix_layout, char compz, lapack_int n,
                                double* d, double* e, double* z, lapack_int ldz,
                                double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );
lapack_int LAPACKE_zstedc_work( int matrix_layout, char compz, lapack_int n,
                                double* d, double* e, lapack_complex_double* z,
                                lapack_int ldz, lapack_complex_double* work,
                                lapack_int lwork, double* rwork,
                                lapack_int lrwork, lapack_int* iwork,
                                lapack_int liwork );

lapack_int LAPACKE_dstein_work( int matrix_layout, lapack_int n, const double* d,
                                const double* e, lapack_int m, const double* w,
                                const lapack_int* iblock,
                                const lapack_int* isplit, double* z,
                                lapack_int ldz, double* work, lapack_int* iwork,
                                lapack_int* ifailv );
lapack_int LAPACKE_zstein_work( int matrix_layout, lapack_int n, const double* d,
                                const double* e, lapack_int m, const double* w,
                                const lapack_int* iblock,
                                const lapack_int* isplit,
                                lapack_complex_double* z, lapack_int ldz,
                                double* work, lapack_int* iwork,
                                lapack_int* ifailv );

lapack_int LAPACKE_dstemr_work( int matrix_layout, char jobz, char range,
                                lapack_int n, double* d, double* e, double vl,
                                double vu, lapack_int il, lapack_int iu,
                                lapack_int* m, double* w, double* z,
                                lapack_int ldz, lapack_int nzc,
                                lapack_int* isuppz, lapack_logical* tryrac,
                                double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );
lapack_int LAPACKE_zstemr_work( int matrix_layout, char jobz, char range,
                                lapack_int n, double* d, double* e, double vl,
                                double vu, lapack_int il, lapack_int iu,
                                lapack_int* m, double* w,
                                lapack_complex_double* z, lapack_int ldz,
                                lapack_int nzc, lapack_int* isuppz,
                                lapack_logical* tryrac, double* work,
                                lapack_int lwork, lapack_int* iwork,
                                lapack_int liwork );

lapack_int LAPACKE_dsteqr_work( int matrix_layout, char compz, lapack_int n,
                                double* d, double* e, double* z, lapack_int ldz,
                                double* work );
lapack_int LAPACKE_zsteqr_work( int matrix_layout, char compz, lapack_int n,
                                double* d, double* e, lapack_complex_double* z,
                                lapack_int ldz, double* work );

lapack_int LAPACKE_dsterf_work( lapack_int n, double* d, double* e );

lapack_int LAPACKE_dstev_work( int matrix_layout, char jobz, lapack_int n,
                               double* d, double* e, double* z, lapack_int ldz,
                               double* work );

lapack_int LAPACKE_dstevd_work( int matrix_layout, char jobz, lapack_int n,
                                double* d, double* e, double* z, lapack_int ldz,
                                double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_dstevr_work( int matrix_layout, char jobz, char range,
                                lapack_int n, double* d, double* e, double vl,
                                double vu, lapack_int il, lapack_int iu,
                                double abstol, lapack_int* m, double* w,
                                double* z, lapack_int ldz, lapack_int* isuppz,
                                double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_dstevx_work( int matrix_layout, char jobz, char range,
                                lapack_int n, double* d, double* e, double vl,
                                double vu, lapack_int il, lapack_int iu,
                                double abstol, lapack_int* m, double* w,
                                double* z, lapack_int ldz, double* work,
                                lapack_int* iwork, lapack_int* ifail );

lapack_int LAPACKE_dsycon_work( int matrix_layout, char uplo, lapack_int n,
                                const double* a, lapack_int lda,
                                const lapack_int* ipiv, double anorm,
                                double* rcond, double* work,
                                lapack_int* iwork );
lapack_int LAPACKE_zsycon_work( int matrix_layout, char uplo, lapack_int n,
                                const lapack_complex_double* a, lapack_int lda,
                                const lapack_int* ipiv, double anorm,
                                double* rcond, lapack_complex_double* work );

lapack_int LAPACKE_dsyev_work( int matrix_layout, char jobz, char uplo,
                               lapack_int n, double* a, lapack_int lda,
                               double* w, double* work, lapack_int lwork );

lapack_int LAPACKE_dsyevd_work( int matrix_layout, char jobz, char uplo,
                                lapack_int n, double* a, lapack_int lda,
                                double* w, double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_dsyevr_work( int matrix_layout, char jobz, char range,
                                char uplo, lapack_int n, double* a,
                                lapack_int lda, double vl, double vu,
                                lapack_int il, lapack_int iu, double abstol,
                                lapack_int* m, double* w, double* z,
                                lapack_int ldz, lapack_int* isuppz,
                                double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_dsyevx_work( int matrix_layout, char jobz, char range,
                                char uplo, lapack_int n, double* a,
                                lapack_int lda, double vl, double vu,
                                lapack_int il, lapack_int iu, double abstol,
                                lapack_int* m, double* w, double* z,
                                lapack_int ldz, double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int* ifail );

lapack_int LAPACKE_dsygv_work( int matrix_layout, lapack_int itype, char jobz,
                               char uplo, lapack_int n, double* a,
                               lapack_int lda, double* b, lapack_int ldb,
                               double* w, double* work, lapack_int lwork );

lapack_int LAPACKE_dsygvd_work( int matrix_layout, lapack_int itype, char jobz,
                                char uplo, lapack_int n, double* a,
                                lapack_int lda, double* b, lapack_int ldb,
                                double* w, double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );

lapack_int LAPACKE_dsygvx_work( int matrix_layout, lapack_int itype, char jobz,
                                char range, char uplo, lapack_int n, double* a,
                                lapack_int lda, double* b, lapack_int ldb,
                                double vl, double vu, lapack_int il,
                                lapack_int iu, double abstol, lapack_int* m,
                                double* w, double* z, lapack_int ldz,
                                double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int* ifail );

lapack_int LAPACKE_zsyr_work( int matrix_layout, char uplo, lapack_int n,
                                  lapack_complex_double alpha,
                                  const lapack_complex_double* x,
                                  lapack_int incx, lapack_complex_double* a,
                                  lapack_int lda );

lapack_int LAPACKE_dsysv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int nrhs, double* a, lapack_int lda,
                               lapack_int* ipiv, double* b, lapack_int ldb,
                               double* work, lapack_int lwork );
lapack_int LAPACKE_zsysv_work( int matrix_layout, char uplo, lapack_int n,
                               lapack_int nrhs, lapack_complex_double* a,
                               lapack_int lda, lapack_int* ipiv,
                               lapack_complex_double* b, lapack_int ldb,
                               lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_dsysvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int nrhs, const double* a,
                                lapack_int lda, double* af, lapack_int ldaf,
                                lapack_int* ipiv, const double* b,
                                lapack_int ldb, double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                double* work, lapack_int lwork,
                                lapack_int* iwork );
lapack_int LAPACKE_zsysvx_work( int matrix_layout, char fact, char uplo,
                                lapack_int n, lapack_int nrhs,
                                const lapack_complex_double* a, lapack_int lda,
                                lapack_complex_double* af, lapack_int ldaf,
                                lapack_int* ipiv,
                                const lapack_complex_double* b, lapack_int ldb,
                                lapack_complex_double* x, lapack_int ldx,
                                double* rcond, double* ferr, double* berr,
                                lapack_complex_double* work, lapack_int lwork,
                                double* rwork );

lapack_int LAPACKE_dsytrd_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda, double* d, double* e,
                                double* tau, double* work, lapack_int lwork );

lapack_int LAPACKE_dsytrf_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda, lapack_int* ipiv,
                                double* work, lapack_int lwork );
lapack_int LAPACKE_zsytrf_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                lapack_int* ipiv, lapack_complex_double* work,
                                lapack_int lwork );

lapack_int LAPACKE_dsytri_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda,
                                const lapack_int* ipiv, double* work );
lapack_int LAPACKE_zsytri_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                const lapack_int* ipiv,
                                lapack_complex_double* work );

lapack_int LAPACKE_dsytrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs, const double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                double* b, lapack_int ldb );
lapack_int LAPACKE_zsytrs_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_int nrhs, const lapack_complex_double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dtbcon_work( int matrix_layout, char norm, char uplo,
                                char diag, lapack_int n, lapack_int kd,
                                const double* ab, lapack_int ldab,
                                double* rcond, double* work,
                                lapack_int* iwork );
lapack_int LAPACKE_ztbcon_work( int matrix_layout, char norm, char uplo,
                                char diag, lapack_int n, lapack_int kd,
                                const lapack_complex_double* ab,
                                lapack_int ldab, double* rcond,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_dtbtrs_work( int matrix_layout, char uplo, char trans,
                                char diag, lapack_int n, lapack_int kd,
                                lapack_int nrhs, const double* ab,
                                lapack_int ldab, double* b, lapack_int ldb );
lapack_int LAPACKE_ztbtrs_work( int matrix_layout, char uplo, char trans,
                                char diag, lapack_int n, lapack_int kd,
                                lapack_int nrhs,
                                const lapack_complex_double* ab,
                                lapack_int ldab, lapack_complex_double* b,
                                lapack_int ldb );

lapack_int LAPACKE_dtpcon_work( int matrix_layout, char norm, char uplo,
                                char diag, lapack_int n, const double* ap,
                                double* rcond, double* work,
                                lapack_int* iwork );
lapack_int LAPACKE_ztpcon_work( int matrix_layout, char norm, char uplo,
                                char diag, lapack_int n,
                                const lapack_complex_double* ap, double* rcond,
                                lapack_complex_double* work, double* rwork );

lapack_int LAPACKE_dtptri_work( int matrix_layout, char uplo, char diag,
                                lapack_int n, double* ap );
lapack_int LAPACKE_ztptri_work( int matrix_layout, char uplo, char diag,
                                lapack_int n, lapack_complex_double* ap );

lapack_int LAPACKE_dtptrs_work( int matrix_layout, char uplo, char trans,
                                char diag, lapack_int n, lapack_int nrhs,
                                const double* ap, double* b, lapack_int ldb );
lapack_int LAPACKE_ztptrs_work( int matrix_layout, char uplo, char trans,
                                char diag, lapack_int n, lapack_int nrhs,
                                const lapack_complex_double* ap,
                                lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_dtrcon_work( int matrix_layout, char norm, char uplo,
                                char diag, lapack_int n, const double* a,
                                lapack_int lda, double* rcond, double* work,
                                lapack_int* iwork );
lapack_int LAPACKE_ztrcon_work( int matrix_layout, char norm, char uplo,
                                char diag, lapack_int n,
                                const lapack_complex_double* a, lapack_int lda,
                                double* rcond, lapack_complex_double* work,
                                double* rwork );

lapack_int LAPACKE_dtrevc3_work(int matrix_layout, char side, char howmny,
                                lapack_logical* select, lapack_int n, const double* t,
                                lapack_int ldt, double* vl, lapack_int ldvl, double* vr,
                                lapack_int ldvr, lapack_int mm, lapack_int* m,
                                double* work, lapack_int lwork);
lapack_int LAPACKE_ztrevc3_work(int matrix_layout, char side, char howmny,
                                const lapack_logical* select, lapack_int n, lapack_complex_double* t,
                                lapack_int ldt, lapack_complex_double* vl, lapack_int ldvl,
                                lapack_complex_double* vr, lapack_int ldvr, lapack_int mm, lapack_int* m,
                                lapack_complex_double* work, lapack_int lwork, double* rwork, lapack_int lrwork);

lapack_int LAPACKE_dtrexc_work( int matrix_layout, char compq, lapack_int n,
                                double* t, lapack_int ldt, double* q,
                                lapack_int ldq, lapack_int* ifst,
                                lapack_int* ilst, double* work );
lapack_int LAPACKE_ztrexc_work( int matrix_layout, char compq, lapack_int n,
                                lapack_complex_double* t, lapack_int ldt,
                                lapack_complex_double* q, lapack_int ldq,
                                lapack_int ifst, lapack_int ilst );

lapack_int LAPACKE_dtrsen_work( int matrix_layout, char job, char compq,
                                const lapack_logical* select, lapack_int n,
                                double* t, lapack_int ldt, double* q,
                                lapack_int ldq, double* wr, double* wi,
                                lapack_int* m, double* s, double* sep,
                                double* work, lapack_int lwork,
                                lapack_int* iwork, lapack_int liwork );
lapack_int LAPACKE_ztrsen_work( int matrix_layout, char job, char compq,
                                const lapack_logical* select, lapack_int n,
                                lapack_complex_double* t, lapack_int ldt,
                                lapack_complex_double* q, lapack_int ldq,
                                lapack_complex_double* w, lapack_int* m,
                                double* s, double* sep,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_dtrsna_work( int matrix_layout, char job, char howmny,
                                const lapack_logical* select, lapack_int n,
                                const double* t, lapack_int ldt,
                                const double* vl, lapack_int ldvl,
                                const double* vr, lapack_int ldvr, double* s,
                                double* sep, lapack_int mm, lapack_int* m,
                                double* work, lapack_int ldwork,
                                lapack_int* iwork );
lapack_int LAPACKE_ztrsna_work( int matrix_layout, char job, char howmny,
                                const lapack_logical* select, lapack_int n,
                                const lapack_complex_double* t, lapack_int ldt,
                                const lapack_complex_double* vl,
                                lapack_int ldvl,
                                const lapack_complex_double* vr,
                                lapack_int ldvr, double* s, double* sep,
                                lapack_int mm, lapack_int* m,
                                lapack_complex_double* work, lapack_int ldwork,
                                double* rwork );

lapack_int LAPACKE_dtrsyl_work( int matrix_layout, char trana, char tranb,
                                lapack_int isgn, lapack_int m, lapack_int n,
                                const double* a, lapack_int lda,
                                const double* b, lapack_int ldb, double* c,
                                lapack_int ldc, double* scale );
lapack_int LAPACKE_ztrsyl_work( int matrix_layout, char trana, char tranb,
                                lapack_int isgn, lapack_int m, lapack_int n,
                                const lapack_complex_double* a, lapack_int lda,
                                const lapack_complex_double* b, lapack_int ldb,
                                lapack_complex_double* c, lapack_int ldc,
                                double* scale );

lapack_int LAPACKE_dtrtri_work( int matrix_layout, char uplo, char diag,
                                lapack_int n, double* a, lapack_int lda );
lapack_int LAPACKE_ztrtri_work( int matrix_layout, char uplo, char diag,
                                lapack_int n, lapack_complex_double* a,
                                lapack_int lda );

lapack_int LAPACKE_dtrtrs_work( int matrix_layout, char uplo, char trans,
                                char diag, lapack_int n, lapack_int nrhs,
                                const double* a, lapack_int lda, double* b,
                                lapack_int ldb );
lapack_int LAPACKE_ztrtrs_work( int matrix_layout, char uplo, char trans,
                                char diag, lapack_int n, lapack_int nrhs,
                                const lapack_complex_double* a, lapack_int lda,
                                lapack_complex_double* b, lapack_int ldb );

lapack_int LAPACKE_zunghr_work( int matrix_layout, lapack_int n, lapack_int ilo,
                                lapack_int ihi, lapack_complex_double* a,
                                lapack_int lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_zunglq_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int k, lapack_complex_double* a,
                                lapack_int lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_zungqr_work( int matrix_layout, lapack_int m, lapack_int n,
                                lapack_int k, lapack_complex_double* a,
                                lapack_int lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_zungtr_work( int matrix_layout, char uplo, lapack_int n,
                                lapack_complex_double* a, lapack_int lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_zunmhr_work( int matrix_layout, char side, char trans,
                                lapack_int m, lapack_int n, lapack_int ilo,
                                lapack_int ihi, const lapack_complex_double* a,
                                lapack_int lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, lapack_int ldc,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_zunmlq_work( int matrix_layout, char side, char trans,
                                lapack_int m, lapack_int n, lapack_int k,
                                const lapack_complex_double* a, lapack_int lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, lapack_int ldc,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_zunmqr_work( int matrix_layout, char side, char trans,
                                lapack_int m, lapack_int n, lapack_int k,
                                const lapack_complex_double* a, lapack_int lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, lapack_int ldc,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_zunmtr_work( int matrix_layout, char side, char uplo,
                                char trans, lapack_int m, lapack_int n,
                                const lapack_complex_double* a, lapack_int lda,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, lapack_int ldc,
                                lapack_complex_double* work, lapack_int lwork );

lapack_int LAPACKE_zupgtr_work( int matrix_layout, char uplo, lapack_int n,
                                const lapack_complex_double* ap,
                                const lapack_complex_double* tau,
                                lapack_complex_double* q, lapack_int ldq,
                                lapack_complex_double* work );

lapack_int LAPACKE_zupmtr_work( int matrix_layout, char side, char uplo,
                                char trans, lapack_int m, lapack_int n,
                                const lapack_complex_double* ap,
                                const lapack_complex_double* tau,
                                lapack_complex_double* c, lapack_int ldc,
                                lapack_complex_double* work );

/* APIs for set/get nancheck flags */
void LAPACKE_set_nancheck( int flag );
int LAPACKE_get_nancheck( void );

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* _LAPACKE_H_ */
