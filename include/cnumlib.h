/****************************************
 *                                      *
 *  C/C++ Numerical Library             *
 *  Version 6.0 (October 10, 2021)      *
 *  (C) 2014-2021 K Technologies        *
 *                                      *
 ****************************************/
#pragma once

#undef _VLARRAY

#if defined(__cplusplus)

#include <cstddef>
#include <complex>
#define doublecomplex	std::complex<double>
#define floatcomplex	std::complex<float>

#else

#include <stddef.h>
#include <complex.h>
#if defined(__clang__) || defined(__GNUC__) || defined(__APPLE__)
#define doublecomplex	double _Complex
#define floatcomplex	float _Complex
#if !defined(_NO_VLARRAY)
#define _VLARRAY
#endif
#else	// Visual C/C++
#define doublecomplex	_Dcomplex
#define floatcomplex	_Fcomplex
#endif

#endif

#if defined(__cplusplus)
extern "C" {
#endif

/*
 * A3. Real arithmetic
 */
extern double d1num(int i);

/*
 * A4. C utility routines for complex arithmetic
 */
#if !defined(__cplusplus)
extern doublecomplex cmplx(double r, double i);
extern doublecomplex cpolar(double rho, double theta);
extern doublecomplex cminus(doublecomplex a);
extern doublecomplex cadd(doublecomplex a, doublecomplex b);
extern doublecomplex cdadd(doublecomplex a, double rb);
extern doublecomplex csub(doublecomplex a, doublecomplex b);
extern doublecomplex cdsub(doublecomplex a, double rb);
extern doublecomplex dcsub(double ra, doublecomplex b);
extern doublecomplex cmul(doublecomplex a, doublecomplex b);
extern doublecomplex cdmul(doublecomplex a, double rb);
extern doublecomplex cdiv(doublecomplex a, doublecomplex b);
extern doublecomplex cddiv(doublecomplex a, double rb);
extern doublecomplex dcdiv(double ra, doublecomplex b);
extern doublecomplex cdpow(doublecomplex a, double rb);
extern doublecomplex cipow(doublecomplex a, int ib);
#endif
#if !defined(_NO_SUB_DEFS)
extern void conj_sub(doublecomplex z, doublecomplex *conjz);
extern void cpolar_sub(double rho, double theta, doublecomplex *z);
extern void cproj_sub(doublecomplex z, doublecomplex *cprojz);
extern void cmplx_sub(double r, double i, doublecomplex *z);
extern void cminus_sub(doublecomplex a, doublecomplex *minus);
extern void cadd_sub(doublecomplex a, doublecomplex b, doublecomplex *add);
extern void cdadd_sub(doublecomplex a, double rb, doublecomplex *add);
extern void csub_sub(doublecomplex a, doublecomplex b, doublecomplex *sub);
extern void cdsub_sub(doublecomplex a, double rb, doublecomplex *sub);
extern void dcsub_sub(double ra, doublecomplex b, doublecomplex *sub);
extern void cmul_sub(doublecomplex a, doublecomplex b, doublecomplex *mul);
extern void cdmul_sub(doublecomplex a, double rb, doublecomplex *mul);
extern void cdiv_sub(doublecomplex a, doublecomplex b, doublecomplex *div);
extern void cddiv_sub(doublecomplex a, double rb, doublecomplex *div);
extern void dcdiv_sub(double ra, doublecomplex b, doublecomplex *div);
extern void cpow_sub(doublecomplex a, doublecomplex b, doublecomplex *pow);
extern void cdpow_sub(doublecomplex a, double rb, doublecomplex *pow);
extern void cipow_sub(doublecomplex a, int ib, doublecomplex *pow);
#endif

/*
 * C1. Special functions (Integer-valued functions)
 */
extern double factorial(unsigned int n);

/*
 * C2. Special functions (Powers, roots, reciprocals)
 */
#if !defined(__cplusplus)
extern doublecomplex ccbrt(doublecomplex z);
#endif
#if !defined(_NO_SUB_DEFS)
extern void ccbrt_sub(doublecomplex z, doublecomplex *cbrtz);
extern void csqrt_sub(doublecomplex z, doublecomplex *sqrtz);
#endif

/*
 * C3. Special functions (Polynomials)
 */
extern double laguerre(unsigned int n, double x);
extern double alaguerre(unsigned int n, unsigned int m, double x);
extern double legendre(unsigned int n, double x);
extern double legendred(unsigned int n, double x);
extern double alegendre(unsigned int n, unsigned int m, double x);
#if !defined(__cplusplus)
extern doublecomplex sharmonic(unsigned int l, int m, double theta, double phi);
#endif
#if !defined(_NO_SUB_DEFS)
extern void sharmonic_sub(unsigned int l, int m, double theta, double phi, doublecomplex *z);
#endif
extern double sharmonicr(unsigned int l, int m, double theta, double phi);
extern double sharmonici(unsigned int l, int m, double theta, double phi);
extern double hermite(unsigned int n, double x);
extern double chebt(unsigned int n, double x);
extern double chebtd(unsigned int n, double x);
extern double chebu(unsigned int n, double x);
extern double chebs(double c[], size_t n, double x);

/*
 * C4. Special functions (Elementary transcendental functions)
 */
#if !defined(__cplusplus)
extern doublecomplex cexpm1(doublecomplex z);
extern doublecomplex clog1p(doublecomplex z);
extern doublecomplex ccot(doublecomplex z);
#endif
#if !defined(_NO_SUB_DEFS)
extern void cexpm1_sub(doublecomplex z, doublecomplex *expm1z);
extern void clog1p_sub(doublecomplex z, doublecomplex *log1pz);
extern void cexp_sub(doublecomplex z, doublecomplex *expz);
extern void clog_sub(doublecomplex z, doublecomplex *logz);
extern void ccos_sub(doublecomplex z, doublecomplex *cosz);
extern void csin_sub(doublecomplex z, doublecomplex *sinz);
extern void ctan_sub(doublecomplex z, doublecomplex *tanz);
extern void ccosh_sub(doublecomplex z, doublecomplex *coshz);
extern void csinh_sub(doublecomplex z, doublecomplex *sinhz);
extern void ctanh_sub(doublecomplex z, doublecomplex *tanhz);
extern void cacos_sub(doublecomplex z, doublecomplex *acosz);
extern void casin_sub(doublecomplex z, doublecomplex *asinz);
extern void catan_sub(doublecomplex z, doublecomplex *atanz);
extern void cacosh_sub(doublecomplex z, doublecomplex *acoshz);
extern void casinh_sub(doublecomplex z, doublecomplex *asinhz);
extern void catanh_sub(doublecomplex z, doublecomplex *atanhz);
extern void ccot_sub(doublecomplex z, doublecomplex *cotz);
#endif
extern double sqrt1pm1(double x);
extern double powm1(double x, double y);
extern double sin_pi(double x);
extern double cos_pi(double x);

/*
 * C5. Special functions (Exponential and logarithmic integrals)
 */
extern double li(double x);
extern double ei(double x);
extern double e1(double x);
extern double en(unsigned int n, double x);
extern double spence(double x);

/*
 * C6. Special functions (Cosine and sine integrals)
 */
extern double ci(double x);
extern double si(double x);
extern double chi(double x);
extern double shi(double x);

/*
 * C7a. Special functions (Gamma functions)
 */
extern double tgamma1pm1(double x);
extern double lgammas(double x, int *sign);
extern double rgamma(double x);
extern double tgammaratio(double a, double b);
extern double tgammadratio(double a, double delta);
#if !defined(__cplusplus)
extern doublecomplex cgamma(doublecomplex z);
extern doublecomplex clgamma(doublecomplex z);
extern doublecomplex crgamma(doublecomplex z);
#endif
#if !defined(_NO_SUB_DEFS)
extern void cgamma_sub(doublecomplex z, doublecomplex *zout);
extern void clgamma_sub(doublecomplex z, doublecomplex *zout);
extern void crgamma_sub(doublecomplex z, doublecomplex *zout);
#endif
extern double poch(double a, double x);
extern double poch1(double a, double x);

/*
 * C7b. Special functions (Beta functions)
 */
extern double beta(double a, double b);
extern double lbeta(double a, double b);
#if !defined(__cplusplus)
extern doublecomplex cbeta(doublecomplex a, doublecomplex b);
extern doublecomplex clbeta(doublecomplex a, doublecomplex b);
#endif
#if !defined(_NO_SUB_DEFS)
extern void cbeta_sub(doublecomplex a, doublecomplex b, doublecomplex *zout);
extern void clbeta_sub(doublecomplex a, doublecomplex b, doublecomplex *zout);
#endif

/*
 * C7c. Special functions (Psi function)
 */
extern double digamma(double x);
extern double trigamma(double x);
extern double polygamma(int n, double x);
#if !defined(__cplusplus)
extern doublecomplex cdigamma(doublecomplex z);
#endif
#if !defined(_NO_SUB_DEFS)
extern void cdigamma_sub(doublecomplex z, doublecomplex *zout);
#endif

/*
 * C7e. Special functions (Incomplete Gamma functions)
 */
extern double gammai(double a, double x);
extern double gammaic(double a, double x);
extern double gammait(double a, double x);
extern double gammap(double a, double x);
extern double gammaq(double a, double x);
extern double gammapi(double a, double p);
extern double gammaqi(double a, double q);
extern double gammapia(double x, double p);
extern double gammaqia(double x, double q);
extern double gammapd(double a, double x);

/*
 * C7f. Special functions (Incomplete Beta function)
 */
extern double betax(double a, double b, double x);
extern double betaxc(double a, double b, double x);
extern double ibeta(double a, double b, double x);
extern double ibetac(double a, double b, double x);
extern double ibetai(double a, double b, double p, double *py);
extern double ibetaci(double a, double b, double q, double *py);
extern double ibetaia(double b, double x, double p);
extern double ibetacia(double b, double x, double q);
extern double ibetaib(double a, double x, double p);
extern double ibetacib(double a, double x, double q);
extern double ibetad(double a, double b, double x);

/*
 * C7g. Special functions (Riemann zeta function)
 */
extern double zeta(double x);

/*
 * C8. Special functions (Error functions)
 */
extern double erfi(double p);
extern double erfci(double q);
extern double dawson(double x);
extern double fresc(double x);
extern double fress(double x);

/*
 * C10a. Special functions (Bessel functions)
 */
extern double besj0(double x);
extern double besj1(double x);
extern double besjn(int n, double x);
extern double besy0(double x);
extern double besy1(double x);
extern double besyn(int n, double x);
extern double besjnu(double nu, double x);
extern double besynu(double nu, double x);
extern double besjnd(int n, double x);
extern double besynd(int n, double x);
extern double besjnud(double nu, double x);
extern double besynud(double nu, double x);
extern double sbesjn(int n, double x);
extern double sbesyn(int n, double x);
extern double sbesjnu(double nu, double x);
extern double sbesynu(double nu, double x);
extern void cbesh(doublecomplex z, double nu, int kode, int m, int n, doublecomplex y[], int *info);
extern void cbesj(doublecomplex z, double nu, int kode, int n, doublecomplex y[], int *info);
extern void cbesy(doublecomplex z, double nu, int kode, int n, doublecomplex y[], doublecomplex work[], int *info);

/*
 * C10b. Special functions (Modified Bessel functions)
 */
extern double besi0(double x);
extern double besi0e(double x);
extern double besi1(double x);
extern double besi1e(double x);
extern double besin(int n, double x);
extern double besk0(double x);
extern double besk0e(double x);
extern double besk1(double x);
extern double besk1e(double x);
extern double beskn(int n, double x);
extern double besinu(double nu, double x);
extern double besknu(double nu, double x);
extern double besind(int n, double x);
extern double besknd(int n, double x);
extern double besinud(double nu, double x);
extern double besknud(double nu, double x);
extern double sbesin(int n, double x);
extern double sbeskn(int n, double x);
extern double sbesinu(double nu, double x);
extern double sbesknu(double nu, double x);
extern void cbesi(doublecomplex z, double nu, int kode, int n, doublecomplex y[], int *info);
extern void cbesk(doublecomplex z, double nu, int kode, int n, doublecomplex y[], int *info);

/*
 * C10d. Special functions (Airy functions)
 */
extern double airyai(double x);
extern double airybi(double x);
extern double airyaid(double x);
extern double airybid(double x);
extern void cairy(doublecomplex z, int id, int kode, doublecomplex *ai, int *info);
extern void cbiry(doublecomplex z, int id, int kode, doublecomplex *bi, int *info);

/*
 * C11. Special functions (Confluent hypergeometric function)
 */
extern double chu(double a, double b, double x);
extern double hyp1f1(double a, double b, double x);
extern double hyp1f1r(double a, double b, double x);
extern double lhyp1f1(double a, double b, double z, int *sign);
extern double hyp2f1(double a, double b, double c, double x);
extern double hyp0f1(double b, double z);
extern double hyp1f0(double a, double z);
extern double hyp2f0(double a1, double a2, double z);

/*
 * C13. Special functions (Jacobi elliptic functions)
 */
extern void jelli(double u, double k, double *sn, double *cn, double *dn);
extern double jsn(double u, double k);
extern double jcn(double u, double k);
extern double jdn(double u, double k);
extern double jns(double u, double k);
extern double jnc(double u, double k);
extern double jnd(double u, double k);
extern double jsc(double u, double k);
extern double jsd(double u, double k);
extern double jdc(double u, double k);
extern double jds(double u, double k);
extern double jcs(double u, double k);
extern double jcd(double u, double k);

/*
 * C14. Special functions (Elliptic integrals)
 */
extern double celli1(double k);
extern double celli2(double k);
extern double celli3(double n, double k);
extern double elli1(double phi, double k);
extern double elli2(double phi, double k);
extern double elli3(double phi, double n, double k);
extern double rc(double x, double y);
extern double rd(double x, double y, double z);
extern double rg(double x, double y, double z);
extern double rf(double x, double y, double z);
extern double rj(double x, double y, double z, double p);
extern double jzeta(double phi, double k);

/*
 * C19. Special functions (other special functions)
 */
extern double dconst(int i);

/*
 * D1a. Elementary vector operations (BLAS 1)
 */
extern void daxpy(int n, double a, double x[], int incx, double y[], int incy);
extern void dcopy(int n, double x[], int incx, double y[], int incy);
extern double ddot(int n, double x[], int incx, double y[], int incy);
extern void drotg(double *a, double *b, double *c, double *s);
extern void drotmg(double *d1, double *d2, double *x1, double y1, double p[]);
extern void drot(int n, double x[], int incx, double y[], int incy, double c, double s);
extern void drotm(int n, double x[], int incx, double y[], int incy, double p[]);
extern void dscal(int n, double a, double x[], int incx);
extern void dswap(int n, double x[], int incx, double y[], int incy);
extern double dasum(int n, double x[], int incx);
extern double dnrm2(int n, double x[], int incx);
extern int idamax(int n, double x[], int incx);
extern int lsame(char ca, char cb);

/*
 * D1a-2. Elementary vector operations (BLAS 1) (Complex)
 */
extern void zaxpy(int n, doublecomplex a, doublecomplex x[], int incx, doublecomplex y[], int incy);
extern void zcopy(int n, doublecomplex x[], int incx, doublecomplex y[], int incy);
#if !defined(__cplusplus)
extern doublecomplex zdotu(int n, doublecomplex x[], int incx, doublecomplex y[], int incy);
extern doublecomplex zdotc(int n, doublecomplex x[], int incx, doublecomplex y[], int incy);
#endif
#if !defined(_NO_SUB_DEFS)
extern void zdotu_sub(doublecomplex *dotu, int n, doublecomplex zx[], int incx, doublecomplex zy[], int incy);
extern void zdotc_sub(doublecomplex *dotc, int n, doublecomplex zx[], int incx, doublecomplex zy[], int incy);
#endif
extern void zrotg(doublecomplex *a, doublecomplex *b, double *c, doublecomplex *s);
extern void zrot(int n, doublecomplex x[], int incx, doublecomplex y[], int incy, double c, doublecomplex s);
extern void zdrot(int n, doublecomplex x[], int incx, doublecomplex y[], int incy, double c, double s);
extern void zdscal(int n, double a, doublecomplex x[], int incx);
extern void zscal(int n, doublecomplex a, doublecomplex x[], int incx);
extern void zswap(int n, doublecomplex x[], int incx, doublecomplex y[], int incy);
extern double dzasum(int n, doublecomplex x[], int incx);
extern double dznrm2(int n, doublecomplex x[], int incx);
extern double dcabs1(doublecomplex z);
extern int izamax(int n, doublecomplex z[], int incz);

/*
 * D1a-3. Elementary vector operations (BLAS 2)
 */
#if defined(_VLARRAY)
extern void dgemv(char trans, int m, int n, double alpha, int lda, double a[][lda], double x[], int incx, double beta, double y[], int incy);
extern void dgbmv(char trans, int m, int n, int kl, int ku, double alpha, int ldab, double ab[][ldab], double x[], int incx, double beta, double y[], int incy);
extern void dsymv(char uplo, int n, double alpha, int lda, double a[][lda], double x[], int incx, double beta, double y[], int incy);
extern void dsbmv(char uplo, int n, int k, double alpha, int ldab, double ab[][ldab], double x[], int incx, double beta, double y[], int incy);
extern void dspmv(char uplo, int n, double alpha, double ap[], double x[], int incx, double beta, double y[], int incy);
extern void dtrmv(char uplo, char trans, char diag, int n, int lda, double a[][lda], double x[], int incx);
extern void dtbmv(char uplo, char trans, char diag, int n, int k, int ldab, double ab[][ldab], double x[], int incx);
extern void dtpmv(char uplo, char trans, char diag, int n, double ap[], double x[], int incx);
extern void dtrsv(char uplo, char trans, char diag, int n, int lda, double a[][lda], double x[], int incx);
extern void dtbsv(char uplo, char trans, char diag, int n, int k, int ldab, double ab[][ldab], double x[], int incx);
extern void dtpsv(char uplo, char trans, char diag, int n, double ap[], double x[], int incx);
extern void dger(int m, int n, double alpha, double x[], int incx, double y[], int incy, int lda, double a[][lda]);
extern void dsyr(char uplo, int n, double alpha, double x[], int incx, int lda, double a[][lda]);
extern void dspr(char uplo, int n, double alpha, double x[], int incx, double ap[]);
extern void dsyr2(char uplo, int n, double alpha, double x[], int incx, double y[], int incy, int lda, double a[][lda]);
extern void dspr2(char uplo, int n, double alpha, double x[], int incx, double y[], int incy, double ap[]);
#else
extern void dgemv(char trans, int m, int n, double alpha, int lda, double a[], double x[], int incx, double beta, double y[], int incy);
extern void dgbmv(char trans, int m, int n, int kl, int ku, double alpha, int ldab, double ab[], double x[], int incx, double beta, double y[], int incy);
extern void dsymv(char uplo, int n, double alpha, int lda, double a[], double x[], int incx, double beta, double y[], int incy);
extern void dsbmv(char uplo, int n, int k, double alpha, int ldab, double ab[], double x[], int incx, double beta, double y[], int incy);
extern void dspmv(char uplo, int n, double alpha, double ap[], double x[], int incx, double beta, double y[], int incy);
extern void dtrmv(char uplo, char trans, char diag, int n, int lda, double a[], double x[], int incx);
extern void dtbmv(char uplo, char trans, char diag, int n, int k, int ldab, double ab[], double x[], int incx);
extern void dtpmv(char uplo, char trans, char diag, int n, double ap[], double x[], int incx);
extern void dtrsv(char uplo, char trans, char diag, int n, int lda, double a[], double x[], int incx);
extern void dtbsv(char uplo, char trans, char diag, int n, int k, int ldab, double ab[], double x[], int incx);
extern void dtpsv(char uplo, char trans, char diag, int n, double ap[], double x[], int incx);
extern void dger(int m, int n, double alpha, double x[], int incx, double y[], int incy, int lda, double a[]);
extern void dsyr(char uplo, int n, double alpha, double x[], int incx, int lda, double a[]);
extern void dspr(char uplo, int n, double alpha, double x[], int incx, double ap[]);
extern void dsyr2(char uplo, int n, double alpha, double x[], int incx, double y[], int incy, int lda, double a[]);
extern void dspr2(char uplo, int n, double alpha, double x[], int incx, double y[], int incy, double ap[]);
#endif

/*
 * D1a-4. Elementary vector operations (BLAS 2) (Complex)
 */
#if defined(_VLARRAY)
extern void zgemv(char trans, int m, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zgbmv(char trans, int m, int n, int kl, int ku, doublecomplex alpha, int ldab, doublecomplex ab[][ldab], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zhemv(char uplo, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zhbmv(char uplo, int n, int k, doublecomplex alpha, int ldab, doublecomplex ab[][ldab], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zhpmv(char uplo, int n, doublecomplex alpha, doublecomplex ap[], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zsymv(char uplo, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zsbmv(char uplo, int n, int k, doublecomplex alpha, int ldab, doublecomplex ab[][ldab], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zspmv(char uplo, int n, doublecomplex alpha, doublecomplex ap[], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void ztrmv(char uplo, char trans, char diag, int n, int lda, doublecomplex a[][lda], doublecomplex x[], int incx);
extern void ztbmv(char uplo, char trans, char diag, int n, int k, int ldab, doublecomplex ab[][ldab], doublecomplex x[], int incx);
extern void ztpmv(char uplo, char trans, char diag, int n, doublecomplex ap[], doublecomplex x[], int incx);
extern void ztrsv(char uplo, char trans, char diag, int n, int lda, doublecomplex a[][lda], doublecomplex x[], int incx);
extern void ztbsv(char uplo, char trans, char diag, int n, int k, int ldab, doublecomplex ab[][ldab], doublecomplex x[], int incx);
extern void ztpsv(char uplo, char trans, char diag, int n, doublecomplex ap[], doublecomplex x[], int incx);
extern void zgeru(int m, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[][lda]);
extern void zgerc(int m, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[][lda]);
extern void zher(char uplo, int n, double alpha, doublecomplex x[], int incx, int lda, doublecomplex a[][lda]);
extern void zhpr(char uplo, int n, double alpha, doublecomplex x[], int incx, doublecomplex ap[]);
extern void zsyr(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, int lda, doublecomplex a[][lda]);
extern void zspr(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex ap[]);
extern void zher2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[][lda]);
extern void zhpr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, doublecomplex ap[]);
extern void zsyr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[][lda]);
extern void zspr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, doublecomplex ap[]);
#else
extern void zgemv(char trans, int m, int n, doublecomplex alpha, int lda, doublecomplex a[], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zgbmv(char trans, int m, int n, int kl, int ku, doublecomplex alpha, int ldab, doublecomplex ab[], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zhemv(char uplo, int n, doublecomplex alpha, int lda, doublecomplex a[], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zhbmv(char uplo, int n, int k, doublecomplex alpha, int ldab, doublecomplex ab[], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zhpmv(char uplo, int n, doublecomplex alpha, doublecomplex ap[], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zsymv(char uplo, int n, doublecomplex alpha, int lda, doublecomplex a[], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zsbmv(char uplo, int n, int k, doublecomplex alpha, int ldab, doublecomplex ab[], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void zspmv(char uplo, int n, doublecomplex alpha, doublecomplex ap[], doublecomplex x[], int incx, doublecomplex beta, doublecomplex y[], int incy);
extern void ztrmv(char uplo, char trans, char diag, int n, int lda, doublecomplex a[], doublecomplex x[], int incx);
extern void ztbmv(char uplo, char trans, char diag, int n, int k, int ldab, doublecomplex ab[], doublecomplex x[], int incx);
extern void ztpmv(char uplo, char trans, char diag, int n, doublecomplex ap[], doublecomplex x[], int incx);
extern void ztrsv(char uplo, char trans, char diag, int n, int lda, doublecomplex a[], doublecomplex x[], int incx);
extern void ztbsv(char uplo, char trans, char diag, int n, int k, int ldab, doublecomplex ab[], doublecomplex x[], int incx);
extern void ztpsv(char uplo, char trans, char diag, int n, doublecomplex ap[], doublecomplex x[], int incx);
extern void zgeru(int m, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[]);
extern void zgerc(int m, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[]);
extern void zher(char uplo, int n, double alpha, doublecomplex x[], int incx, int lda, doublecomplex a[]);
extern void zhpr(char uplo, int n, double alpha, doublecomplex x[], int incx, doublecomplex ap[]);
extern void zsyr(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, int lda, doublecomplex a[]);
extern void zspr(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex ap[]);
extern void zher2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[]);
extern void zhpr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, doublecomplex ap[]);
extern void zsyr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[]);
extern void zspr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, doublecomplex ap[]);
#endif

/*
 * D1b. Elementary matrix operations (BLAS 3)
 */
#if defined(_VLARRAY)
extern void dgemm(char transa, char transb, int m, int n, int k, double alpha, int lda, double a[][lda], int ldb, double b[][ldb], double beta, int ldc, double c[][ldc]);
extern void dsymm(char side, char uplo, int m, int n, double alpha, int lda, double a[][lda], int ldb, double b[][ldb], double beta, int ldc, double c[][ldc]);
extern void dsyrk(char uplo, char trans, int n, int k, double alpha, int lda, double a[][lda], double beta, int ldc, double c[][ldc]);
extern void dsyr2k(char uplo, char trans, int n, int k, double alpha, int lda, double a[][lda], int ldb, double b[][ldb], double beta, int ldc, double c[][ldc]);
extern void dtrmm(char side, char uplo, char transa, char diag, int m, int n, double alpha, int lda, double a[][lda], int ldb, double b[][ldb]);
extern void dtrsm(char side, char uplo, char transa, char diag, int m, int n, double alpha, int lda, double a[][lda], int ldb, double b[][ldb]);
#else
extern void dgemm(char transa, char transb, int m, int n, int k, double alpha, int lda, double a[], int ldb, double b[], double beta, int ldc, double c[]);
extern void dsymm(char side, char uplo, int m, int n, double alpha, int lda, double a[], int ldb, double b[], double beta, int ldc, double c[]);
extern void dsyrk(char uplo, char trans, int n, int k, double alpha, int lda, double a[], double beta, int ldc, double c[]);
extern void dsyr2k(char uplo, char trans, int n, int k, double alpha, int lda, double a[], int ldb, double b[], double beta, int ldc, double c[]);
extern void dtrmm(char side, char uplo, char transa, char diag, int m, int n, double alpha, int lda, double a[], int ldb, double b[]);
extern void dtrsm(char side, char uplo, char transa, char diag, int m, int n, double alpha, int lda, double a[], int ldb, double b[]);
#endif

/*
 * D1b-2. Elementary matrix operations (BLAS 3) (Complex)
 */
#if defined(_VLARRAY)
extern void zgemm(char transa, char transb, int m, int n, int k, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex beta, int ldc, doublecomplex c[][ldc]);
extern void zsymm(char side, char uplo, int m, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex beta, int ldc, doublecomplex c[][ldc]);
extern void zhemm(char side, char uplo, int m, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex beta, int ldc, doublecomplex c[][ldc]);
extern void zsyrk(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[][lda], doublecomplex beta, int ldc, doublecomplex c[][ldc]);
extern void zherk(char uplo, char trans, int n, int k, double alpha, int lda, doublecomplex a[][lda], double beta, int ldc, doublecomplex c[][ldc]);
extern void zsyr2k(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex beta, int ldc, doublecomplex c[][ldc]);
extern void zher2k(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double beta, int ldc, doublecomplex c[][ldc]);
extern void ztrmm(char side, char uplo, char transa, char diag, int m, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb]);
extern void ztrsm(char side, char uplo, char transa, char diag, int m, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb]);
#else
extern void zgemm(char transa, char transb, int m, int n, int k, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex beta, int ldc, doublecomplex c[]);
extern void zsymm(char side, char uplo, int m, int n, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex beta, int ldc, doublecomplex c[]);
extern void zhemm(char side, char uplo, int m, int n, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex beta, int ldc, doublecomplex c[]);
extern void zsyrk(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[], doublecomplex beta, int ldc, doublecomplex c[]);
extern void zherk(char uplo, char trans, int n, int k, double alpha, int lda, doublecomplex a[], double beta, int ldc, doublecomplex c[]);
extern void zsyr2k(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex beta, int ldc, doublecomplex c[]);
void zher2k(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[], double beta, int ldc, doublecomplex c[]);
extern void ztrmm(char side, char uplo, char transa, char diag, int m, int n, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[]);
extern void ztrsm(char side, char uplo, char transa, char diag, int m, int n, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[]);
#endif

/*
 * D1b-3. Elementary matrix operations (LAPACK auxiliary routines)
 */
#if defined(_VLARRAY)
extern double dlange(char norm, int m, int n, int lda, double a[][lda], double work[]);
extern double dlangb(char norm, int n, int kl, int ku, int ldab, double ab[][ldab], double work[]);
extern double dlangt(char norm, int n, double dl[], double d[], double du[]);
extern double dlansy(char norm, char uplo, int n, int lda, double a[][lda], double work[]);
extern double dlansb(char norm, char uplo, int n, int k, int ldab, double ab[][ldab], double work[]);
extern double dlansp(char norm, char uplo, int n, double ap[], double work[]);
extern double dlanst(char norm, int n, double d[], double e[]);
extern double dlantr(char norm, char uplo, char diag, int m, int n, int lda, double a[][lda], double work[]);
#else
extern double dlange(char norm, int m, int n, int lda, double a[], double work[]);
extern double dlangb(char norm, int n, int kl, int ku, int ldab, double ab[], double work[]);
extern double dlangt(char norm, int n, double dl[], double d[], double du[]);
extern double dlansy(char norm, char uplo, int n, int lda, double a[], double work[]);
extern double dlansb(char norm, char uplo, int n, int k, int ldab, double ab[], double work[]);
extern double dlansp(char norm, char uplo, int n, double ap[], double work[]);
extern double dlanst(char norm, int n, double d[], double e[]);
extern double dlantr(char norm, char uplo, char diag, int m, int n, int lda, double a[], double work[]);
#endif

/*
 * D1b-4. Elementary matrix operations (LAPACK auxiliary routines) (Complex)
 */
#if defined(_VLARRAY)
extern double zlange(char norm, int m, int n, int lda, doublecomplex a[][lda], double work[]);
extern double zlangb(char norm, int n, int kl, int ku, int ldab, doublecomplex ab[][ldab], double work[]);
extern double zlangt(char norm, int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[]);
extern double zlansy(char norm, char uplo, int n, int lda, doublecomplex a[][lda], double work[]);
extern double zlansb(char norm, char uplo, int n, int k, int ldab, doublecomplex ab[][ldab], double work[]);
extern double zlansp(char norm, char uplo, int n, doublecomplex ap[], double work[]);
extern double zlanhe(char norm, char uplo, int n, int lda, doublecomplex a[][lda], double work[]);
extern double zlanhb(char norm, char uplo, int n, int k, int ldab, doublecomplex ab[][ldab], double work[]);
extern double zlanhp(char norm, char uplo, int n, doublecomplex ap[], double work[]);
extern double zlanht(char norm, int n, double d[], doublecomplex e[]);
extern double zlantr(char norm, char uplo, char diag, int m, int n, int lda, doublecomplex a[][lda], double work[]);
#else
extern double zlange(char norm, int m, int n, int lda, doublecomplex a[], double work[]);
extern double zlangb(char norm, int n, int kl, int ku, int ldab, doublecomplex ab[], double work[]);
extern double zlangt(char norm, int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[]);
extern double zlansy(char norm, char uplo, int n, int lda, doublecomplex a[], double work[]);
extern double zlansb(char norm, char uplo, int n, int k, int ldab, doublecomplex ab[], double work[]);
extern double zlansp(char norm, char uplo, int n, doublecomplex ap[], double work[]);
extern double zlanhe(char norm, char uplo, int n, int lda, doublecomplex a[], double work[]);
extern double zlanhb(char norm, char uplo, int n, int k, int ldab, doublecomplex ab[], double work[]);
extern double zlanhp(char norm, char uplo, int n, doublecomplex ap[], double work[]);
extern double zlanht(char norm, int n, double d[], doublecomplex e[]);
extern double zlantr(char norm, char uplo, char diag, int m, int n, int lda, doublecomplex a[], double work[]);
#endif

/*
 * D2a. Solution of systems of linear equations (real nonsymmetric matrices)
 */
#if defined(_VLARRAY)
extern void dgesv(int n, int nrhs, int lda, double a[][lda], int ipiv[], int ldb, double b[][ldb], int *info);
extern void dgetrf(int m, int n, int lda, double a[][lda], int ipiv[], int *info);
extern void dgetrs(char trans, int n, int nrhs, int lda, double a[][lda], int ipiv[], int ldb, double b[][ldb], int *info);
extern void dgetri(int n, int lda, double a[][lda], int ipiv[], double work[], int lwork, int *info);
extern void dgesvx(char fact, char trans, int n, int nrhs, int lda, double a[][lda], int ldaf, double af[][ldaf], int ipiv[], char *equed, double r[], double c[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dsgesv(int n, int nrhs, int lda, double a[][lda], int ipiv[], int ldb, double b[][ldb], int ldx, double x[][ldx], double work[], float swork[], int *iter, int *info);
extern void dgecon(char norm, int n, int lda, double a[][lda], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void dgbsv(int n, int kl, int ku, int nrhs, int ldab, double ab[][ldab], int ipiv[], int ldb, double b[][ldb], int *info);
extern void dgbtrf(int m, int n, int kl, int ku, int ldab, double ab[][ldab], int ipiv[], int *info);
extern void dgbtrs(char trans, int n, int kl, int ku, int nrhs, int ldab, double ab[][ldab], int ipiv[], int ldb, double b[][ldb], int *info);
extern void dgbsvx(char fact, char trans, int n, int kl, int ku, int nrhs, int ldab, double ab[][ldab], int ldafb, double afb[][ldafb], int ipiv[], char *equed, double r[], double c[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dgbcon(char norm, int n, int kl, int ku, int ldab, double ab[][ldab], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void dgtsv(int n, int nrhs, double dl[], double d[], double du[], int ldb, double b[][ldb], int *info);
extern void dgttrf(int n, double dl[], double d[], double du[], double du2[], int ipiv[], int *info);
extern void dgttrs(char trans, int n, int nrhs, double dl[], double d[], double du[], double du2[], int ipiv[], int ldb, double b[][ldb], int *info);
extern void dgtsvx(char fact, char trans, int n, int nrhs, double dl[], double d[], double du[], double dlf[], double df[], double duf[], double du2[], int ipiv[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dgtcon(char norm, int n, double dl[], double d[], double du[], double du2[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
#else
extern void dgesv(int n, int nrhs, int lda, double a[], int ipiv[], int ldb, double b[], int *info);
extern void dgetrf(int m, int n, int lda, double a[], int ipiv[], int *info);
extern void dgetrs(char trans, int n, int nrhs, int lda, double a[], int ipiv[], int ldb, double b[], int *info);
extern void dgetri(int n, int lda, double a[], int ipiv[], double work[], int lwork, int *info);
extern void dgesvx(char fact, char trans, int n, int nrhs, int lda, double a[], int ldaf, double af[], int ipiv[], char *equed, double r[], double c[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dsgesv(int n, int nrhs, int lda, double a[], int ipiv[], int ldb, double b[], int ldx, double x[], double work[], float swork[], int *iter, int *info);
extern void dgecon(char norm, int n, int lda, double a[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void dgbsv(int n, int kl, int ku, int nrhs, int ldab, double ab[], int ipiv[], int ldb, double b[], int *info);
extern void dgbtrf(int m, int n, int kl, int ku, int ldab, double ab[], int ipiv[], int *info);
extern void dgbtrs(char trans, int n, int kl, int ku, int nrhs, int ldab, double ab[], int ipiv[], int ldb, double b[], int *info);
extern void dgbsvx(char fact, char trans, int n, int kl, int ku, int nrhs, int ldab, double ab[], int ldafb, double afb[], int ipiv[], char *equed, double r[], double c[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dgbcon(char norm, int n, int kl, int ku, int ldab, double ab[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void dgtsv(int n, int nrhs, double dl[], double d[], double du[], int ldb, double b[], int *info);
extern void dgttrf(int n, double dl[], double d[], double du[], double du2[], int ipiv[], int *info);
extern void dgttrs(char trans, int n, int nrhs, double dl[], double d[], double du[], double du2[], int ipiv[], int ldb, double b[], int *info);
extern void dgtsvx(char fact, char trans, int n, int nrhs, double dl[], double d[], double du[], double dlf[], double df[], double duf[], double du2[], int ipiv[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dgtcon(char norm, int n, double dl[], double d[], double du[], double du2[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
#endif

/*
 * D2a3. Solution of systems of linear equations (real triangular matrices)
 */
#if defined(_VLARRAY)
extern void dtrtrs(char uplo, char trans, char diag, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], int *info);
extern void dtrtri(char uplo, char diag, int n, int lda, double a[][lda], int *info);
extern void dtrcon(char norm, char uplo, char diag, int n, int lda, double a[][lda], double *rcond, double work[], int iwork[], int *info);
extern void dtptrs(char uplo, char trans, char diag, int n, int nrhs, double ap[], int ldb, double b[][ldb], int *info);
extern void dtptri(char uplo, char diag, int n, double ap[], int *info);
extern void dtpcon(char norm, char uplo, char diag, int n, double ap[], double *rcond, double work[], int iwork[], int *info);
extern void dtbtrs(char uplo, char trans, char diag, int n, int kd, int nrhs, int ldab, double ab[][ldab], int ldb, double b[][ldb], int *info);
extern void dtbcon(char norm, char uplo, char diag, int n, int kd, int ldab, double ab[][ldab], double *rcond, double work[], int iwork[], int *info);
#else
extern void dtrtrs(char uplo, char trans, char diag, int n, int nrhs, int lda, double a[], int ldb, double b[], int *info);
extern void dtrtri(char uplo, char diag, int n, int lda, double a[], int *info);
extern void dtrcon(char norm, char uplo, char diag, int n, int lda, double a[], double *rcond, double work[], int iwork[], int *info);
extern void dtptrs(char uplo, char trans, char diag, int n, int nrhs, double ap[], int ldb, double b[], int *info);
extern void dtptri(char uplo, char diag, int n, double ap[], int *info);
extern void dtpcon(char norm, char uplo, char diag, int n, double ap[], double *rcond, double work[], int iwork[], int *info);
extern void dtbtrs(char uplo, char trans, char diag, int n, int kd, int nrhs, int ldab, double ab[], int ldb, double b[], int *info);
extern void dtbcon(char norm, char uplo, char diag, int n, int kd, int ldab, double ab[], double *rcond, double work[], int iwork[], int *info);
#endif

/*
 * D2b1a. Solution of systems of linear equations (real symmetric indefinite matrices)
 */
#if defined(_VLARRAY)
extern void dsysv(char uplo, int n, int nrhs, int lda, double a[][lda], int ipiv[], int ldb, double b[][ldb], double work[], int lwork, int *info);
extern void dsytrf(char uplo, int n, int lda, double a[][lda], int ipiv[], double work[], int lwork, int *info);
extern void dsytrs(char uplo, int n, int nrhs, int lda, double a[][lda], int ipiv[], int ldb, double b[][ldb], int *info);
extern void dsytri(char uplo, int n, int lda, double a[][lda], int ipiv[], double work[], int *info);
extern void dsysvx(char fact, char uplo, int n, int nrhs, int lda, double a[][lda], int ldaf, double af[][ldaf], int ipiv[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int lwork, int iwork[], int *info);
extern void dsycon(char uplo, int n, int lda, double a[][lda], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void dspsv(char uplo, int n, int nrhs, double ap[], int ipiv[], int ldb, double b[][ldb], int *info);
extern void dsptrf(char uplo, int n, double ap[], int ipiv[], int *info);
extern void dsptrs(char uplo, int n, int nrhs, double ap[], int ipiv[], int ldb, double b[][ldb], int *info);
extern void dsptri(char uplo, int n, double ap[], int ipiv[], double work[], int *info);
extern void dspsvx(char fact, char uplo, int n, int nrhs, double ap[], double afp[], int ipiv[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dspcon(char uplo, int n, double ap[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
#else
extern void dsysv(char uplo, int n, int nrhs, int lda, double a[], int ipiv[], int ldb, double b[], double work[], int lwork, int *info);
extern void dsytrf(char uplo, int n, int lda, double a[], int ipiv[], double work[], int lwork, int *info);
extern void dsytrs(char uplo, int n, int nrhs, int lda, double a[], int ipiv[], int ldb, double b[], int *info);
extern void dsytri(char uplo, int n, int lda, double a[], int ipiv[], double work[], int *info);
extern void dsysvx(char fact, char uplo, int n, int nrhs, int lda, double a[], int ldaf, double af[], int ipiv[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int lwork, int iwork[], int *info);
extern void dsycon(char uplo, int n, int lda, double a[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void dspsv(char uplo, int n, int nrhs, double ap[], int ipiv[], int ldb, double b[], int *info);
extern void dsptrf(char uplo, int n, double ap[], int ipiv[], int *info);
extern void dsptrs(char uplo, int n, int nrhs, double ap[], int ipiv[], int ldb, double b[], int *info);
extern void dsptri(char uplo, int n, double ap[], int ipiv[], double work[], int *info);
extern void dspsvx(char fact, char uplo, int n, int nrhs, double ap[], double afp[], int ipiv[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dspcon(char uplo, int n, double ap[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
#endif

/*
 * D2b1b. Solution of systems of linear equations (real symmetric positive definite matrices)
 */
#if defined(_VLARRAY)
extern void dposv(char uplo, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], int *info);
extern void dpotrf(char uplo, int n, int lda, double a[][lda], int *info);
extern void dpotrs(char uplo, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], int *info);
extern void dpotri(char uplo, int n, int lda, double a[][lda], int *info);
extern void dposvx(char fact, char uplo, int n, int nrhs, int lda, double a[][lda], int ldaf, double af[][ldaf], char *equed, double s[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dsposv(char uplo, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], int ldx, double x[][ldx], double work[], float swork[], int *iter, int *info);
extern void dpocon(char uplo, int n, int lda, double a[][lda], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void dppsv(char uplo, int n, int nrhs, double ap[], int ldb, double b[][ldb], int *info);
extern void dpptrf(char uplo, int n, double ap[], int *info);
extern void dpptrs(char uplo, int n, int nrhs, double ap[], int ldb, double b[][ldb], int *info);
extern void dpptri(char uplo, int n, double ap[], int *info);
extern void dppsvx(char fact, char uplo, int n, int nrhs, double ap[], double afp[], char *equed, double s[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dppcon(char uplo, int n, double ap[], double anorm, double *rcond, double work[], int iwork[], int *info);
#else
extern void dposv(char uplo, int n, int nrhs, int lda, double a[], int ldb, double b[], int *info);
extern void dpotrf(char uplo, int n, int lda, double a[], int *info);
extern void dpotrs(char uplo, int n, int nrhs, int lda, double a[], int ldb, double b[], int *info);
extern void dpotri(char uplo, int n, int lda, double a[], int *info);
extern void dposvx(char fact, char uplo, int n, int nrhs, int lda, double a[], int ldaf, double af[], char *equed, double s[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dsposv(char uplo, int n, int nrhs, int lda, double a[], int ldb, double b[], int ldx, double x[], double work[], float swork[], int *iter, int *info);
extern void dpocon(char uplo, int n, int lda, double a[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void dppsv(char uplo, int n, int nrhs, double ap[], int ldb, double b[], int *info);
extern void dpptrf(char uplo, int n, double ap[], int *info);
extern void dpptrs(char uplo, int n, int nrhs, double ap[], int ldb, double b[], int *info);
extern void dpptri(char uplo, int n, double ap[], int *info);
extern void dppsvx(char fact, char uplo, int n, int nrhs, double ap[], double afp[], char *equed, double s[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void dppcon(char uplo, int n, double ap[], double anorm, double *rcond, double work[], int iwork[], int *info);
#endif

/*
 * D2b2. Solution of systems of linear equations (real symmetric positive definite band matrices)
 */
#if defined(_VLARRAY)
extern void dpbsv(char uplo, int n, int kd, int nrhs, int ldab, double ab[][ldab], int ldb, double b[][ldb], int *info);
extern void dpbtrf(char uplo, int n, int kd, int ldab, double ab[][ldab], int *info);
extern void dpbtrs(char uplo, int n, int kd, int nrhs, int ldab, double ab[][ldab], int ldb, double b[][ldb], int *info);
extern void dpbsvx(char fact, char uplo, int n, int kd, int nrhs, int ldab, double ab[][ldab], int ldafb, double afb[][ldafb], char *equed, double s[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[],  double berr[], double work[], int iwork[], int *info);
extern void dpbcon(char uplo, int n, int kd, int ldab, double ab[][ldab], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void dptsv(int n, int nrhs, double d[], double e[], int ldb, double b[][ldb], int *info);
extern void dpttrf(int n, double d[], double e[], int *info);
extern void dpttrs(int n, int nrhs, double d[], double e[], int ldb, double b[][ldb], int *info);
extern void dptsvx(char fact, int n, int nrhs, double d[], double e[], double df[], double ef[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int *info);
extern void dptcon(int n, double d[], double e[], double anorm, double *rcond, double work[], int *info);
#else
extern void dpbsv(char uplo, int n, int kd, int nrhs, int ldab, double ab[], int ldb, double b[], int *info);
extern void dpbtrf(char uplo, int n, int kd, int ldab, double ab[], int *info);
extern void dpbtrs(char uplo, int n, int kd, int nrhs, int ldab, double ab[], int ldb, double b[], int *info);
extern void dpbsvx(char fact, char uplo, int n, int kd, int nrhs, int ldab, double ab[], int ldafb, double afb[], char *equed, double s[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[],  double berr[], double work[], int iwork[], int *info);
extern void dpbcon(char uplo, int n, int kd, int ldab, double ab[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void dptsv(int n, int nrhs, double d[], double e[], int ldb, double b[], int *info);
extern void dpttrf(int n, double d[], double e[], int *info);
extern void dpttrs(int n, int nrhs, double d[], double e[], int ldb, double b[], int *info);
extern void dptsvx(char fact, int n, int nrhs, double d[], double e[], double df[], double ef[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int *info);
extern void dptcon(int n, double d[], double e[], double anorm, double *rcond, double work[], int *info);
#endif

/*
 * D2c. Solution of systems of linear equations (complex non-Hermitian matrices)
 */
#if defined(_VLARRAY)
extern void zgesv(int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void zgetrf(int m, int n, int lda, doublecomplex a[][lda], int ipiv[], int *info);
extern void zgetrs(char trans, int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void zgetri(int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void zgesvx(char fact, char trans, int n, int nrhs, int lda, doublecomplex a[][lda], int ldaf, doublecomplex af[][ldaf], int ipiv[], char *equed, double r[], double c[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zcgesv(int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], doublecomplex work[], floatcomplex swork[], double rwork[], int *iter, int *info);
extern void zgecon(char norm, int n, int lda, doublecomplex a[][lda], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void zgbsv(int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[][ldab], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void zgbtrf(int m, int n, int kl, int ku, int ldab, doublecomplex ab[][ldab], int ipiv[], int *info);
extern void zgbtrs(char trans, int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[][ldab], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void zgbsvx(char fact, char trans, int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[][ldab], int ldafb, doublecomplex afb[][ldafb], int ipiv[], char *equed, double r[], double c[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zgbcon(char norm, int n, int kl, int ku, int ldab, doublecomplex ab[][ldab], int ipiv[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void zgtsv(int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], int ldb, doublecomplex b[][ldb], int *info);
extern void zgttrf(int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], int *info);
extern void zgttrs(char trans, int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void zgtsvx(char fact, char trans, int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex dlf[], doublecomplex df[], doublecomplex duf[], doublecomplex du2[], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zgtcon(char norm, int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void zsysv(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], doublecomplex work[], int lwork, int *info);
extern void zsytrf(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void zsytrs(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void zsytri(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex work[], int *info);
extern void zsysvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldaf, doublecomplex af[][ldaf], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zsycon(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void zspsv(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void zsptrf(char uplo, int n, doublecomplex ap[], int ipiv[], int *info);
extern void zsptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void zsptri(char uplo, int n, doublecomplex ap[], int ipiv[], doublecomplex work[], int *info);
extern void zspsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zspcon(char uplo, int n, doublecomplex ap[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
#else
extern void zgesv(int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void zgetrf(int m, int n, int lda, doublecomplex a[], int ipiv[], int *info);
extern void zgetrs(char trans, int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void zgetri(int n, int lda, doublecomplex a[], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void zgesvx(char fact, char trans, int n, int nrhs, int lda, doublecomplex a[], int ldaf, doublecomplex af[], int ipiv[], char *equed, double r[], double c[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zcgesv(int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], doublecomplex work[], floatcomplex swork[], double rwork[], int *iter, int *info);
extern void zgecon(char norm, int n, int lda, doublecomplex a[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void zgbsv(int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void zgbtrf(int m, int n, int kl, int ku, int ldab, doublecomplex ab[], int ipiv[], int *info);
extern void zgbtrs(char trans, int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void zgbsvx(char fact, char trans, int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[], int ldafb, doublecomplex afb[], int ipiv[], char *equed, double r[], double c[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zgbcon(char norm, int n, int kl, int ku, int ldab, doublecomplex ab[], int ipiv[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void zgtsv(int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], int ldb, doublecomplex b[], int *info);
extern void zgttrf(int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], int *info);
extern void zgttrs(char trans, int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void zgtsvx(char fact, char trans, int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex dlf[], doublecomplex df[], doublecomplex duf[], doublecomplex du2[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zgtcon(char norm, int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void zsysv(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], doublecomplex work[], int lwork, int *info);
extern void zsytrf(char uplo, int n, int lda, doublecomplex a[], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void zsytrs(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void zsytri(char uplo, int n, int lda, doublecomplex a[], int ipiv[], doublecomplex work[], int *info);
extern void zsysvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldaf, doublecomplex af[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zsycon(char uplo, int n, int lda, doublecomplex a[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void zspsv(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void zsptrf(char uplo, int n, doublecomplex ap[], int ipiv[], int *info);
extern void zsptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void zsptri(char uplo, int n, doublecomplex ap[], int ipiv[], doublecomplex work[], int *info);
extern void zspsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zspcon(char uplo, int n, doublecomplex ap[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
#endif

/*
 * D2c3. Solution of systems of linear equations (complex triangular matrices)
 */
#if defined(_VLARRAY)
extern void ztrtrs(char uplo, char trans, char diag, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *info);
extern void ztrtri(char uplo, char diag, int n, int lda, doublecomplex a[][lda], int *info);
extern void ztrcon(char norm, char uplo, char diag, int n, int lda, doublecomplex a[][lda], double *rcond, doublecomplex work[], double rwork[], int *info);
extern void ztptrs(char uplo, char trans, char diag, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[][ldb], int *info);
extern void ztptri(char uplo, char diag, int n, doublecomplex ap[], int *info);
extern void ztpcon(char norm, char uplo, char diag, int n, doublecomplex ap[], double *rcond, doublecomplex work[], double rwork[], int *info);
extern void ztbtrs(char uplo, char trans, char diag, int n, int kd, int nrhs, int ldab, doublecomplex ab[][ldab], int ldb, doublecomplex b[][ldb], int *info);
extern void ztbcon(char norm, char uplo, char diag, int n, int kd, int ldab, doublecomplex ab[][ldab], double *rcond, doublecomplex work[], double rwork[], int *info);
#else
extern void ztrtrs(char uplo, char trans, char diag, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *info);
extern void ztrtri(char uplo, char diag, int n, int lda, doublecomplex a[], int *info);
extern void ztrcon(char norm, char uplo, char diag, int n, int lda, doublecomplex a[], double *rcond, doublecomplex work[], double rwork[], int *info);
extern void ztptrs(char uplo, char trans, char diag, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[], int *info);
extern void ztptri(char uplo, char diag, int n, doublecomplex ap[], int *info);
extern void ztpcon(char norm, char uplo, char diag, int n, doublecomplex ap[], double *rcond, doublecomplex work[], double rwork[], int *info);
extern void ztbtrs(char uplo, char trans, char diag, int n, int kd, int nrhs, int ldab, doublecomplex ab[], int ldb, doublecomplex b[], int *info);
extern void ztbcon(char norm, char uplo, char diag, int n, int kd, int ldab, doublecomplex ab[], double *rcond, doublecomplex work[], double rwork[], int *info);
#endif

/*
 * D2d1a. Solution of systems of linear equations (Hermitian matrices)
 */
#if defined(_VLARRAY)
extern void zhesv(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], doublecomplex work[], int lwork, int *info);
extern void zhetrf(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void zhetrs(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void zhetri(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex work[], int *info);
extern void zhesvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldaf, doublecomplex af[][ldaf], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zhecon(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void zhpsv(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void zhptrf(char uplo, int n, doublecomplex ap[], int ipiv[], int *info);
extern void zhptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void zhptri(char uplo, int n, doublecomplex ap[], int ipiv[], doublecomplex work[], int *info);
extern void zhpsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zhpcon(char uplo, int n, doublecomplex ap[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
#else
extern void zhesv(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], doublecomplex work[], int lwork, int *info);
extern void zhetrf(char uplo, int n, int lda, doublecomplex a[], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void zhetrs(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void zhetri(char uplo, int n, int lda, doublecomplex a[], int ipiv[], doublecomplex work[], int *info);
extern void zhesvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldaf, doublecomplex af[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zhecon(char uplo, int n, int lda, doublecomplex a[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void zhpsv(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void zhptrf(char uplo, int n, doublecomplex ap[], int ipiv[], int *info);
extern void zhptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void zhptri(char uplo, int n, doublecomplex ap[], int ipiv[], doublecomplex work[], int *info);
extern void zhpsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zhpcon(char uplo, int n, doublecomplex ap[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
#endif

/*
 * D2d1b. Solution of systems of linear equations (Hermitian positive definite matrices)
 */
#if defined(_VLARRAY)
extern void zposv(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *info);
extern void zpotrf(char uplo, int n, int lda, doublecomplex a[][lda], int *info);
extern void zpotrs(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *info);
extern void zpotri(char uplo, int n, int lda, doublecomplex a[][lda], int *info);
extern void zposvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldaf, doublecomplex af[][ldaf], char *equed, double s[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zcposv(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], doublecomplex work[], floatcomplex swork[], double rwork[], int *iter, int *info);
extern void zpocon(char uplo, int n, int lda, doublecomplex a[][lda], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void zppsv(char uplo, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[][ldb], int *info);
extern void zpptrf(char uplo, int n, doublecomplex ap[], int *info);
extern void zpptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[][ldb], int *info);
extern void zpptri(char uplo, int n, doublecomplex ap[], int *info);
extern void zppsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], char *equed, double s[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zppcon(char uplo, int n, doublecomplex ap[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
#else
extern void zposv(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *info);
extern void zpotrf(char uplo, int n, int lda, doublecomplex a[], int *info);
extern void zpotrs(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *info);
extern void zpotri(char uplo, int n, int lda, doublecomplex a[], int *info);
extern void zposvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldaf, doublecomplex af[], char *equed, double s[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zcposv(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], doublecomplex work[], floatcomplex swork[], double rwork[], int *iter, int *info);
extern void zpocon(char uplo, int n, int lda, doublecomplex a[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void zppsv(char uplo, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[], int *info);
extern void zpptrf(char uplo, int n, doublecomplex ap[], int *info);
extern void zpptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[], int *info);
extern void zpptri(char uplo, int n, doublecomplex ap[], int *info);
extern void zppsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], char *equed, double s[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zppcon(char uplo, int n, doublecomplex ap[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
#endif

/*
 * D2d2. Solution of systems of linear equations (Hermitian positive definite band matrices)
 */
#if defined(_VLARRAY)
extern void zpbsv(char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[][ldab], int ldb, doublecomplex b[][ldb], int *info);
extern void zpbtrf(char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], int *info);
extern void zpbtrs(char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[][ldab], int ldb, doublecomplex b[][ldb], int *info);
extern void zpbsvx(char fact, char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[][ldab], int ldafb, doublecomplex afb[][ldafb], char *equed, double s[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[],  double berr[], doublecomplex work[], double rwork[], int *info);
extern void zpbcon(char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void zptsv(int n, int nrhs, double d[], doublecomplex e[], int ldb, doublecomplex b[][ldb], int *info);
extern void zpttrf(int n, double d[], doublecomplex e[], int *info);
extern void zpttrs(char uplo, int n, int nrhs, double d[], doublecomplex e[], int ldb, doublecomplex b[][ldb], int *info);
extern void zptsvx(char fact, int n, int nrhs, double d[], doublecomplex e[], double df[], doublecomplex ef[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zptcon(int n, double d[], doublecomplex e[], double anorm, double *rcond, double rwork[], int *info);
#else
extern void zpbsv(char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[], int ldb, doublecomplex b[], int *info);
extern void zpbtrf(char uplo, int n, int kd, int ldab, doublecomplex ab[], int *info);
extern void zpbtrs(char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[], int ldb, doublecomplex b[], int *info);
extern void zpbsvx(char fact, char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[], int ldafb, doublecomplex afb[], char *equed, double s[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[],  double berr[], doublecomplex work[], double rwork[], int *info);
extern void zpbcon(char uplo, int n, int kd, int ldab, doublecomplex ab[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void zptsv(int n, int nrhs, double d[], doublecomplex e[], int ldb, doublecomplex b[], int *info);
extern void zpttrf(int n, double d[], doublecomplex e[], int *info);
extern void zpttrs(char uplo, int n, int nrhs, double d[], doublecomplex e[], int ldb, doublecomplex b[], int *info);
extern void zptsvx(char fact, int n, int nrhs, double d[], doublecomplex e[], double df[], doublecomplex ef[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void zptcon(int n, double d[], doublecomplex e[], double anorm, double *rcond, double rwork[], int *info);
#endif

/*
 * D4a1. Real symmetric matrix eigenvalue problems
 */
#if defined(_VLARRAY)
extern void dsyev(char jobz, char uplo, int n,  int lda, double a[][lda], double w[], double work[], int lwork, int *info);
extern void dsyevd(char jobz, char uplo, int n, int lda, double a[][lda], double w[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dsyevr(char jobz, char range, char uplo, int n, int lda, double a[][lda], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], int isuppz[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dsyevx(char jobz, char range, char uplo, int n, int lda, double a[][lda], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int ifail[], int *info);
extern void dspev(char jobz, char uplo, int n, double ap[], double w[], int ldz, double z[][ldz], double work[], int *info);
extern void dspevd(char jobz, char uplo, int n, double ap[], double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dspevx(char jobz, char range, char uplo, int n, double ap[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void dsbev(char jobz, char uplo, int n, int kd, int ldab, double ab[][ldab], double w[], int ldz, double z[][ldz], double work[], int *info);
extern void dsbevd(char jobz, char uplo, int n, int kd, int ldab, double ab[][ldab], double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dsbevx(char jobz, char range, char uplo, int n, int kd, int ldab, double ab[][ldab], int ldq, double q[][ldq], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void dsbtrd(char vect, char uplo, int n, int kd, int ldab, double ab[][ldab], double d[], double e[], int ldq, double q[][ldq], double work[], int *info);
extern void dstev(char jobz, int n, double d[], double e[], int ldz, double z[][ldz], double work[], int *info);
extern void dstevd(char jobz, int n, double d[], double e[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dstevr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], int isuppz[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dstevx(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void dsytrd(char uplo, int n, int lda, double a[][lda], double d[], double e[], double tau[], double work[], int lwork, int *info);
extern void dorgtr(char uplo, int n, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void dormtr(char side, char uplo, char trans, int m, int n, int lda, double a[][lda], double tau[], int ldc, double c[][ldc], double work[], int lwork, int *info);
extern void dopgtr(char uplo, int n, double ap[], double tau[], int ldq, double q[][ldq],  double work[], int *info);
extern void dopmtr(char side, char uplo, char trans, int m, int n, double ap[], double tau[], int ldc, double c[][ldc], double work[], int *info);
extern void dsteqr(char compz, int n, double d[], double e[], int ldz, double z[][ldz], double work[], int *info);
extern void dstedc(char compz, int n, double d[], double e[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dstemr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, int *m, double w[], int ldz, double z[][ldz], int nzc, int isuppz[], int *tryrac, double work[], int lwork, int iwork[], int liwork, int *info);
extern void dstein(int n, double d[], double e[], int m, double w[], int iblock[], int isplit[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void dpteqr(char compz, int n, double d[], double e[], int ldz, double z[][ldz], double work[], int *info);
#else
extern void dsyev(char jobz, char uplo, int n,  int lda, double a[], double w[], double work[], int lwork, int *info);
extern void dsyevd(char jobz, char uplo, int n, int lda, double a[], double w[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dsyevr(char jobz, char range, char uplo, int n, int lda, double a[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], int isuppz[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dsyevx(char jobz, char range, char uplo, int n, int lda, double a[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int lwork, int iwork[], int ifail[], int *info);
extern void dspev(char jobz, char uplo, int n, double ap[], double w[], int ldz, double z[], double work[], int *info);
extern void dspevd(char jobz, char uplo, int n, double ap[], double w[], int ldz, double z[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dspevx(char jobz, char range, char uplo, int n, double ap[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
extern void dsbev(char jobz, char uplo, int n, int kd, int ldab, double ab[], double w[], int ldz, double z[], double work[], int *info);
extern void dsbevd(char jobz, char uplo, int n, int kd, int ldab, double ab[], double w[], int ldz, double z[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dsbevx(char jobz, char range, char uplo, int n, int kd, int ldab, double ab[], int ldq, double q[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
extern void dsbtrd(char vect, char uplo, int n, int kd, int ldab, double ab[], double d[], double e[], int ldq, double q[], double work[], int *info);
extern void dstev(char jobz, int n, double d[], double e[], int ldz, double z[], double work[], int *info);
extern void dstevd(char jobz, int n, double d[], double e[], int ldz, double z[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dstevr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], int isuppz[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dstevx(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
extern void dsytrd(char uplo, int n, int lda, double a[], double d[], double e[], double tau[], double work[], int lwork, int *info);
extern void dorgtr(char uplo, int n, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void dormtr(char side, char uplo, char trans, int m, int n, int lda, double a[], double tau[], int ldc, double c[], double work[], int lwork, int *info);
extern void dopgtr(char uplo, int n, double ap[], double tau[], int ldq, double q[],  double work[], int *info);
extern void dopmtr(char side, char uplo, char trans, int m, int n, double ap[], double tau[], int ldc, double c[], double work[], int *info);
extern void dsteqr(char compz, int n, double d[], double e[], int ldz, double z[], double work[], int *info);
extern void dstedc(char compz, int n, double d[], double e[], int ldz, double z[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dstemr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, int *m, double w[], int ldz, double z[], int nzc, int isuppz[], int *tryrac, double work[], int lwork, int iwork[], int liwork, int *info);
extern void dstein(int n, double d[], double e[], int m, double w[], int iblock[], int isplit[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
extern void dpteqr(char compz, int n, double d[], double e[], int ldz, double z[], double work[], int *info);
#endif
extern void dsterf(int n, double d[], double e[], int *info);
extern void dstebz(char range, char order, int n, double vl, double vu, int il, int iu, double abstol, double d[], double e[], int *m, int *nsplit, double w[], int iblock[], int isplit[], double work[], int iwork[], int *info);
extern void dsptrd(char uplo, int n, double ap[], double d[], double e[], double tau[], int *info);
extern void ddisna(char job, int m, int n, double d[], double sep[], int *info);

/*
 * D4a2. Real nonsymmetric matrix eigenvalue problems
 */
#if defined(_VLARRAY)
extern void dgeev(char jobvl, char jobvr, int n, int lda, double a[][lda], double wr[], double wi[], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], double work[], int lwork, int *info);
extern void dgeevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, double a[][lda], double wr[], double wi[], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], int *ilo, int *ihi, double scale[], double *abnrm, double rconde[], double rcondv[], double work[], int lwork, int iwork[], int *info);
extern void dgehrd(int n, int ilo, int ihi, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void dgebal(char job, int n, int lda, double a[][lda], int *ilo, int *ihi, double scale[], int *info);
extern void dgebak(char job, char side, int n, int ilo, int ihi, double scale[], int m, int ldv, double v[][ldv], int *info);
extern void dorghr(int n, int ilo, int ihi, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void dormhr(char side, char trans, int m, int n, int ilo, int ihi, int lda, double a[][lda], double tau[], int ldc, double c[][ldc], double work[], int lwork, int *info);
extern void dhseqr(char job, char compz, int n, int ilo, int ihi, int ldh, double h[][ldh], double wr[], double wi[], int ldz, double z[][ldz], double work[], int lwork, int *info);
extern void dhsein(char side, char eigsrc, char initv, int select[], int n, int ldh, double h[][ldh], double wr[], double wi[], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], int mm, int *m, double work[], int ifaill[], int ifailr[], int *info);
extern void dtrevc3(char side, char howmny, int select[], int n, int ldt, double t[][ldt], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], int mm, int *m, double work[], int lwork, int *info);
extern void dtrexc(char compq, int n, int ldt, double t[][ldt], int ldq, double q[][ldq], int *ifst, int *ilst, double work[], int *info);
extern void dtrsyl(char transa, char transb, int isgn, int m, int n, int lda, double a[][lda], int ldb, double b[][ldb], int ldc, double c[][ldc], double *scale, int *info);
extern void dtrsna(char job, char howmny, int select[], int n, int ldt, double t[][ldt], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], double s[], double sep[], int mm, int *m, double work[], int lwork, int iwork[], int liwork, int *info);
extern void dtrsen(char job, char compq, int select[], int n, int ldt, double t[][ldt], int ldq, double q[][ldq], double wr[], double wi[], int *m, double *s, double *sep, double work[], int lwork, int iwork[], int liwork, int *info);
extern void dgees(char jobvs, char sort, int (*select)(double, double), int n, int lda, double a[][lda], int *sdim, double wr[], double wi[], int ldvs, double vs[][ldvs], double work[], int lwork, int bwork[], int *info);
extern void dgeesx(char jobvs, char sort, int (*select)(double, double), char sense, int n, int lda, double a[][lda], int *sdim, double wr[], double wi[], int ldvs, double vs[][ldvs], double *rconde, double *rcondv, double work[], int lwork, int iwork[], int liwork, int bwork[], int *info);
extern void dgees_r(char jobvs, char sort, int n, int lda, double a[][lda], int *sdim, double wr[], double wi[], int ldvs, double vs[][ldvs], double work[], int lwork, int bwork[], int *info, int *irev);
extern void dgeesx_r(char jobvs, char sort, char sense, int n, int lda, double a[][lda], int *sdim, double wr[], double wi[], int ldvs, double vs[][ldvs], double *rconde, double *rcondv, double work[], int lwork, int iwork[], int liwork, int bwork[], int *info, int *irev);
#else
extern void dgeev(char jobvl, char jobvr, int n, int lda, double a[], double wr[], double wi[], int ldvl, double vl[], int ldvr, double vr[], double work[], int lwork, int *info);
extern void dgeevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, double a[], double wr[], double wi[], int ldvl, double vl[], int ldvr, double vr[], int *ilo, int *ihi, double scale[], double *abnrm, double rconde[], double rcondv[], double work[], int lwork, int iwork[], int *info);
extern void dgehrd(int n, int ilo, int ihi, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void dgebal(char job, int n, int lda, double a[], int *ilo, int *ihi, double scale[], int *info);
extern void dgebak(char job, char side, int n, int ilo, int ihi, double scale[], int m, int ldv, double v[], int *info);
extern void dorghr(int n, int ilo, int ihi, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void dormhr(char side, char trans, int m, int n, int ilo, int ihi, int lda, double a[], double tau[], int ldc, double c[], double work[], int lwork, int *info);
extern void dhseqr(char job, char compz, int n, int ilo, int ihi, int ldh, double h[], double wr[], double wi[], int ldz, double z[], double work[], int lwork, int *info);
extern void dhsein(char side, char eigsrc, char initv, int select[], int n, int ldh, double h[], double wr[], double wi[], int ldvl, double vl[], int ldvr, double vr[], int mm, int *m, double work[], int ifaill[], int ifailr[], int *info);
extern void dtrevc3(char side, char howmny, int select[], int n, int ldt, double t[], int ldvl, double vl[], int ldvr, double vr[], int mm, int *m, double work[], int lwork, int *info);
extern void dtrexc(char compq, int n, int ldt, double t[], int ldq, double q[], int *ifst, int *ilst, double work[], int *info);
extern void dtrsyl(char transa, char transb, int isgn, int m, int n, int lda, double a[], int ldb, double b[], int ldc, double c[], double *scale, int *info);
extern void dtrsna(char job, char howmny, int select[], int n, int ldt, double t[], int ldvl, double vl[], int ldvr, double vr[], double s[], double sep[], int mm, int *m, double work[], int lwork, int iwork[], int liwork, int *info);
extern void dtrsen(char job, char compq, int select[], int n, int ldt, double t[], int ldq, double q[], double wr[], double wi[], int *m, double *s, double *sep, double work[], int lwork, int iwork[], int liwork, int *info);
extern void dgees(char jobvs, char sort, int (*select)(double, double), int n, int lda, double a[], int *sdim, double wr[], double wi[], int ldvs, double vs[], double work[], int lwork, int bwork[], int *info);
extern void dgeesx(char jobvs, char sort, int (*select)(double, double), char sense, int n, int lda, double a[], int *sdim, double wr[], double wi[], int ldvs, double vs[], double *rconde, double *rcondv, double work[], int lwork, int iwork[], int liwork, int bwork[], int *info);
extern void dgees_r(char jobvs, char sort, int n, int lda, double a[], int *sdim, double wr[], double wi[], int ldvs, double vs[], double work[], int lwork, int bwork[], int *info, int *irev);
extern void dgeesx_r(char jobvs, char sort, char sense, int n, int lda, double a[], int *sdim, double wr[], double wi[], int ldvs, double vs[], double *rconde, double *rcondv, double work[], int lwork, int iwork[], int liwork, int bwork[], int *info, int *irev);
#endif

/*
 * D4a3. Hermitian matrix eigenvalue problems
 */
#if defined(_VLARRAY)
extern void zheev(char jobz, char uplo, int n,  int lda, doublecomplex a[][lda], double w[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zheevd(char jobz, char uplo, int n, int lda, doublecomplex a[][lda], double w[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zheevr(char jobz, char range, char uplo, int n, int lda, doublecomplex a[][lda], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], int isuppz[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zheevx(char jobz, char range, char uplo, int n, int lda, doublecomplex a[][lda], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int iwork[], int ifail[], int *info);
extern void zhetrd(char uplo, int n, int lda, doublecomplex a[][lda], double d[], double e[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zungtr(char uplo, int n, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zunmtr(char side, char uplo, char trans, int m, int n, int lda, doublecomplex a[][lda], doublecomplex tau[], int ldc, doublecomplex c[][ldc], doublecomplex work[], int lwork, int *info);
extern void zupgtr(char uplo, int n, doublecomplex ap[], doublecomplex tau[], int ldq, doublecomplex q[][ldq], doublecomplex work[], int *info);
extern void zupmtr(char side, char uplo, char trans, int m, int n, doublecomplex ap[], doublecomplex tau[], int ldc, doublecomplex c[][ldc], doublecomplex work[], int *info);
extern void zsteqr(char compz, int n, double d[], double e[], int ldz, doublecomplex z[][ldz], double work[], int *info);
extern void zstedc(char compz, int n, double d[], double e[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zstemr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, int *m, double w[], int ldz, doublecomplex z[][ldz], int nzc, int isuppz[], int *tryrac, double work[], int lwork, int iwork[], int liwork, int *info);
extern void zstein(int n, double d[], double e[], int m, double w[], int iblock[], int isplit[], int ldz, doublecomplex z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void zpteqr(char compz, int n, double d[], double e[], int ldz, doublecomplex z[][ldz], double work[], int *info);
extern void zhpev(char jobz, char uplo, int n, doublecomplex ap[], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int *info);
extern void zhpevd(char jobz, char uplo, int n, doublecomplex ap[], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zhpevx(char jobz, char range, char uplo, int n, doublecomplex ap[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void zhbev(char jobz, char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int *info);
extern void zhbevd(char jobz, char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zhbevx(char jobz, char range, char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], int ldq, doublecomplex q[][ldq], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void zhbtrd(char vect, char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], double d[], double e[], int ldq, doublecomplex q[][ldq], doublecomplex work[], int *info);
#else
extern void zheev(char jobz, char uplo, int n,  int lda, doublecomplex a[], double w[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zheevd(char jobz, char uplo, int n, int lda, doublecomplex a[], double w[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zheevr(char jobz, char range, char uplo, int n, int lda, doublecomplex a[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], int isuppz[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zheevx(char jobz, char range, char uplo, int n, int lda, doublecomplex a[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int iwork[], int ifail[], int *info);
extern void zhetrd(char uplo, int n, int lda, doublecomplex a[], double d[], double e[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zungtr(char uplo, int n, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zunmtr(char side, char uplo, char trans, int m, int n, int lda, doublecomplex a[], doublecomplex tau[], int ldc, doublecomplex c[], doublecomplex work[], int lwork, int *info);
extern void zupgtr(char uplo, int n, doublecomplex ap[], doublecomplex tau[], int ldq, doublecomplex q[], doublecomplex work[], int *info);
extern void zupmtr(char side, char uplo, char trans, int m, int n, doublecomplex ap[], doublecomplex tau[], int ldc, doublecomplex c[], doublecomplex work[], int *info);
extern void zsteqr(char compz, int n, double d[], double e[], int ldz, doublecomplex z[], double work[], int *info);
extern void zstedc(char compz, int n, double d[], double e[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zstemr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, int *m, double w[], int ldz, doublecomplex z[], int nzc, int isuppz[], int *tryrac, double work[], int lwork, int iwork[], int liwork, int *info);
extern void zstein(int n, double d[], double e[], int m, double w[], int iblock[], int isplit[], int ldz, doublecomplex z[], double work[], int iwork[], int ifail[], int *info);
extern void zpteqr(char compz, int n, double d[], double e[], int ldz, doublecomplex z[], double work[], int *info);
extern void zhpev(char jobz, char uplo, int n, doublecomplex ap[], double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int *info);
extern void zhpevd(char jobz, char uplo, int n, doublecomplex ap[], double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zhpevx(char jobz, char range, char uplo, int n, doublecomplex ap[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void zhbev(char jobz, char uplo, int n, int kd, int ldab, doublecomplex ab[], double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int *info);
extern void zhbevd(char jobz, char uplo, int n, int kd, int ldab, doublecomplex ab[], double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zhbevx(char jobz, char range, char uplo, int n, int kd, int ldab, doublecomplex ab[], int ldq, doublecomplex q[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void zhbtrd(char vect, char uplo, int n, int kd, int ldab, doublecomplex ab[], double d[], double e[], int ldq, doublecomplex q[], doublecomplex work[], int *info);
#endif
extern void zhptrd(char uplo, int n, doublecomplex ap[], double d[], double e[], doublecomplex tau[], int *info);
/*extern void ddisna(char job, int m, int n, double d[], double sep[], int *info);*/

/*
 * D4a4. Complex nonsymmetric matrix eigenvalue problems
 */
#if defined(_VLARRAY)
extern void zgeev(char jobvl, char jobvr, int n, int lda, doublecomplex a[][lda], doublecomplex w[], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgeevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, doublecomplex a[][lda], doublecomplex w[], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], int *ilo, int *ihi, double scale[], double *abnrm, double rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgehrd(int n, int ilo, int ihi, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zgebal(char job, int n, int lda, doublecomplex a[][lda], int *ilo, int *ihi, double scale[], int *info);
extern void zgebak(char job, char side, int n, int ilo, int ihi, double scale[], int m, int ldv, doublecomplex v[][ldv], int *info);
extern void zunghr(int n, int ilo, int ihi, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zunmhr(char side, char trans, int m, int n, int ilo, int ihi, int lda, doublecomplex a[][lda], doublecomplex tau[], int ldc, doublecomplex c[][ldc], doublecomplex work[], int lwork, int *info);
extern void zhseqr(char job, char compz, int n, int ilo, int ihi, int ldh, doublecomplex h[][ldh], doublecomplex w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, int *info);
extern void zhsein(char side, char eigsrc, char initv, int select[], int n, int ldh, doublecomplex h[][ldh], doublecomplex w[], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], int mm, int *m, doublecomplex work[], double rwork[], int ifaill[], int ifailr[], int *info);
extern void ztrevc3(char side, char howmny, int select[], int n, int ldt, doublecomplex t[][ldt], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], int mm, int *m, doublecomplex work[], int lwork, double rwork[], int lrwork, int *info);
extern void ztrexc(char compq, int n, int ldt, doublecomplex t[][ldt], int ldq, doublecomplex q[][ldq], int ifst, int ilst, int *info);
extern void ztrsyl(char transa, char transb, int isgn, int m, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int ldc, doublecomplex c[][ldc], double *scale, int *info);
extern void ztrsna(char job, char howmny, int select[], int n, int ldt, doublecomplex t[][ldt], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], double s[], double sep[], int mm, int *m, doublecomplex work[], int lwork, double rwork[], int *info);
extern void ztrsen(char job, char compq, int select[], int n, int ldt, doublecomplex t[][ldt], int ldq, doublecomplex q[][ldq], doublecomplex w[], int *m, double *s, double *sep, doublecomplex work[], int lwork, int *info);
extern void zgees(char jobvs, char sort, int (*select)(doublecomplex), int n, int lda, doublecomplex a[][lda], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[][ldvs], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void zgeesx(char jobvs, char sort, int (*select)(doublecomplex), char sense, int n, int lda, doublecomplex a[][lda], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[][ldvs], double *rconde, double *rcondv, doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void zgees_r(char jobvs, char sort, int n, int lda, doublecomplex a[][lda], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[][ldvs], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
extern void zgeesx_r(char jobvs, char sort, char sense, int n, int lda, doublecomplex a[][lda], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[][ldvs], double *rconde, double *rcondv, doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
#else
extern void zgeev(char jobvl, char jobvr, int n, int lda, doublecomplex a[], doublecomplex w[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgeevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, doublecomplex a[], doublecomplex w[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], int *ilo, int *ihi, double scale[], double *abnrm, double rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgehrd(int n, int ilo, int ihi, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zgebal(char job, int n, int lda, doublecomplex a[], int *ilo, int *ihi, double scale[], int *info);
extern void zgebak(char job, char side, int n, int ilo, int ihi, double scale[], int m, int ldv, doublecomplex v[], int *info);
extern void zunghr(int n, int ilo, int ihi, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zunmhr(char side, char trans, int m, int n, int ilo, int ihi, int lda, doublecomplex a[], doublecomplex tau[], int ldc, doublecomplex c[], doublecomplex work[], int lwork, int *info);
extern void zhseqr(char job, char compz, int n, int ilo, int ihi, int ldh, doublecomplex h[], doublecomplex w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, int *info);
extern void zhsein(char side, char eigsrc, char initv, int select[], int n, int ldh, doublecomplex h[], doublecomplex w[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], int mm, int *m, doublecomplex work[], double rwork[], int ifaill[], int ifailr[], int *info);

extern void ztrevc3(char side, char howmny, int select[], int n, int ldt, doublecomplex t[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], int mm, int *m, doublecomplex work[], int lwork, double rwork[], int lrwork, int *info);
extern void ztrexc(char compq, int n, int ldt, doublecomplex t[], int ldq, doublecomplex q[], int ifst, int ilst, int *info);
extern void ztrsyl(char transa, char transb, int isgn, int m, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], int ldc, doublecomplex c[], double *scale, int *info);

extern void ztrsna(char job, char howmny, int select[], int n, int ldt, doublecomplex t[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], double s[], double sep[], int mm, int *m, doublecomplex work[], int lwork, double rwork[], int *info);
extern void ztrsen(char job, char compq, int select[], int n, int ldt, doublecomplex t[], int ldq, doublecomplex q[], doublecomplex w[], int *m, double *s, double *sep, doublecomplex work[], int lwork, int *info);

extern void zgees(char jobvs, char sort, int (*select)(doublecomplex), int n, int lda, doublecomplex a[], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void zgeesx(char jobvs, char sort, int (*select)(doublecomplex), char sense, int n, int lda, doublecomplex a[], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[], double *rconde, double *rcondv, doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void zgees_r(char jobvs, char sort, int n, int lda, doublecomplex a[], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
extern void zgeesx_r(char jobvs, char sort, char sense, int n, int lda, doublecomplex a[], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[], double *rconde, double *rcondv, doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
#endif

/*
 * D4b1. Real symmetric generalized matrix eigenvalue problems
 */
#if defined(_VLARRAY)
extern void dsygv(int itype, char jobz, char uplo, int n, int lda, double a[][lda], int ldb, double b[][ldb], double w[], double work[], int lwork, int *info);
extern void dsygvd(int itype, char jobz, char uplo, int n, int lda, double a[][lda], int ldb, double b[][ldb], double w[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dsygvx(int itype, char jobz, char range, char uplo, int n, int lda, double a[][lda], int ldb, double b[][ldb], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int ifail[], int *info);
extern void dspgv(int itype, char jobz, char uplo, int n, double ap[], double bp[], double w[], int ldz, double z[][ldz], double work[], int *info);
extern void dspgvd(int itype, char jobz, char uplo, int n, double ap[], double bp[], double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dspgvx(int itype, char jobz, char range, char uplo, int n, double ap[], double bp[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void dsbgv(char jobz, char uplo, int n, int ka, int kb, int ldab, double ab[][ldab], int ldbb, double bb[][ldbb], double w[], int ldz, double z[][ldz], double work[], int *info);
extern void dsbgvd(char jobz, char uplo, int n, int ka, int kb, int ldab, double ab[][ldab], int ldbb, double bb[][ldbb], double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dsbgvx(char jobz, char range, char uplo, int n, int ka, int kb, int ldab, double ab[][ldab], int ldbb, double bb[][ldbb], int ldq, double q[][ldq], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
#else
extern void dsygv(int itype, char jobz, char uplo, int n, int lda, double a[], int ldb, double b[], double w[], double work[], int lwork, int *info);
extern void dsygvd(int itype, char jobz, char uplo, int n, int lda, double a[], int ldb, double b[], double w[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dsygvx(int itype, char jobz, char range, char uplo, int n, int lda, double a[], int ldb, double b[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int lwork, int iwork[], int ifail[], int *info);
extern void dspgv(int itype, char jobz, char uplo, int n, double ap[], double bp[], double w[], int ldz, double z[], double work[], int *info);
extern void dspgvd(int itype, char jobz, char uplo, int n, double ap[], double bp[], double w[], int ldz, double z[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dspgvx(int itype, char jobz, char range, char uplo, int n, double ap[], double bp[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
extern void dsbgv(char jobz, char uplo, int n, int ka, int kb, int ldab, double ab[], int ldbb, double bb[], double w[], int ldz, double z[], double work[], int *info);
extern void dsbgvd(char jobz, char uplo, int n, int ka, int kb, int ldab, double ab[], int ldbb, double bb[], double w[], int ldz, double z[], double work[], int lwork, int iwork[], int liwork, int *info);
extern void dsbgvx(char jobz, char range, char uplo, int n, int ka, int kb, int ldab, double ab[], int ldbb, double bb[], int ldq, double q[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
#endif

/*
 * D4b2. Real generalized matrix eigenvalue problems
 */
#if defined(_VLARRAY)
extern void dggev(char jobvl, char jobvr, int n, int lda, double a[][lda], int ldb, double b[][ldb], double alphar[], double alphai[], double beta[], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], double work[], int lwork, int *info);
extern void dggevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, double a[][lda], int ldb, double b[][ldb], double alphar[], double alphai[], double beta[], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], int *ilo, int *ihi, double lscale[], double rscale[], double *abnrm, double *bbnrm, double rconde[], double rcondv[], double work[], int lwork, int iwork[], int bwork[], int *info);
extern void dgges(char jobvsl, char jobvsr, char sort, int (*selctg)(double, double, double), int n, int lda, double a[][lda], int ldb, double b[][ldb], int *sdim, double alphar[], double alphai[], double beta[], int ldvsl, double vsl[][ldvsl], int ldvsr, double vsr[][ldvsr], double work[], int lwork, int bwork[], int *info);
extern void dggesx(char jobvsl, char jobvsr, char sort, int (*selctg)(double, double, double), char sense, int n, int lda, double a[][lda], int ldb, double b[][ldb], int *sdim, double alphar[], double alphai[], double beta[], int ldvsl, double vsl[][ldvsl], int ldvsr, double vsr[][ldvsr], double rconde[], double rcondv[], double work[], int lwork, int iwork[], int liwork, int bwork[], int *info);
extern void dgges_r(char jobvsl, char jobvsr, char sort, int n, int lda, double a[][lda], int ldb, double b[][ldb], int *sdim, double alphar[], double alphai[], double beta[], int ldvsl, double vsl[][ldvsl], int ldvsr, double vsr[][ldvsr], double work[], int lwork, int bwork[], int *info, int *irev);
extern void dggesx_r(char jobvsl, char jobvsr, char sort, char sense, int n, int lda, double a[][lda], int ldb, double b[][ldb], int *sdim, double alphar[], double alphai[], double beta[], int ldvsl, double vsl[][ldvsl], int ldvsr, double vsr[][ldvsr], double rconde[], double rcondv[], double work[], int lwork, int iwork[], int liwork, int bwork[], int *info, int *irev);
#else
extern void dggev(char jobvl, char jobvr, int n, int lda, double a[], int ldb, double b[], double alphar[], double alphai[], double beta[], int ldvl, double vl[], int ldvr, double vr[], double work[], int lwork, int *info);
extern void dggevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, double a[], int ldb, double b[], double alphar[], double alphai[], double beta[], int ldvl, double vl[], int ldvr, double vr[], int *ilo, int *ihi, double lscale[], double rscale[], double *abnrm, double *bbnrm, double rconde[], double rcondv[], double work[], int lwork, int iwork[], int bwork[], int *info);
extern void dgges(char jobvsl, char jobvsr, char sort, int (*selctg)(double, double, double), int n, int lda, double a[], int ldb, double b[], int *sdim, double alphar[], double alphai[], double beta[], int ldvsl, double vsl[], int ldvsr, double vsr[], double work[], int lwork, int bwork[], int *info);
extern void dggesx(char jobvsl, char jobvsr, char sort, int (*selctg)(double, double, double), char sense, int n, int lda, double a[], int ldb, double b[], int *sdim, double alphar[], double alphai[], double beta[], int ldvsl, double vsl[], int ldvsr, double vsr[], double rconde[], double rcondv[], double work[], int lwork, int iwork[], int liwork, int bwork[], int *info);
extern void dgges_r(char jobvsl, char jobvsr, char sort, int n, int lda, double a[], int ldb, double b[], int *sdim, double alphar[], double alphai[], double beta[], int ldvsl, double vsl[], int ldvsr, double vsr[], double work[], int lwork, int bwork[], int *info, int *irev);
extern void dggesx_r(char jobvsl, char jobvsr, char sort, char sense, int n, int lda, double a[], int ldb, double b[], int *sdim, double alphar[], double alphai[], double beta[], int ldvsl, double vsl[], int ldvsr, double vsr[], double rconde[], double rcondv[], double work[], int lwork, int iwork[], int liwork, int bwork[], int *info, int *irev);
#endif

/*
 * D4b3. Complex Hermitian generalized matrix eigenvalue problems
 */
#if defined(_VLARRAY)
extern void zhegv(int itype, char jobz, char uplo, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double w[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zhegvd(int itype, char jobz, char uplo, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double w[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zhegvx(int itype, char jobz, char range, char uplo, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int iwork[], int ifail[], int *info);
extern void zhpgv(int itype, char jobz, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int *info);
extern void zhpgvd(int itype, char jobz, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zhpgvx(int itype, char jobz, char range, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void zhbgv(char jobz, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[][ldab], int ldbb, doublecomplex bb[][ldbb], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int *info);
extern void zhbgvd(char jobz, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[][ldab], int ldbb, doublecomplex bb[][ldbb], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zhbgvx(char jobz, char range, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[][ldab], int ldbb, doublecomplex bb[][ldbb], int ldq, doublecomplex q[][ldq], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
#else
extern void zhegv(int itype, char jobz, char uplo, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], double w[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zhegvd(int itype, char jobz, char uplo, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], double w[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zhegvx(int itype, char jobz, char range, char uplo, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int iwork[], int ifail[], int *info);
extern void zhpgv(int itype, char jobz, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int *info);
extern void zhpgvd(int itype, char jobz, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zhpgvx(int itype, char jobz, char range, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void zhbgv(char jobz, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[], int ldbb, doublecomplex bb[], double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int *info);
extern void zhbgvd(char jobz, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[], int ldbb, doublecomplex bb[], double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zhbgvx(char jobz, char range, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[], int ldbb, doublecomplex bb[], int ldq, doublecomplex q[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
#endif

/*
 * D4b4. Complex generalized matrix eigenvalue problems
 */
#if defined(_VLARRAY)
extern void zggev(char jobvl, char jobvr, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex alpha[], doublecomplex beta[], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zggevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex alpha[], doublecomplex beta[], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], int *ilo, int *ihi, double lscale[], double rscale[], double *abnrm, double *bbnrm, double rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int bwork[], int *info);
extern void zgges(char jobvsl, char jobvsr, char sort, int (*selctg)(doublecomplex, doublecomplex), int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *sdim, doublecomplex alpha[], doublecomplex beta[], int ldvsl, doublecomplex vsl[][ldvsl], int ldvsr, doublecomplex vsr[][ldvsr], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void zggesx(char jobvsl, char jobvsr, char sort, int (*selctg)(doublecomplex, doublecomplex), char sense, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *sdim, doublecomplex alpha[], doublecomplex beta[], int ldvsl, doublecomplex vsl[][ldvsl], int ldvsr, doublecomplex vsr[][ldvsr], double rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int liwork, int bwork[], int *info);
extern void zgges_r(char jobvsl, char jobvsr, char sort, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *sdim, doublecomplex alpha[], doublecomplex beta[], int ldvsl, doublecomplex vsl[][ldvsl], int ldvsr, doublecomplex vsr[][ldvsr], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
extern void zggesx_r(char jobvsl, char jobvsr, char sort, char sense, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *sdim, doublecomplex alpha[], doublecomplex beta[], int ldvsl, doublecomplex vsl[][ldvsl], int ldvsr, doublecomplex vsr[][ldvsr], double rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int liwork, int bwork[], int *info, int *irev);
#else
extern void zggev(char jobvl, char jobvr, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex alpha[], doublecomplex beta[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zggevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex alpha[], doublecomplex beta[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], int *ilo, int *ihi, double lscale[], double rscale[], double *abnrm, double *bbnrm, double rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int bwork[], int *info);
extern void zgges(char jobvsl, char jobvsr, char sort, int (*selctg)(doublecomplex, doublecomplex), int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *sdim, doublecomplex alpha[], doublecomplex beta[], int ldvsl, doublecomplex vsl[], int ldvsr, doublecomplex vsr[], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void zggesx(char jobvsl, char jobvsr, char sort, int (*selctg)(doublecomplex, doublecomplex), char sense, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *sdim, doublecomplex alpha[], doublecomplex beta[], int ldvsl, doublecomplex vsl[], int ldvsr, doublecomplex vsr[], double rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int liwork, int bwork[], int *info);
extern void zgges_r(char jobvsl, char jobvsr, char sort, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *sdim, doublecomplex alpha[], doublecomplex beta[], int ldvsl, doublecomplex vsl[], int ldvsr, doublecomplex vsr[], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
extern void zggesx_r(char jobvsl, char jobvsr, char sort, char sense, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *sdim, doublecomplex alpha[], doublecomplex beta[], int ldvsl, doublecomplex vsl[], int ldvsr, doublecomplex vsr[], double rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int liwork, int bwork[], int *info, int *irev);
#endif

/*
 * D5. QR factorization
 */
#if defined(_VLARRAY)
extern void dgeqp3(int m, int n, int lda, double a[][lda], int jpvt[], double tau[], double work[], int lwork, int *info);
extern void dgeqrf(int m, int n, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void dorgqr(int m, int n, int k, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void dormqr(char side, char trans, int m, int n, int k, int lda, double a[][lda], double tau[], int ldc, double c[][ldc], double work[], int lwork, int *info);
extern void dgelqf(int m, int n, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void dorglq(int m, int n, int k, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void dormlq(char side, char trans, int m, int n, int k, int lda, double a[][lda], double tau[], int ldc, double c[][ldc], double work[], int lwork, int *info);
#else
extern void dgeqp3(int m, int n, int lda, double a[], int jpvt[], double tau[], double work[], int lwork, int *info);
extern void dgeqrf(int m, int n, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void dorgqr(int m, int n, int k, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void dormqr(char side, char trans, int m, int n, int k, int lda, double a[], double tau[], int ldc, double c[], double work[], int lwork, int *info);
extern void dgelqf(int m, int n, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void dorglq(int m, int n, int k, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void dormlq(char side, char trans, int m, int n, int k, int lda, double a[], double tau[], int ldc, double c[], double work[], int lwork, int *info);
#endif

/*
 * D5-2. QR factorization (complex matrices)
 */
#if defined(_VLARRAY)
extern void zgeqp3(int m, int n, int lda, doublecomplex a[][lda], int jpvt[], doublecomplex tau[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgeqrf(int m, int n, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zungqr(int m, int n, int k, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zunmqr(char side, char trans, int m, int n, int k, int lda, doublecomplex a[][lda], doublecomplex tau[], int ldc, doublecomplex c[][ldc], doublecomplex work[], int lwork, int *info);
extern void zgelqf(int m, int n, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zunglq(int m, int n, int k, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zunmlq(char side, char trans, int m, int n, int k, int lda, doublecomplex a[][lda], doublecomplex tau[], int ldc, doublecomplex c[][ldc], doublecomplex work[], int lwork, int *info);
#else
extern void zgeqp3(int m, int n, int lda, doublecomplex a[], int jpvt[], doublecomplex tau[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgeqrf(int m, int n, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zungqr(int m, int n, int k, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zunmqr(char side, char trans, int m, int n, int k, int lda, doublecomplex a[], doublecomplex tau[], int ldc, doublecomplex c[], doublecomplex work[], int lwork, int *info);
extern void zgelqf(int m, int n, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zunglq(int m, int n, int k, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void zunmlq(char side, char trans, int m, int n, int k, int lda, doublecomplex a[], doublecomplex tau[], int ldc, doublecomplex c[], doublecomplex work[], int lwork, int *info);
#endif

/*
 * D6. Singular value decomposition (SVD)
 */
#if defined(_VLARRAY)
extern void dgesvd(char jobu, char jobvt, int m, int n, int lda, double a[][lda], double s[], int ldu, double u[][ldu], int ldvt, double vt[][ldvt], double work[], int lwork, int *info);
extern void dgesvdx(char jobu, char jobvt, char range, int m, int n, int lda, double a[][lda], double vl, double vu, int il, int iu, int *ns, double s[], int ldu, double u[][ldu], int ldvt, double vt[][ldvt], double work[], int lwork, int iwork[], int *info);
extern void dgesdd(char jobz, int m, int n, int lda, double a[][lda], double s[], int ldu, double u[][ldu], int ldvt, double vt[][ldvt], double work[], int lwork, int iwork[], int *info);
extern void dgesvdq(char joba, char jobp, char jobr, char jobu, char jobv, int m, int n, int lda, double a[][lda], double s[], int ldu, double u[][ldu], int ldv, double v[][ldv], int *numrank, double work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void dgejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, int m, int n, int lda, double a[][lda], double sva[], int ldu, double u[][ldu], int ldv, double v[][ldv], double work[], int lwork, int iwork[], int *info);
extern void dggsvd3(char jobu, char jobv, char jobq, int m, int n, int p, int *k, int *l, int lda, double a[][lda], int ldb, double b[][ldb], double alpha[], double beta[], int ldu, double u[][ldu], int ldv, double v[][ldv], int ldq, double q[][ldq], double work[], int lwork, int iwork[], int *info);
#else
extern void dgesvd(char jobu, char jobvt, int m, int n, int lda, double a[], double s[], int ldu, double u[], int ldvt, double vt[], double work[], int lwork, int *info);
extern void dgesvdx(char jobu, char jobvt, char range, int m, int n, int lda, double a[], double vl, double vu, int il, int iu, int *ns, double s[], int ldu, double u[], int ldvt, double vt[], double work[], int lwork, int iwork[], int *info);
extern void dgesdd(char jobz, int m, int n, int lda, double a[], double s[], int ldu, double u[], int ldvt, double vt[], double work[], int lwork, int iwork[], int *info);
extern void dgesvdq(char joba, char jobp, char jobr, char jobu, char jobv, int m, int n, int lda, double a[], double s[], int ldu, double u[], int ldv, double v[], int *numrank, double work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void dgejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, int m, int n, int lda, double a[], double sva[], int ldu, double u[], int ldv, double v[], double work[], int lwork, int iwork[], int *info);
extern void dggsvd3(char jobu, char jobv, char jobq, int m, int n, int p, int *k, int *l, int lda, double a[], int ldb, double b[], double alpha[], double beta[], int ldu, double u[], int ldv, double v[], int ldq, double q[], double work[], int lwork, int iwork[], int *info);
#endif

/*
 * D6-2. Singular value decomposition (SVD) (complex matrices)
 */
#if defined(_VLARRAY)
extern void zgesvd(char jobu, char jobvt, int m, int n, int lda, doublecomplex a[][lda], double s[], int ldu, doublecomplex u[][ldu], int ldvt, doublecomplex vt[][ldvt], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgesvdx(char jobu, char jobvt, char range, int m, int n, int lda, doublecomplex a[][lda], double vl, double vu, int il, int iu, int *ns, double s[], int ldu, doublecomplex u[][ldu], int ldvt, doublecomplex vt[][ldvt], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void zgesdd(char jobz, int m, int n, int lda, doublecomplex a[][lda], double s[], int ldu, doublecomplex u[][ldu], int ldvt, doublecomplex vt[][ldvt], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void zgesvdq(char joba, char jobp, char jobr, char jobu, char jobv, int m, int n, int lda, doublecomplex a[][lda], double s[], int ldu, doublecomplex u[][ldu], int ldv, doublecomplex v[][ldv], int *numrank, doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zgejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, int m, int n, int lda, doublecomplex a[][lda], double sva[], int ldu, doublecomplex u[][ldu], int ldv, doublecomplex v[][ldv], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int *info);
extern void zggsvd3(char jobu, char jobv, char jobq, int m, int n, int p, int *k, int *l, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double alpha[], double beta[], int ldu, doublecomplex u[][ldu], int ldv, doublecomplex v[][ldv], int ldq, doublecomplex q[][ldq], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
#else
extern void zgesvd(char jobu, char jobvt, int m, int n, int lda, doublecomplex a[], double s[], int ldu, doublecomplex u[], int ldvt, doublecomplex vt[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgesvdx(char jobu, char jobvt, char range, int m, int n, int lda, doublecomplex a[], double vl, double vu, int il, int iu, int *ns, double s[], int ldu, doublecomplex u[], int ldvt, doublecomplex vt[], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void zgesdd(char jobz, int m, int n, int lda, doublecomplex a[], double s[], int ldu, doublecomplex u[], int ldvt, doublecomplex vt[], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void zgesvdq(char joba, char jobp, char jobr, char jobu, char jobv, int m, int n, int lda, doublecomplex a[], double s[], int ldu, doublecomplex u[], int ldv, doublecomplex v[], int *numrank, doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int liwork, int *info);
extern void zgejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, int m, int n, int lda, doublecomplex a[], double sva[], int ldu, doublecomplex u[], int ldv, doublecomplex v[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int *info);
extern void zggsvd3(char jobu, char jobv, char jobq, int m, int n, int p, int *k, int *l, int lda, doublecomplex a[], int ldb, doublecomplex b[], double alpha[], double beta[], int ldu, doublecomplex u[], int ldv, doublecomplex v[], int ldq, doublecomplex q[], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
#endif

/*
 * D9a. Singular, overdetermined or underdetermined systems of linear equations without constraints
 */
#if defined(_VLARRAY)
extern void dgels(char trans, int m, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], double work[], int lwork, int *info);
extern void dgelsy(int m, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], int jpvt[], double rcond, int *rank, double work[], int lwork, int *info);
extern void dgelss(int m, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], double s[], double rcond, int *rank, double work[], int lwork, int *info);
extern void dgetsls(char trans, int m, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], double work[], int lwork, int *info);
extern void dgelsd(int m, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], double s[], double rcond, int *rank, double work[], int lwork, int iwork[], int *info);
extern void dgecov(int job, int n, int lda, double a[][lda], double ci[], int *info);
extern void dgecovy(int job, int n, int lda, double a[][lda], int ipiv[], double ci[], int iwork[], int *info);
extern void dgecovs(int job, int n, int lda, double a[][lda], double s[], double ci[], double work[], int *info);
#else
extern void dgels(char trans, int m, int n, int nrhs, int lda, double a[], int ldb, double b[], double work[], int lwork, int *info);
extern void dgelsy(int m, int n, int nrhs, int lda, double a[], int ldb, double b[], int jpvt[], double rcond, int *rank, double work[], int lwork, int *info);
extern void dgelss(int m, int n, int nrhs, int lda, double a[], int ldb, double b[], double s[], double rcond, int *rank, double work[], int lwork, int *info);
extern void dgetsls(char trans, int m, int n, int nrhs, int lda, double a[], int ldb, double b[], double work[], int lwork, int *info);
extern void dgelsd(int m, int n, int nrhs, int lda, double a[], int ldb, double b[], double s[], double rcond, int *rank, double work[], int lwork, int iwork[], int *info);
extern void dgecov(int job, int n, int lda, double a[], double ci[], int *info);
extern void dgecovy(int job, int n, int lda, double a[], int ipiv[], double ci[], int iwork[], int *info);
extern void dgecovs(int job, int n, int lda, double a[], double s[], double ci[], double work[], int *info);
#endif

/*
 * D9a-2. Singular, overdetermined or underdetermined systems of linear equations without constraints (complex matrices)
 */
#if defined(_VLARRAY)
extern void zgels(char trans, int m, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex work[], int lwork, int *info);
extern void zgelsy(int m, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int jpvt[], double rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgelss(int m, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double s[], double rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgetsls(char trans, int m, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex work[], int lwork, int *info);
extern void zgelsd(int m, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double s[], double rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void zgecov(int job, int n, int lda, doublecomplex a[][lda], doublecomplex ci[], int *info);
extern void zgecovy(int job, int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex ci[], int iwork[], int *info);
extern void zgecovs(int job, int n, int lda, doublecomplex a[][lda], double s[], doublecomplex ci[], doublecomplex work[], int *info);
#else
extern void zgels(char trans, int m, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex work[], int lwork, int *info);
extern void zgelsy(int m, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], int jpvt[], double rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgelss(int m, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], double s[], double rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int *info);
extern void zgetsls(char trans, int m, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex work[], int lwork, int *info);
extern void zgelsd(int m, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], double s[], double rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void zgecov(int job, int n, int lda, doublecomplex a[], doublecomplex ci[], int *info);
extern void zgecovy(int job, int n, int lda, doublecomplex a[], int ipiv[], doublecomplex ci[], int iwork[], int *info);
extern void zgecovs(int job, int n, int lda, doublecomplex a[], double s[], doublecomplex ci[], doublecomplex work[], int *info);
#endif

/*
 * D9b. Singular, overdetermined or underdetermined systems of linear equations with constraints
 */
#if defined(_VLARRAY)
extern void dgglse(int m, int n, int p, int lda, double a[][lda], int ldb, double b[][ldb], double c[], double d[], double x[], double work[], int lwork, int *info);
extern void dggglm(int n, int m, int p, int lda, double a[][lda], int ldb, double b[][ldb], double d[], double x[], double y[], double work[], int lwork, int *info);
#else
extern void dgglse(int m, int n, int p, int lda, double a[], int ldb, double b[], double c[], double d[], double x[], double work[], int lwork, int *info);
extern void dggglm(int n, int m, int p, int lda, double a[], int ldb, double b[], double d[], double x[], double y[], double work[], int lwork, int *info);
#endif

/*
 * D9b-2. Singular, overdetermined or underdetermined systems of complex linear equations with constraints (complex matrices)
 */
#if defined(_VLARRAY)
extern void zgglse(int m, int n, int p, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex c[], doublecomplex d[], doublecomplex x[], doublecomplex work[], int lwork, int *info);
extern void zggglm(int n, int m, int p, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex d[], doublecomplex x[], doublecomplex y[], doublecomplex work[], int lwork, int *info);
#else
extern void zgglse(int m, int n, int p, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex c[], doublecomplex d[], doublecomplex x[], doublecomplex work[], int lwork, int *info);
extern void zggglm(int n, int m, int p, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex d[], doublecomplex x[], doublecomplex y[], doublecomplex work[], int lwork, int *info);
#endif

/*
 * E. Interpolation
 */
extern void polint(int n, double x[], double y[], double c[], int *info);
extern void polyvl(int nder, double xx, double *yfit, double yp[], int n, double x[], double c[], double work[], int *info);
extern void polcof(double xx, int n, double x[], double c[], double d[], double work[]);
extern double fitlag(double x, int n, double a[], double f[], int m, double eps, int *info);
extern void pchim(int n, double x[], double f[], double d[], int incfd, int *info);
extern void pchic(int ic[], double vc[], double sw, int n, double x[], double f[], double d[], int incfd, double work[], int lwork, int *info);
extern void pchse(int n, double x[], double f[], double d[], int incfd, double work[], int lwork, int *info);
extern void pchsp(int ic[], double vc[], int n, double x[], double f[], double d[], int incfd, double work[], int lwork, int *info);
extern void pchfe(int n, double x[], double f[], double d[], int incfd, int skip, int ne, double xe[], double fe[], int *info);
extern void pchfd(int n, double x[], double f[], double d[], int incfd, int skip, int ne, double xe[], double fe[], double de[], int *info);
extern void chfev(double x1, double x2, double f1, double f2, double d1, double d2, int ne, double xe[], double fe[], int next[], int *info);
extern void chfdv(double x1, double x2, double f1, double f2, double d1, double d2, int ne, double xe[], double fe[], double de[], int next[], int *info);
extern void pchbs(int n, double x[], double f[], double d[], int incfd, int knotyp, int *nknots, double t[], double bcoef[], int *ndim, int *kord, int *info);
extern void pchcm(int n, double x[], double f[], double d[], int incfd, int skip, int ismon[], int *info);

extern void bint4(double x[], double y[], int ndata, int ibcl, int ibcr, double fbcl, double fbcr, int kntopt, double t[], double bcoef[], int *n, int *k, double work[], int *info);
extern void bintk(double x[], double y[], double t[], int n, int k, double bcoef[], double q[], double work[], int *info);
extern double bvalue(double t[], double a[], int n, int k, int ideriv, double x, int *inbv, double work[], int *info);
#if defined(_VLARRAY)
extern double ppvalu(int ldc, double c[][ldc], double xi[], int lxi, int k, int ideriv, double x, int *inppv, int *info);
extern void bsplpp(double t[], double a[], int n, int k, int ldc, double c[][ldc], double xi[], int *lxi, double work[], int *info);
extern void bsplvd(double t[], int k, int nderiv, double x, int ileft, int ldvnikx, double vnikx[][ldvnikx], double work[], int *info);
extern void banfac(int n, int kl, int ku, int ldab, double ab[][ldab], int *info);
extern void banslv(int n, int kl, int ku, int ldab, double ab[][ldab], double b[], int *info);
#else
extern double ppvalu(int ldc, double c[], double xi[], int lxi, int k, int ideriv, double x, int *inppv, int *info);
extern void bsplpp(double t[], double a[], int n, int k, int ldc, double c[], double xi[], int *lxi, double work[], int *info);
extern void bsplvd(double t[], int k, int nderiv, double x, int ileft, int ldvnikx, double vnikx[], double work[], int *info);
extern void banfac(int n, int kl, int ku, int ldab, double ab[], int *info);
extern void banslv(int n, int kl, int ku, int ldab, double ab[], double b[], int *info);
#endif
extern void bsplvn(double t[], int jhigh, int k, int index, double x, int ileft, double vnikx[], double work[], int *iwork, int *info);
extern void bspldr(double t[], double a[], int n, int k, int nderiv, double ad[], int *info);
extern void bsplev(double t[], double ad[], int n, int k, int nderiv, double x, int *inev, double svalue[], double work[], int *info);
extern void interv(double xt[], int lxt, double x, int *ilo, int *ileft, int *info);

extern double pchia(int n, double x[], double f[], double d[], int incfd, int skip, double a, double b, int *info);
extern double pchid(int n, double x[], double f[], double d[], int incfd, int skip, int ia, int ib, int *info);
extern void bsqad(double t[], double bcoef[], int n, int k, double x1, double x2, double *bquad, double work[], int *info);
extern void bfqad(double (*f)(double), double t[], double bcoef[], int n, int k, int id, double x1, double x2, double tol, double *quad, double work[], int *info);
extern void bfqad_r(double t[], double bcoef[], int n, int k, int id, double x1, double x2, double tol, double *quad, double work[], int *info, double *xx, double yy, int *irev);
#if defined(_VLARRAY)
extern void ppqad(int ldc, double c[][ldc], double xi[], int lxi, int k, double x1, double x2, double *pquad, int *info);
extern void pfqad(double (*f)(double), int ldc, double c[][ldc], double xi[], int lxi, int k, int id, double x1, double x2, double tol, double *quad, int *info);
extern void pfqad_r(int ldc, double c[][ldc], double xi[], int lxi, int k, int id, double x1, double x2, double tol, double *quad, int *info, double *xx, double yy, int *irev);
#else
extern void ppqad(int ldc, double c[], double xi[], int lxi, int k, double x1, double x2, double *pquad, int *info);
extern void pfqad(double (*f)(double), int ldc, double c[], double xi[], int lxi, int k, int id, double x1, double x2, double tol, double *quad, int *info);
extern void pfqad_r(int ldc, double c[], double xi[], int lxi, int k, int id, double x1, double x2, double tol, double *quad, int *info, double *xx, double yy, int *irev);
#endif

/*
 * F1a. Roots of polynomials
 */
extern void cpzero(int n, doublecomplex a[], doublecomplex r[], int iflag, int maxiter, int *iter, double s[], doublecomplex work[], int *info);
extern void rpzero(int n, double a[], doublecomplex r[], int iflag, int maxiter, int *iter, double s[], doublecomplex work[], int *info);
extern void rpzero2(int n, double a[], double rr[], double ri[], int iflag, int maxiter, int *iter, double s[], double work[], int *info);
extern void cpqr79(int n, doublecomplex a[], doublecomplex r[], doublecomplex work[], int lwork, int *info);
extern void rpqr79(int n, double a[], doublecomplex r[], double work[], int lwork, int *info);
extern void dka(int n, doublecomplex a[], doublecomplex r[], int maxiter, int *iter, doublecomplex work[], double rwork[], int *info);

/*
 * F1b. Solution of single general nonlinear equation
 */
extern void dfzero(double (*f)(double), double *b, double *c, double r, double re, double ae, int *info);
extern void dfzero_r(double *b, double *c, double r, double re, double ae, int *info, double *xx, double yy, int *irev);

/*
 * F2. Solution of a system of nonlinear equations
 */
#if defined(_VLARRAY)
extern void hybrj(void (*fcn)(int, double *, double *, int, double (*)[*], int *), int n, double x[], double fvec[], double xtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, double work[], int lwork, int *info);
extern void hybrj1(void (*fcn)(int, double *, double *, int, double (*)[*], int *), int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info);
#else
extern void hybrj(void (*fcn)(int, double *, double *, int, double *, int *), int n, double x[], double fvec[], double xtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, double work[], int lwork, int *info);
extern void hybrj1(void (*fcn)(int, double *, double *, int, double *, int *), int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info);
#endif
extern void hybrd(void (*fcn)(int, double *, double *, int *), int n, double x[], double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, double work[], int lwork, int *info);
extern void hybrd1(void (*fcn)(int, double *, double *, int *), int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info);
#if defined(_VLARRAY)
extern void chkder(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double xp[], double fvecp[], int mode, double err[], int *info);
#else
extern void chkder(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double xp[], double fvecp[], int mode, double err[], int *info);
#endif
extern void sos(void (*fnc)(int, double *, int, double *), int neq, double x[], double rtolx, double atolx, double tolf, int nprint, int maxiter, int *iter, double work[], int lwork, int iwork[], int liwork, int *info);

#if defined(_VLARRAY)
extern void hybrj_r(int n, double x[], double fvec[], double xtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, double work[], int lwork, int *info, double xx[], double yy[], int ldyypd, double yypd[][ldyypd], int *irev);
extern void hybrj1_r(int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info, double xx[], double yy[], int ldyypd, double yypd[][ldyypd], int *irev);
#else
extern void hybrj_r(int n, double x[], double fvec[], double xtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, double work[], int lwork, int *info, double xx[], double yy[], int ldyypd, double yypd[], int *irev);
extern void hybrj1_r(int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info, double xx[], double yy[], int ldyypd, double yypd[], int *irev);
#endif
extern void hybrd_r(int n, double x[], double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void hybrd1_r(int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void sos_r(int neq, double x[], double rtolx, double atolx, double tolf, int nprint, int maxiter, int *iter, double work[], int lwork, int iwork[], int liwork, int *info, int *k, double xx[], double *yy, int *irev);

/*
 * G1a. Unconstrained optimization of a general univariate function
 */
extern double dfmin(double ax, double bx, double (*f)(double), double tol);
extern void dfmin_r(double ax, double bx, double tol, double *xx, double yy, int *irev);

/*
 * G1b. Unconstrained optimization of a general multivariate function
 */
#if defined(_VLARRAY)
extern void optif9(int n, double x[], void (*fcn)(int, double *, double *), void (*d1fcn)(int, double *, double *), void (*d2fcn)(int, double *, int, double [][*]), double typsiz[], double fscale, int method, int iexp, int ndigit, int maxiter, int iagflg, int iahflg, double dlt, double gradtl, double stepmx, double steptl, void (*result)(int, double *, double, double *, int, double [][*], double *, int, int), int iresult, double xpls[], double *fpls, double gpls[], int *iter, double work[], int lwork, int *info);
#else
extern void optif9(int n, double x[], void (*fcn)(int, double *, double *), void (*d1fcn)(int, double *, double *), void (*d2fcn)(int, double *, int, double *), double typsiz[], double fscale, int method, int iexp, int ndigit, int maxiter, int iagflg, int iahflg, double dlt, double gradtl, double stepmx, double steptl, void (*result)(int, double *, double, double *, int, double *, double *, int, int), int iresult, double xpls[], double *fpls, double gpls[], int *iter, double work[], int lwork, int *info);
#endif
extern void optif0(int n, double x[], void (*fcn)(int, double *, double *), double xpls[], double *fpls, double work[], int lwork, int *info);
extern void mng(int n, double x[], void (*calcf)(int, double *, int *, double *), void (*calcg)(int, double *, int *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);
extern void mnf(int n, double x[], void (*calcf)(int, double *, int *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);
extern void mnh(int n, double x[], void (*calcf)(int, double *, int *, double *), void (*calcgh)(int, double *, int *, double *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);
extern void ivset(int alg, double v[], int lv, int iv[], int liv);
extern void subplex(int n, double x[], void (*f)(int, double *, double *), double tol, int maxeval, double scale[], double *fx, int *neval, int nsmin, int nsmax, int icont, int nfstop, double fstop, int minf, double alpha, double beta, double gamma, double delta, double psi, double omega, int irepl, int ifxsw, double bonus, double work[], int lwork, int iwork[], int liwork, int *info);

#if defined(_VLARRAY)
extern void optif9_r(int n, double x[], double typsiz[], double fscale, int method, int iexp, int ndigit, int maxiter, int iagflg, int iahflg, double dlt, double gradtl, double stepmx, double steptl, int iresult, double xpls[], double *fpls, double gpls[], int *iter, double work[], int lwork, int *info, double xx[], double *yy, double yyp[], int ldyyp2, double yyp2[][ldyyp2], int *irev);
#else
extern void optif9_r(int n, double x[], double typsiz[], double fscale, int method, int iexp, int ndigit, int maxiter, int iagflg, int iahflg, double dlt, double gradtl, double stepmx, double steptl, int iresult, double xpls[], double *fpls, double gpls[], int *iter, double work[], int lwork, int *info, double xx[], double *yy, double yyp[], int ldyyp2, double yyp2[], int *irev);
#endif
extern void optif0_r(int n, double x[], double xpls[], double *fpls, double work[], int lwork, int *info, double xx[], double yy, int *irev);
extern void mng_r(int n, double x[], double d[], double v[], int lv, int iv[], int liv, int *info, double *yy, double yyp[], int *irev);
extern void mnf_r(int n, double x[], double d[], double v[], int lv, int iv[], int liv, int *info, double *yy, int *irev);
extern void mnh_r(int n, double x[], double d[], double v[], int lv, int iv[], int liv, int *info, double *yy, double yyp[], double yypd[], int *irev);
extern void subplex_r(int n, double x[], double tol, int maxeval, double scale[], double *fx, int *neval, int nsmin, int nsmax, int icont, int nfstop, double fstop, int minf, double alpha, double beta, double gamma, double delta, double psi, double omega, int irepl, int ifxsw, double bonus, double work[], int lwork, int iwork[], int liwork, int *info, double yy, int *irev);

/*
 * G2. Constrained optimization
 */
extern void mngb(int n, double x[], double b[][2], void (*calcf)(int, double *, int *, double *), void (*calcg)(int, double *, int *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);
extern void mnfb(int n, double x[], double b[][2], void (*calcf)(int, double *, int *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);
extern void mnhb(int n, double x[], double b[][2], void (*calcf)(int, double *, int *, double *), void (*calcgh)(int, double *, int *, double *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);

extern void mngb_r(int n, double x[], double b[][2], double d[], double v[], int lv, int iv[], int liv, int *info, double *yy, double yyp[], int *irev);
extern void mnfb_r(int n, double x[], double b[][2], double d[], double v[], int lv, int iv[], int liv, int *info, double *yy, int *irev);
extern void mnhb_r(int n, double x[], double b[][2], double d[], double v[], int lv, int iv[], int liv, int *info, double *yy, double yyp[], double yypd[], int *irev);

/*
 * H2a1a. 1-D finite interval quadrature, integrand available via user-defined procedure
 */
extern void qk15(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void qk21(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void qk31(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void qk41(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void qk51(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void qk61(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void qng(double (*f)(double), double a, double b, double epsabs, double epsrel, double *result, double *abserr, int *neval, int *info);
extern void qag(double (*f)(double), double a, double b, double epsabs, double epsrel, int key, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info);
extern void qags(double (*f)(double), double a, double b, double epsabs, double epsrel, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info);
extern void defint(double (*f)(double), double a, double b, double eps, int l, double *result, int *neval, int *info);

extern void qk15_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void qk21_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void qk31_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void qk41_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void qk51_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void qk61_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void qng_r(double a, double b, double epsabs, double epsrel, double *result, double *abserr, int *neval, int *info, double *xx, double yy, int *irev);
extern void qag_r(double a, double b, double epsabs, double epsrel, int key, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info, double *xx, double yy, int *irev);
extern void qags_r(double a, double b, double epsabs, double epsrel, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info, double *xx, double yy, int *irev);
extern void defint_r(double a, double b, double eps, int l, double *result, int *neval, int *info, double *xx, double yy, int *irev);

/*
 * H2a1b. 1-D finite interval quadrature, integrand available on a grid
 */
extern void avint(int n, double x[], double y[], double a, double b, double *result, int *info);

/*
 * H2a2a. 1-D finite interval quadrature, special integrand, integrand available via user-defined procedure
 */
extern void qagp(double (*f)(double), double a, double b, int npts, double points[], double epsabs, double epsrel, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int liwork, int *info);
extern void qawc(double (*f)(double), double a, double b, double c, double epsabs, double epsrel, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info);
extern void qaws(double (*f)(double), double a, double b, double alfa, double beta, int integr, double epsabs, double epsrel, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info);
extern void qawo(double (*f)(double), double a, double b, double omega, int integr, double epsabs, double epsrel, int limit, int maxp1, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int liwork, int *info);

extern void qagp_r(double a, double b, int npts, double points[], double epsabs, double epsrel, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int liwork, int *info, double *xx, double yy, int *irev);
extern void qawc_r(double a, double b, double c, double epsabs, double epsrel, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info, double *xx, double yy, int *irev);
extern void qaws_r(double a, double b, double alfa, double beta, int integr, double epsabs, double epsrel, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info, double *xx, double yy, int *irev);
extern void qawo_r(double a, double b, double omega, int integr, double epsabs, double epsrel, int limit, int maxp1, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int liwork, int *info, double *xx, double yy, int *irev);

/*
 * H2a3a. 1-D semi-infinite interval quadrature, integrand available via
user-defined procedure
 */
extern void qk15i(double (*f)(double), double bound, int inf, double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void qagi(double (*f)(double), double bound, int inf, double epsabs, double epsrel, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info);
extern void qawf(double (*f)(double), double a, double omega, int integr, double epsabs, int limlst, int limit, int maxp1, double *result, double *abserr, int *neval, int *lst, double work[], int lwork, int iwork[], int liwork, int *info);
extern void dehint(double (*f)(double), double a, double eps, double *result, int *neval, int *l, int *info);
extern void deoint(double (*f)(double), double a, double omega, int iw, double eps, double *result, int *neval, double *err, int *info);

extern void qk15i_r(double bound, int inf, double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void qagi_r(double bound, int inf, double epsabs, double epsrel, int limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info, double *xx, double yy, int *irev);
extern void qawf_r(double a, double omega, int integr, double epsabs, int limlst, int limit, int maxp1, double *result, double *abserr, int *neval, int *lst, double work[], int lwork, int iwork[], int liwork, int *info, double *xx, double yy, int *irev);
extern void dehint_r(double a, double eps, double *result, int *neval, int *l, int *info, double *xx, double yy, int *irev);
extern void deoint_r(double a, double omega, int iw, double eps, double *result, int *neval, double *err, int *info, double *xx, double yy, int *irev);

/*
 * H2a4. 1-D infinite interval quadrature
 */
extern void deiint(double (*f)(double), double eps, double *result, int *neval, int *info);
extern void deiint_r(double eps, double *result, int *neval, int *info, double *xx, double yy, int *irev);

/*
 * I1a1. Nonstiff initial value problem for ordinary differential equations (ODEs)
 */
extern void derkf(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void derkf_int(int n, double t, double y[], double work[]);
extern void dopri5(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info);
extern double contd5(int ii, double x, double rcont[], int icont[]);
extern void dverk(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double tol, double c[], int lc, double work[], int lwork, int *info);
extern void dverk_int(int n, double t, double y[], double c[], double work[]);
extern void dop853(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info);
extern double contd8(int ii, double x, double rcont[], int icont[]);
extern void deabm(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, int itstop, double work[], int lwork, int iwork[], int liwork, int *info);
extern void odex(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info);
extern double contx1(int ii, double x, double rcont[], int icont[]);
extern void doprin(int n, void (*f)(int, double, double *, double *), double *t, double y[], double yp[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double *, double *, int, int *), int iout, double work[], int lwork, int iwork[], int liwork, int *info);
extern void odex2(int n, void (*f)(int, double, double *, double *), double *t, double y[], double yp[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double, double *, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info);
extern double contx2(int ii, double x, double rcont[], int icont[]);
extern void retard(int n, void (*f)(int, double, double *, double *, double *, int *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info);
extern double ylag(int ii, double x, double (*phi)(int, double), double rcont[], int icont[]);

extern void derkf_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *xx, double yy[], double yyp[], int *irev);
extern void dopri5_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int *irtrn, int *irev);
extern void dverk_r(int n, double *t, double y[], double tout, double tol, double c[], int lc, double work[], int lwork, int *info, double *tt, double yy[], double yyp[], int *irev);
extern void dop853_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int *irtrn, int *irev);
extern void deabm_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, int itstop, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yyp[], int *irev);
extern void odex_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int *irtrn, int *irev);
extern void doprin_r(int n, double *t, double y[], double yp[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yypp[], int *irtrn, int *irev);
extern void odex2_r(int n, double *t, double y[], double yp[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int *irtrn, int *irev);
extern void retard_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int *irtrn, int *irev);
extern void ylag_r(int ii, double x, double rcont[], int icont[], double *yy, int *irev);

/*
 * I1a2. Stiff initial value problem for ordinary differential equations (ODEs)
 */
extern void debdf(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, int itstop, void (*djac)(int, double, double *, int, double *), int idjac, double work[], int lwork, int iwork[], int liwork, int *info);
extern void radau5(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*jac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*mas)(int, int, double *), int imas, int mlmas, int mumas, void (*solout)(int, double, double, double *, int, double *, int *), int iout, double work[], int lwork, int iwork[], int liwork, double cont[], int lcont, int *info);
extern double contr5(int ii, double x, double cont[]);
extern void radaup(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*jac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*mas)(int, int, double *), int imas, int mlmas, int mumas, void (*solout)(int, double, double, double *, int, double *, int *), int iout, double work[], int lwork, int iwork[], int liwork, double cont[], int lcont, int *info);
extern double contrp(int ii, double x, double cont[]);
extern void radau(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*jac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*mas)(int, int, double *), int imas, int mlmas, int mumas, void (*solout)(int, double, double, double *, int, double *, int *), int iout, double work[], int lwork, int iwork[], int liwork, double cont[], int lcont, int *info);
extern double contra(int ii, double x, double cont[]);
extern void rodas(int n, void (*f)(int, double, double *, double *), int ifcn, double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*jac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*dfx)(int, double, double *, double *), int idfx, void (*mas)(int, int, double *), int imas, int mlmas, int mumas, void (*solout)(int, double, double, double *, int, double *, int *), int iout, double work[], int lwork, int iwork[], int liwork, double cont[], int lcont, int *info);
extern double contro(int ii, double x, double cont[]);
extern void seulex(int n, void (*f)(int, double, double *, double *), int ifcn, double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*jac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*mas)(int, int, double *), int imas, int mlmas, int mumas, void (*solout)(int, double, double, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info);
extern double contex(int ii, double x, double rcont[], int icont[]);
extern void dassl(int n, void (*res)(int, double, double *, double *, double *, int *), double *t, double y[], double yp[], double tout, int iopt[], double *rtol, double *atol, void (*jac)(int, double, double *, double *, int, double *, double), double work[], int lwork, int iwork[], int liwork, int *info);

extern void debdf_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, int itstop, int idjac, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irev);
extern void radau5_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int imas, int mlmas, int mumas, int iout, double work[], int lwork, int iwork[], int liwork, double cont[], int lcont, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irtrn, int *irev);
extern void radaup_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int imas, int mlmas, int mumas, int iout, double work[], int lwork, int iwork[], int liwork, double cont[], int lcont, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irtrn, int *irev);
extern void radau_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int imas, int mlmas, int mumas, int iout, double work[], int lwork, int iwork[], int liwork, double cont[], int lcont, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irtrn, int *irev);
extern void rodas_r(int n, int ifcn, double *t, double y[], double tout, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int idfx, int imas, int mlmas, int mumas, int iout, double work[], int lwork, int iwork[], int liwork, double cont[], int lcont, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irtrn, int *irev);
extern void seulex_r(int n, int ifcn, double *t, double y[], double tout, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int imas, int mlmas, int mumas, int iout, double work[], int lwork, int iwork[], int liwork, double rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irtrn, int *irev);
extern void dassl_r(int n, double *t, double y[], double yp[], double tout, int iopt[], double *rtol, double *atol, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yyp[], int ldyypd, double yypd[], double *cj, int *ires, int *irev);

/*
 * J1a1. One-dimensional real fast Fourier transforms
 */
extern void rfft1f(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void rfft1b(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void rfft1i(int n, double wsave[], int lwsave, int *info);
extern void rfftmf(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void rfftmb(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void rfftmi(int n, double wsave[], int lwsave, int *info);

/*
 * J1a2. One-dimensional complex fast Fourier transforms
 */
extern void cfft1f(int n, int inc, doublecomplex c[], int lc, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cfft1b(int n, int inc, doublecomplex c[], int lc, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cfft1i(int n, double wsave[], int lwsave, int *info);
extern void cfftmf(int lot, int jump, int n, int inc, doublecomplex c[], int lc, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cfftmb(int lot, int jump, int n, int inc, doublecomplex c[], int lc, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cfftmi(int n, double wsave[], int lwsave, int *info);

/*
 * J1a3. One-dimensional trigonometric fast Fourier transforms
 */
extern void sint1f(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void sint1b(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void sint1i(int n, double wsave[], int lwsave, int *info);
extern void sintmf(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void sintmb(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void sintmi(int n, double wsave[], int lwsave, int *info);
extern void cost1f(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cost1b(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cost1i(int n, double wsave[], int lwsave, int *info);
extern void costmf(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void costmb(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void costmi(int n, double wsave[], int lwsave, int *info);

extern void sinq1f(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void sinq1b(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void sinq1i(int n, double wsave[], int lwsave, int *info);
extern void sinqmf(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void sinqmb(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void sinqmi(int n, double wsave[], int lwsave, int *info);
extern void cosq1f(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cosq1b(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cosq1i(int n, double wsave[], int lwsave, int *info);
extern void cosqmf(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cosqmb(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cosqmi(int n, double wsave[], int lwsave, int *info);

/*
 * J1b. Multidimensional fast Fourier transforms
 */
#if defined(_VLARRAY)
extern void rfft2f(int l, int m, int ldr, double r[][ldr], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void rfft2b(int l, int m, int ldr, double r[][ldr], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void rfft2c(int l, int m, int ldr, double r[][ldr], int ldc, doublecomplex c[][ldc], int *info);
#else
extern void rfft2f(int l, int m, int ldr, double r[], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void rfft2b(int l, int m, int ldr, double r[], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void rfft2c(int l, int m, int ldr, double r[], int ldc, doublecomplex c[], int *info);
#endif
extern void rfft2i(int l, int m, double wsave[], int lwsave, int *info);
#if defined(_VLARRAY)
extern void cfft2f(int l, int m, int ldc, doublecomplex c[][ldc], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cfft2b(int l, int m, int ldc, doublecomplex c[][ldc], double wsave[], int lwsave, double work[], int lwork, int *info);
#else
extern void cfft2f(int l, int m, int ldc, doublecomplex c[], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void cfft2b(int l, int m, int ldc, doublecomplex c[], double wsave[], int lwsave, double work[], int lwork, int *info);
#endif
extern void cfft2i(int l, int m, double wsave[], int lwsave, int *info);

/*
 * K1b1. Unconstrained nonlinear least squares approximation
 */
#if defined(_VLARRAY)
extern void lmder(void (*fcn)(int, int, double *, double *, int, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info);
extern void lmder1(void (*fcn)(int, int, double *, double *, int, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double tol, int ipvt[], double work[], int lwork, int *info);
extern void lmstr(void (*fcn)(int, int, double *, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info);
extern void lmstr1(void (*fcn)(int, int, double *, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double tol, int ipvt[], double work[], int lwork, int *info);
extern void lmdif(void (*fcn)(int, int, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, int ipvt[], double work[], int lwork, int *info);
extern void lmdif1(void (*fcn)(int, int, double *, double *, int *), int m, int n, double x[], double fvec[], double tol, double work[], int lwork, int iwork[], int *info);
extern void covar(int n, int ldr, double r[][ldr], int ipvt[], double tol, double work[], int *info);
#else
extern void lmder(void (*fcn)(int, int, double *, double *, int, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info);
extern void lmder1(void (*fcn)(int, int, double *, double *, int, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double tol, int ipvt[], double work[], int lwork, int *info);
extern void lmstr(void (*fcn)(int, int, double *, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info);
extern void lmstr1(void (*fcn)(int, int, double *, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double tol, int ipvt[], double work[], int lwork, int *info);
extern void lmdif(void (*fcn)(int, int, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, int ipvt[], double work[], int lwork, int *info);
extern void lmdif1(void (*fcn)(int, int, double *, double *, int *), int m, int n, double x[], double fvec[], double tol, double work[], int lwork, int iwork[], int *info);
extern void covar(int n, int ldr, double r[], int ipvt[], double tol, double work[], int *info);
#endif
extern void n2g(int n, int p, double x[], void (*calcr)(int, int, double *, int *, double *), void (*calcj)(int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);
extern void n2f(int n, int p, double x[], void (*calcr)(int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);
extern void n2p(int n, int md, int p, double x[], void (*calcr)(int, int, int, int *, int, double *, int *, double *), void (*calcj)(int, int, int, int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);

#if defined(_VLARRAY)
extern void lmder_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void lmder1_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double tol, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void lmstr_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], double yypr[], int *irev);
extern void lmstr1_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double tol, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], double yypr[], int *irev);
extern void lmdif_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
#else
extern void lmder_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void lmder1_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double tol, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void lmstr_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], double yypr[], int *irev);
extern void lmstr1_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double tol, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], double yypr[], int *irev);
extern void lmdif_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
#endif
extern void lmdif1_r(int m, int n, double x[], double fvec[], double tol, double work[], int lwork, int iwork[], int *info, double xx[], double yy[], int *irev);
#if defined(_VLARRAY)
extern void n2g_r(int n, int p, double x[], double v[], int lv, int iv[], int liv, int *info, double yy[], int ldyyp, double yyp[][ldyyp], int *irev);
extern void n2p_r(int n, int md, int p, double x[], double v[], int lv, int iv[], int liv, int *info, int *m1, int *m2, double yy[], int ldyyp, double yyp[][ldyyp], int *irev);
#else
extern void n2g_r(int n, int p, double x[], double v[], int lv, int iv[], int liv, int *info, double yy[], int ldyyp, double yyp[], int *irev);
extern void n2p_r(int n, int md, int p, double x[], double v[], int lv, int iv[], int liv, int *info, int *m1, int *m2, double yy[], int ldyyp, double yyp[], int *irev);
#endif
extern void n2f_r(int n, int p, double x[], double v[], int lv, int iv[], int liv, int *info, double yy[], int *irev);

/*
 * K1b2. Constrained nonlinear least squares approximation
 */
extern void n2gb(int n, int p, double x[], double b[][2], void (*calcr)(int, int, double *, int *, double *), void (*calcj)(int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);
extern void n2fb(int n, int p, double x[], double b[][2], void (*calcr)(int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);
extern void n2pb(int n, int nd, int p, double x[], double b[][2], void (*calcr)(int, int, int, int *, int, double *, int *, double *), void (*calcj)(int, int, int, int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int liv, int *info);

#if defined(_VLARRAY)
extern void n2gb_r(int n, int p, double x[], double b[][2], double v[], int lv, int iv[], int liv, int *info, double yy[], int ldyyp, double yyp[][ldyyp], int *irev);
extern void n2pb_r(int n, int nd, int p, double x[], double b[][2], double v[], int lv, int iv[], int liv, int *info, int *n1, int *n2, double yy[], int ldyyp, double yyp[][ldyyp], int *irev);
#else
extern void n2gb_r(int n, int p, double x[], double b[][2], double v[], int lv, int iv[], int liv, int *info, double yy[], int ldyyp, double yyp[], int *irev);
extern void n2pb_r(int n, int nd, int p, double x[], double b[][2], double v[], int lv, int iv[], int liv, int *info, int *n1, int *n2, double yy[], int ldyyp, double yyp[], int *irev);
#endif
extern void n2fb_r(int n, int p, double x[], double b[][2], double v[], int lv, int iv[], int liv, int *info, double yy[], int *irev);

/*
 * L6a5. Exponential random numbers
 */
extern double zigexp(long (*iurand)(void), double (*durand)(void), double theta);
extern void init_zigexp(int n);
extern void zigexp_r(double theta, double *r, int ii, double yy, int *irev);
extern void init_zigexp_r(int n);

/*
 * L6a7. Gamma random numbers
 */
extern double rgama(double (*durand)(void), double (*dnrand)(void), double alpha, double beta);
extern void rgama_r(double alpha, double beta, double *r, double yy, int *irev);

/*
 * L6a14. Normal random numbers
 */
extern double zignorm(long (*iurand)(void), double (*durand)(void), double mu, double sigma);
extern void init_zignorm(int n);
extern void zignorm_r(double mu, double sigma, double *r, int ii, double yy, int *irev);
extern void init_zignorm_r(int n);

/*
 * L6a21. Uniform random numbers
 */
extern void init_genrand(unsigned long s);
extern void init_by_array(unsigned long init_key[], int key_length);
extern unsigned long genrand_int32(void);
extern long genrand_int31(void);
extern double genrand_real1(void);
extern double genrand_real2(void);
extern double genrand_real3(void);
extern double genrand_res53(void);
extern void init_genrand64(unsigned long long s);
extern void init_by_array64(unsigned long long initkey[], unsigned long long keylength);
extern unsigned long long genrand64_int64(void);
extern long long genrand64_int63(void);
extern double genrand64_real1(void);
extern double genrand64_real2(void);
extern double genrand64_real3(void);

extern void ran_start(long seed);
extern void ran_array(long aa[], int n);
extern long ran_arr_next(void);
extern void ranf_start(long seed);
extern void ranf_array(double aa[], int n);
extern double ranf_arr_next(void);

extern void srand48(long seed);
extern unsigned short *seed48(unsigned short xseed[]);
extern void lcong48(unsigned short p[]);
extern double drand48(void);
extern double erand48(unsigned short xseed[]);
extern long lrand48(void);
extern long nrand48(unsigned short xseed[]);
extern long mrand48(void);
extern long jrand48(unsigned short xseed[]);

/*
 * R. Service routines
 */
extern double dlamch(char cmach);
extern float slamch(char cmach);
extern double d1mach(int i);
extern float r1mach(int i);
extern int i1mach(int i);

/*
 * Z. Others
 */
#if defined(_VLARRAY)
extern void dlatme(int n, char dist, int iseed[], double d[], int mode, double cond, double dmax, char ei[], char rsign, char upper, char sim, double ds[], int modes, double conds, int kl, int ku, double anorm, int lda, double a[][lda], double work[], int *info);
extern void dlatmr(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, char rsign, char grade, double dl[], int model, double condl, double dr[], int moder, double condr, char pivtng, int ipivot[], int kl, int ku, double sparse, double anorm, char pack, int lda, double a[][lda], int iwork[], int *info);
extern void dlatms(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int kl, int ku, char pack, int lda, double a[][lda], double work[], int *info);
extern void dlatmt(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int rank, int kl, int ku, char pack, int lda, double a[][lda], double work[], int *info);
extern void zlatme(int n, char dist, int iseed[], doublecomplex d[], int mode, double cond, doublecomplex dmax, char rsign, char upper, char sim, double ds[], int modes, double conds, int kl, int ku, double anorm, int lda, doublecomplex a[][lda], doublecomplex work[], int *info);
extern void zlatmr(int m, int n, char dist, int iseed[], char sym, doublecomplex d[], int mode, double cond, doublecomplex dmax, char rsign, char grade, doublecomplex dl[], int model, double condl, doublecomplex dr[], int moder, double condr, char pivtng, int ipivot[], int kl, int ku, double sparse, double anorm, char pack, int lda, doublecomplex a[][lda], int iwork[], int *info);
extern void zlatms(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int kl, int ku, char pack, int lda, doublecomplex a[][lda], doublecomplex work[], int *info);
extern void zlatmt(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int rank, int kl, int ku, char pack, int lda, doublecomplex a[][lda], doublecomplex work[], int *info);
#else
extern void dlatme(int n, char dist, int iseed[], double d[], int mode, double cond, double dmax, char ei[], char rsign, char upper, char sim, double ds[], int modes, double conds, int kl, int ku, double anorm, int lda, double a[], double work[], int *info);
extern void dlatmr(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, char rsign, char grade, double dl[], int model, double condl, double dr[], int moder, double condr, char pivtng, int ipivot[], int kl, int ku, double sparse, double anorm, char pack, int lda, double a[], int iwork[], int *info);
extern void dlatms(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int kl, int ku, char pack, int lda, double a[], double work[], int *info);
extern void dlatmt(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int rank, int kl, int ku, char pack, int lda, double a[], double work[], int *info);
extern void zlatme(int n, char dist, int iseed[], doublecomplex d[], int mode, double cond, doublecomplex dmax, char rsign, char upper, char sim, double ds[], int modes, double conds, int kl, int ku, double anorm, int lda, doublecomplex a[], doublecomplex work[], int *info);
extern void zlatmr(int m, int n, char dist, int iseed[], char sym, doublecomplex d[], int mode, double cond, doublecomplex dmax, char rsign, char grade, doublecomplex dl[], int model, double condl, doublecomplex dr[], int moder, double condr, char pivtng, int ipivot[], int kl, int ku, double sparse, double anorm, char pack, int lda, doublecomplex a[], int iwork[], int *info);
extern void zlatms(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int kl, int ku, char pack, int lda, doublecomplex a[], doublecomplex work[], int *info);
extern void zlatmt(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int rank, int kl, int ku, char pack, int lda, doublecomplex a[], doublecomplex work[], int *info);
#endif

#undef _VLARRAY

#if defined(__cplusplus)
}
#endif
