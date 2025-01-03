/****************************************
 *                                      *
 *  XLPack Numerical Library            *
 *  Version 7.0 (January 31, 2023)      *
 *  (C) 2014-2023  K Technologies       *
 *                                      *
 ****************************************/
#pragma once

#include "cnumlib_complex.h"

#if !defined(_CNUMLIB_NO_MANGLING)
#include "cnumlib_mangling.h"
#endif

#if defined(__cplusplus)
extern "C" {
#endif

/*
 * A3. Real arithmetic
 */
extern double _d1num(int i);

/*
 * A4. C utility routines for complex arithmetic
 */
#if !defined(__cplusplus)
extern doublecomplex _cmplx(double r, double i);
extern doublecomplex _cpolar(double rho, double theta);
extern doublecomplex _cminus(doublecomplex a);
extern doublecomplex _cadd(doublecomplex a, doublecomplex b);
extern doublecomplex _cdadd(doublecomplex a, double rb);
extern doublecomplex _csub(doublecomplex a, doublecomplex b);
extern doublecomplex _cdsub(doublecomplex a, double rb);
extern doublecomplex _dcsub(double ra, doublecomplex b);
extern doublecomplex _cmul(doublecomplex a, doublecomplex b);
extern doublecomplex _cdmul(doublecomplex a, double rb);
extern doublecomplex _cdiv(doublecomplex a, doublecomplex b);
extern doublecomplex _cddiv(doublecomplex a, double rb);
extern doublecomplex _dcdiv(double ra, doublecomplex b);
extern doublecomplex _cdpow(doublecomplex a, double rb);
extern doublecomplex _cipow(doublecomplex a, int ib);
#endif
#if !defined(_NO_SUB_DEFS)
extern void _conj_sub(doublecomplex z, doublecomplex *conjz);
extern void _cpolar_sub(double rho, double theta, doublecomplex *z);
extern void _cproj_sub(doublecomplex z, doublecomplex *cprojz);
extern void _cmplx_sub(double r, double i, doublecomplex *z);
extern void _cminus_sub(doublecomplex a, doublecomplex *minus);
extern void _cadd_sub(doublecomplex a, doublecomplex b, doublecomplex *add);
extern void _cdadd_sub(doublecomplex a, double rb, doublecomplex *add);
extern void _csub_sub(doublecomplex a, doublecomplex b, doublecomplex *sub);
extern void _cdsub_sub(doublecomplex a, double rb, doublecomplex *sub);
extern void _dcsub_sub(double ra, doublecomplex b, doublecomplex *sub);
extern void _cmul_sub(doublecomplex a, doublecomplex b, doublecomplex *mul);
extern void _cdmul_sub(doublecomplex a, double rb, doublecomplex *mul);
extern void _cdiv_sub(doublecomplex a, doublecomplex b, doublecomplex *div);
extern void _cddiv_sub(doublecomplex a, double rb, doublecomplex *div);
extern void _dcdiv_sub(double ra, doublecomplex b, doublecomplex *div);
extern void _cpow_sub(doublecomplex a, doublecomplex b, doublecomplex *pow);
extern void _cdpow_sub(doublecomplex a, double rb, doublecomplex *pow);
extern void _cipow_sub(doublecomplex a, int ib, doublecomplex *pow);
#endif

/*
 * C1. Special functions (Integer-valued functions)
 */
extern double _factorial(unsigned int n);

/*
 * C2. Special functions (Powers, roots, reciprocals)
 */
#if !defined(__cplusplus)
extern doublecomplex _ccbrt(doublecomplex z);
#endif
#if !defined(_NO_SUB_DEFS)
extern void _ccbrt_sub(doublecomplex z, doublecomplex *cbrtz);
extern void _csqrt_sub(doublecomplex z, doublecomplex *sqrtz);
#endif

/*
 * C3. Special functions (Polynomials)
 */
extern double _laguerre(unsigned int n, double x);
extern double _alaguerre(unsigned int n, unsigned int m, double x);
extern double _legendre(unsigned int n, double x);
extern double _legendred(unsigned int n, double x);
extern double _alegendre(unsigned int n, unsigned int m, double x);
#if !defined(__cplusplus)
extern doublecomplex _sharmonic(unsigned int l, int m, double theta, double phi);
#endif
#if !defined(_NO_SUB_DEFS)
extern void _sharmonic_sub(unsigned int l, int m, double theta, double phi, doublecomplex *z);
#endif
extern double _sharmonicr(unsigned int l, int m, double theta, double phi);
extern double _sharmonici(unsigned int l, int m, double theta, double phi);
extern double _hermite(unsigned int n, double x);
extern double _chebt(unsigned int n, double x);
extern double _chebtd(unsigned int n, double x);
extern double _chebu(unsigned int n, double x);
extern double _chebs(double c[], size_t n, double x);
extern double _gegenbauer(unsigned int n, double lambda, double x);
extern double _gegenbauerd1(unsigned int n, double lambda, double x);
extern double _gegenbauerd(unsigned int n, double lambda, double x, unsigned int k);
extern double _jacobi(unsigned int n, double alpha, double beta, double x);
extern double _jacobid1(unsigned int n, double alpha, double beta, double x);
extern double _jacobid2(unsigned int n, double alpha, double beta, double x);
extern double _jacobid(unsigned int n, double alpha, double beta, double x, unsigned int k);

/*
 * C4. Special functions (Elementary transcendental functions)
 */
#if !defined(__cplusplus)
extern doublecomplex _cexpm1(doublecomplex z);
extern doublecomplex _clog1p(doublecomplex z);
extern doublecomplex _ccot(doublecomplex z);
#endif
#if !defined(_NO_SUB_DEFS)
extern void _cexpm1_sub(doublecomplex z, doublecomplex *expm1z);
extern void _clog1p_sub(doublecomplex z, doublecomplex *log1pz);
extern void _cexp_sub(doublecomplex z, doublecomplex *expz);
extern void _clog_sub(doublecomplex z, doublecomplex *logz);
extern void _ccos_sub(doublecomplex z, doublecomplex *cosz);
extern void _csin_sub(doublecomplex z, doublecomplex *sinz);
extern void _ctan_sub(doublecomplex z, doublecomplex *tanz);
extern void _ccosh_sub(doublecomplex z, doublecomplex *coshz);
extern void _csinh_sub(doublecomplex z, doublecomplex *sinhz);
extern void _ctanh_sub(doublecomplex z, doublecomplex *tanhz);
extern void _cacos_sub(doublecomplex z, doublecomplex *acosz);
extern void _casin_sub(doublecomplex z, doublecomplex *asinz);
extern void _catan_sub(doublecomplex z, doublecomplex *atanz);
extern void _cacosh_sub(doublecomplex z, doublecomplex *acoshz);
extern void _casinh_sub(doublecomplex z, doublecomplex *asinhz);
extern void _catanh_sub(doublecomplex z, doublecomplex *atanhz);
extern void _ccot_sub(doublecomplex z, doublecomplex *cotz);
#endif
extern double _sqrt1pm1(double x);
extern double _powm1(double x, double y);
extern double _sin_pi(double x);
extern double _cos_pi(double x);

/*
 * C5. Special functions (Exponential and logarithmic integrals)
 */
extern double _li(double x);
extern double _ei(double x);
extern double _e1(double x);
extern double _en(unsigned int n, double x);
extern double _spence(double x);

/*
 * C6. Special functions (Cosine and _sine integrals)
 */
extern double _ci(double x);
extern double _si(double x);
extern double _chi(double x);
extern double _shi(double x);

/*
 * C7a. Special functions (Gamma functions)
 */
extern double _tgamma1pm1(double x);
extern double _lgammas(double x, int *sign);
extern double _rgamma(double x);
extern double _tgammaratio(double a, double b);
extern double _tgammadratio(double a, double delta);
#if !defined(__cplusplus)
extern doublecomplex _cgamma(doublecomplex z);
extern doublecomplex _clgamma(doublecomplex z);
extern doublecomplex _crgamma(doublecomplex z);
#endif
#if !defined(_NO_SUB_DEFS)
extern void _cgamma_sub(doublecomplex z, doublecomplex *zout);
extern void _clgamma_sub(doublecomplex z, doublecomplex *zout);
extern void _crgamma_sub(doublecomplex z, doublecomplex *zout);
#endif
extern double _poch(double a, double x);
extern double _poch1(double a, double x);

/*
 * C7b. Special functions (Beta functions)
 */
extern double _beta(double a, double b);
extern double _lbeta(double a, double b);
#if !defined(__cplusplus)
extern doublecomplex _cbeta(doublecomplex a, doublecomplex b);
extern doublecomplex _clbeta(doublecomplex a, doublecomplex b);
#endif
#if !defined(_NO_SUB_DEFS)
extern void _cbeta_sub(doublecomplex a, doublecomplex b, doublecomplex *zout);
extern void _clbeta_sub(doublecomplex a, doublecomplex b, doublecomplex *zout);
#endif

/*
 * C7c. Special functions (Psi function)
 */
extern double _digamma(double x);
extern double _trigamma(double x);
extern double _polygamma(int n, double x);
#if !defined(__cplusplus)
extern doublecomplex _cdigamma(doublecomplex z);
#endif
#if !defined(_NO_SUB_DEFS)
extern void _cdigamma_sub(doublecomplex z, doublecomplex *zout);
#endif

/*
 * C7e. Special functions (Incomplete Gamma functions)
 */
extern double _gammai(double a, double x);
extern double _gammaic(double a, double x);
extern double _gammait(double a, double x);
extern double _gammap(double a, double x);
extern double _gammaq(double a, double x);
extern double _gammapi(double a, double p);
extern double _gammaqi(double a, double q);
extern double _gammapia(double x, double p);
extern double _gammaqia(double x, double q);
extern double _gammapd(double a, double x);

/*
 * C7f. Special functions (Incomplete Beta function)
 */
extern double _betax(double a, double b, double x);
extern double _betaxc(double a, double b, double x);
extern double _ibeta(double a, double b, double x);
extern double _ibetac(double a, double b, double x);
extern double _ibetai(double a, double b, double p, double *py);
extern double _ibetaci(double a, double b, double q, double *py);
extern double _ibetaia(double b, double x, double p);
extern double _ibetacia(double b, double x, double q);
extern double _ibetaib(double a, double x, double p);
extern double _ibetacib(double a, double x, double q);
extern double _ibetad(double a, double b, double x);

/*
 * C7g. Special functions (Riemann _zeta function)
 */
extern double _zeta(double x);

/*
 * C8. Special functions (Error functions)
 */
extern double _erfi(double p);
extern double _erfci(double q);
extern double _dawson(double x);
extern double _fresc(double x);
extern double _fress(double x);

/*
 * C10a. Special functions (Bessel functions)
 */
extern double _besj0(double x);
extern double _besj1(double x);
extern double _besjn(int n, double x);
extern double _besy0(double x);
extern double _besy1(double x);
extern double _besyn(int n, double x);
extern double _besjnu(double nu, double x);
extern double _besynu(double nu, double x);
extern double _besjnd(int n, double x);
extern double _besynd(int n, double x);
extern double _besjnud(double nu, double x);
extern double _besynud(double nu, double x);
extern double _sbesjn(int n, double x);
extern double _sbesyn(int n, double x);
extern double _sbesjnu(double nu, double x);
extern double _sbesynu(double nu, double x);
extern void _cbesh(doublecomplex z, double nu, int kode, int m, int n, doublecomplex y[], int *info);
extern void _cbesj(doublecomplex z, double nu, int kode, int n, doublecomplex y[], int *info);
extern void _cbesy(doublecomplex z, double nu, int kode, int n, doublecomplex y[], doublecomplex work[], int *info);

/*
 * C10b. Special functions (Modified Bessel functions)
 */
extern double _besi0(double x);
extern double _besi0e(double x);
extern double _besi1(double x);
extern double _besi1e(double x);
extern double _besin(int n, double x);
extern double _besk0(double x);
extern double _besk0e(double x);
extern double _besk1(double x);
extern double _besk1e(double x);
extern double _beskn(int n, double x);
extern double _besinu(double nu, double x);
extern double _besknu(double nu, double x);
extern double _besind(int n, double x);
extern double _besknd(int n, double x);
extern double _besinud(double nu, double x);
extern double _besknud(double nu, double x);
extern double _sbesin(int n, double x);
extern double _sbeskn(int n, double x);
extern double _sbesinu(double nu, double x);
extern double _sbesknu(double nu, double x);
extern void _cbesi(doublecomplex z, double nu, int kode, int n, doublecomplex y[], int *info);
extern void _cbesk(doublecomplex z, double nu, int kode, int n, doublecomplex y[], int *info);

/*
 * C10d. Special functions (Airy functions)
 */
extern double _airyai(double x);
extern double _airybi(double x);
extern double _airyaid(double x);
extern double _airybid(double x);
extern void _cairy(doublecomplex z, int id, int kode, doublecomplex *ai, int *info);
extern void _cbiry(doublecomplex z, int id, int kode, doublecomplex *bi, int *info);

/*
 * C11. Special functions (Confluent hypergeometric function)
 */
extern double _chu(double a, double b, double x);
extern double _hyp1f1(double a, double b, double x);
extern double _hyp1f1r(double a, double b, double x);
extern double _lhyp1f1(double a, double b, double z, int *sign);
extern double _hyp2f1(double a, double b, double c, double x);
extern double _hyp0f1(double b, double z);
extern double _hyp1f0(double a, double z);
extern double _hyp2f0(double a1, double a2, double z);
extern double _hyppfq(unsigned int p, unsigned int q, double a[], double b[], double z, double *abserr);

/*
 * C13. Special functions (Jacobi elliptic functions)
 */
extern void _jelli(double u, double k, double *sn, double *cn, double *dn);
extern double _jsn(double u, double k);
extern double _jcn(double u, double k);
extern double _jdn(double u, double k);
extern double _jns(double u, double k);
extern double _jnc(double u, double k);
extern double _jnd(double u, double k);
extern double _jsc(double u, double k);
extern double _jsd(double u, double k);
extern double _jdc(double u, double k);
extern double _jds(double u, double k);
extern double _jcs(double u, double k);
extern double _jcd(double u, double k);
extern double _jtheta1(double x, double q);
extern double _jtheta1t(double x, double tau);
extern double _jtheta2(double x, double q);
extern double _jtheta2t(double x, double tau);
extern double _jtheta3(double x, double q);
extern double _jtheta3t(double x, double tau);
extern double _jtheta3m1(double x, double q);
extern double _jtheta3m1t(double x, double tau);
extern double _jtheta4(double x, double q);
extern double _jtheta4t(double x, double tau);
extern double _jtheta4m1(double x, double q);
extern double _jtheta4m1t(double x, double tau);

/*
 * C14. Special functions (Elliptic integrals)
 */
extern double _celli1(double k);
extern double _celli2(double k);
extern double _celli3(double n, double k);
extern double _elli1(double phi, double k);
extern double _elli2(double phi, double k);
extern double _elli3(double phi, double n, double k);
extern double _rc(double x, double y);
extern double _rd(double x, double y, double z);
extern double _rg(double x, double y, double z);
extern double _rf(double x, double y, double z);
extern double _rj(double x, double y, double z, double p);
extern double _jzeta(double phi, double k);
extern double _hlambda(double phi, double k);

/*
 * C19. Special functions (other special functions)
 */
extern double _dconst(int i);

/*
 * D1a. Elementary vector operations (BLAS 1)
 */
extern void _daxpy(int n, double a, double x[], int incx, double y[], int incy);
extern void _dcopy(int n, double x[], int incx, double y[], int incy);
extern double _ddot(int n, double x[], int incx, double y[], int incy);
extern void _drotg(double *a, double *b, double *c, double *s);
extern void _drotmg(double *d1, double *d2, double *x1, double y1, double p[]);
extern void _drot(int n, double x[], int incx, double y[], int incy, double c, double s);
extern void _drotm(int n, double x[], int incx, double y[], int incy, double p[]);
extern void _dscal(int n, double a, double x[], int incx);
extern void _dswap(int n, double x[], int incx, double y[], int incy);
extern double _dasum(int n, double x[], int incx);
extern double _dnrm2(int n, double x[], int incx);
extern int _idamax(int n, double x[], int incx);
extern int _lsame(char ca, char cb);

/*
 * D1a-2. Elementary vector operations (BLAS 1) (Complex)
 */
extern void _zaxpy(int n, doublecomplex a, doublecomplex x[], int incx, doublecomplex y[], int incy);
extern void _zcopy(int n, doublecomplex x[], int incx, doublecomplex y[], int incy);
#if !defined(__cplusplus)
extern doublecomplex _zdotu(int n, doublecomplex x[], int incx, doublecomplex y[], int incy);
extern doublecomplex _zdotc(int n, doublecomplex x[], int incx, doublecomplex y[], int incy);
#endif
#if !defined(_NO_SUB_DEFS)
extern void _zdotu_sub(doublecomplex *dotu, int n, doublecomplex zx[], int incx, doublecomplex zy[], int incy);
extern void _zdotc_sub(doublecomplex *dotc, int n, doublecomplex zx[], int incx, doublecomplex zy[], int incy);
#endif
extern void _zrotg(doublecomplex *a, doublecomplex *b, double *c, doublecomplex *s);
extern void _zrot(int n, doublecomplex x[], int incx, doublecomplex y[], int incy, double c, doublecomplex s);
extern void _zdrot(int n, doublecomplex x[], int incx, doublecomplex y[], int incy, double c, double s);
extern void _zdscal(int n, double a, doublecomplex x[], int incx);
extern void _zscal(int n, doublecomplex a, doublecomplex x[], int incx);
extern void _zswap(int n, doublecomplex x[], int incx, doublecomplex y[], int incy);
extern double _dzasum(int n, doublecomplex x[], int incx);
extern double _dznrm2(int n, doublecomplex x[], int incx);
extern double _dcabs1(doublecomplex z);
extern int _izamax(int n, doublecomplex z[], int incz);

/*
 * D1a-3. Elementary vector operations (BLAS 2)
 */
#if defined(_VLARRAY)
extern void _dgemv(char trans, int m, int n, double alpha, int lda, double a[][lda], double x[], int incx, double _beta, double y[], int incy);
extern void _dgbmv(char trans, int m, int n, int kl, int ku, double alpha, int ldab, double ab[][ldab], double x[], int incx, double _beta, double y[], int incy);
extern void _dsymv(char uplo, int n, double alpha, int lda, double a[][lda], double x[], int incx, double _beta, double y[], int incy);
extern void _dsbmv(char uplo, int n, int k, double alpha, int ldab, double ab[][ldab], double x[], int incx, double _beta, double y[], int incy);
extern void _dspmv(char uplo, int n, double alpha, double ap[], double x[], int incx, double _beta, double y[], int incy);
extern void _dtrmv(char uplo, char trans, char diag, int n, int lda, double a[][lda], double x[], int incx);
extern void _dtbmv(char uplo, char trans, char diag, int n, int k, int ldab, double ab[][ldab], double x[], int incx);
extern void _dtpmv(char uplo, char trans, char diag, int n, double ap[], double x[], int incx);
extern void _dtrsv(char uplo, char trans, char diag, int n, int lda, double a[][lda], double x[], int incx);
extern void _dtbsv(char uplo, char trans, char diag, int n, int k, int ldab, double ab[][ldab], double x[], int incx);
extern void _dtpsv(char uplo, char trans, char diag, int n, double ap[], double x[], int incx);
extern void _dger(int m, int n, double alpha, double x[], int incx, double y[], int incy, int lda, double a[][lda]);
extern void _dsyr(char uplo, int n, double alpha, double x[], int incx, int lda, double a[][lda]);
extern void _dspr(char uplo, int n, double alpha, double x[], int incx, double ap[]);
extern void _dsyr2(char uplo, int n, double alpha, double x[], int incx, double y[], int incy, int lda, double a[][lda]);
extern void _dspr2(char uplo, int n, double alpha, double x[], int incx, double y[], int incy, double ap[]);
#else
extern void _dgemv(char trans, int m, int n, double alpha, int lda, double a[], double x[], int incx, double _beta, double y[], int incy);
extern void _dgbmv(char trans, int m, int n, int kl, int ku, double alpha, int ldab, double ab[], double x[], int incx, double _beta, double y[], int incy);
extern void _dsymv(char uplo, int n, double alpha, int lda, double a[], double x[], int incx, double _beta, double y[], int incy);
extern void _dsbmv(char uplo, int n, int k, double alpha, int ldab, double ab[], double x[], int incx, double _beta, double y[], int incy);
extern void _dspmv(char uplo, int n, double alpha, double ap[], double x[], int incx, double _beta, double y[], int incy);
extern void _dtrmv(char uplo, char trans, char diag, int n, int lda, double a[], double x[], int incx);
extern void _dtbmv(char uplo, char trans, char diag, int n, int k, int ldab, double ab[], double x[], int incx);
extern void _dtpmv(char uplo, char trans, char diag, int n, double ap[], double x[], int incx);
extern void _dtrsv(char uplo, char trans, char diag, int n, int lda, double a[], double x[], int incx);
extern void _dtbsv(char uplo, char trans, char diag, int n, int k, int ldab, double ab[], double x[], int incx);
extern void _dtpsv(char uplo, char trans, char diag, int n, double ap[], double x[], int incx);
extern void _dger(int m, int n, double alpha, double x[], int incx, double y[], int incy, int lda, double a[]);
extern void _dsyr(char uplo, int n, double alpha, double x[], int incx, int lda, double a[]);
extern void _dspr(char uplo, int n, double alpha, double x[], int incx, double ap[]);
extern void _dsyr2(char uplo, int n, double alpha, double x[], int incx, double y[], int incy, int lda, double a[]);
extern void _dspr2(char uplo, int n, double alpha, double x[], int incx, double y[], int incy, double ap[]);
#endif

/*
 * D1a-4. Elementary vector operations (BLAS 2) (Complex)
 */
#if defined(_VLARRAY)
extern void _zgemv(char trans, int m, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zgbmv(char trans, int m, int n, int kl, int ku, doublecomplex alpha, int ldab, doublecomplex ab[][ldab], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zhemv(char uplo, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zhbmv(char uplo, int n, int k, doublecomplex alpha, int ldab, doublecomplex ab[][ldab], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zhpmv(char uplo, int n, doublecomplex alpha, doublecomplex ap[], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zsymv(char uplo, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zsbmv(char uplo, int n, int k, doublecomplex alpha, int ldab, doublecomplex ab[][ldab], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zspmv(char uplo, int n, doublecomplex alpha, doublecomplex ap[], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _ztrmv(char uplo, char trans, char diag, int n, int lda, doublecomplex a[][lda], doublecomplex x[], int incx);
extern void _ztbmv(char uplo, char trans, char diag, int n, int k, int ldab, doublecomplex ab[][ldab], doublecomplex x[], int incx);
extern void _ztpmv(char uplo, char trans, char diag, int n, doublecomplex ap[], doublecomplex x[], int incx);
extern void _ztrsv(char uplo, char trans, char diag, int n, int lda, doublecomplex a[][lda], doublecomplex x[], int incx);
extern void _ztbsv(char uplo, char trans, char diag, int n, int k, int ldab, doublecomplex ab[][ldab], doublecomplex x[], int incx);
extern void _ztpsv(char uplo, char trans, char diag, int n, doublecomplex ap[], doublecomplex x[], int incx);
extern void _zgeru(int m, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[][lda]);
extern void _zgerc(int m, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[][lda]);
extern void _zher(char uplo, int n, double alpha, doublecomplex x[], int incx, int lda, doublecomplex a[][lda]);
extern void _zhpr(char uplo, int n, double alpha, doublecomplex x[], int incx, doublecomplex ap[]);
extern void _zsyr(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, int lda, doublecomplex a[][lda]);
extern void _zspr(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex ap[]);
extern void _zher2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[][lda]);
extern void _zhpr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, doublecomplex ap[]);
extern void _zsyr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[][lda]);
extern void _zspr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, doublecomplex ap[]);
#else
extern void _zgemv(char trans, int m, int n, doublecomplex alpha, int lda, doublecomplex a[], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zgbmv(char trans, int m, int n, int kl, int ku, doublecomplex alpha, int ldab, doublecomplex ab[], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zhemv(char uplo, int n, doublecomplex alpha, int lda, doublecomplex a[], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zhbmv(char uplo, int n, int k, doublecomplex alpha, int ldab, doublecomplex ab[], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zhpmv(char uplo, int n, doublecomplex alpha, doublecomplex ap[], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zsymv(char uplo, int n, doublecomplex alpha, int lda, doublecomplex a[], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zsbmv(char uplo, int n, int k, doublecomplex alpha, int ldab, doublecomplex ab[], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _zspmv(char uplo, int n, doublecomplex alpha, doublecomplex ap[], doublecomplex x[], int incx, doublecomplex _beta, doublecomplex y[], int incy);
extern void _ztrmv(char uplo, char trans, char diag, int n, int lda, doublecomplex a[], doublecomplex x[], int incx);
extern void _ztbmv(char uplo, char trans, char diag, int n, int k, int ldab, doublecomplex ab[], doublecomplex x[], int incx);
extern void _ztpmv(char uplo, char trans, char diag, int n, doublecomplex ap[], doublecomplex x[], int incx);
extern void _ztrsv(char uplo, char trans, char diag, int n, int lda, doublecomplex a[], doublecomplex x[], int incx);
extern void _ztbsv(char uplo, char trans, char diag, int n, int k, int ldab, doublecomplex ab[], doublecomplex x[], int incx);
extern void _ztpsv(char uplo, char trans, char diag, int n, doublecomplex ap[], doublecomplex x[], int incx);
extern void _zgeru(int m, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[]);
extern void _zgerc(int m, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[]);
extern void _zher(char uplo, int n, double alpha, doublecomplex x[], int incx, int lda, doublecomplex a[]);
extern void _zhpr(char uplo, int n, double alpha, doublecomplex x[], int incx, doublecomplex ap[]);
extern void _zsyr(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, int lda, doublecomplex a[]);
extern void _zspr(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex ap[]);
extern void _zher2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[]);
extern void _zhpr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, doublecomplex ap[]);
extern void _zsyr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, int lda, doublecomplex a[]);
extern void _zspr2(char uplo, int n, doublecomplex alpha, doublecomplex x[], int incx, doublecomplex y[], int incy, doublecomplex ap[]);
#endif

/*
 * D1b. Elementary matrix operations (BLAS 3)
 */
#if defined(_VLARRAY)
extern void _dgemm(char transa, char transb, int m, int n, int k, double alpha, int lda, double a[][lda], int ldb, double b[][ldb], double _beta, int ldc, double c[][ldc]);
extern void _dsymm(char _side, char uplo, int m, int n, double alpha, int lda, double a[][lda], int ldb, double b[][ldb], double _beta, int ldc, double c[][ldc]);
extern void _dsyrk(char uplo, char trans, int n, int k, double alpha, int lda, double a[][lda], double _beta, int ldc, double c[][ldc]);
extern void _dsyr2k(char uplo, char trans, int n, int k, double alpha, int lda, double a[][lda], int ldb, double b[][ldb], double _beta, int ldc, double c[][ldc]);
extern void _dtrmm(char _side, char uplo, char transa, char diag, int m, int n, double alpha, int lda, double a[][lda], int ldb, double b[][ldb]);
extern void _dtrsm(char _side, char uplo, char transa, char diag, int m, int n, double alpha, int lda, double a[][lda], int ldb, double b[][ldb]);
#else
extern void _dgemm(char transa, char transb, int m, int n, int k, double alpha, int lda, double a[], int ldb, double b[], double _beta, int ldc, double c[]);
extern void _dsymm(char _side, char uplo, int m, int n, double alpha, int lda, double a[], int ldb, double b[], double _beta, int ldc, double c[]);
extern void _dsyrk(char uplo, char trans, int n, int k, double alpha, int lda, double a[], double _beta, int ldc, double c[]);
extern void _dsyr2k(char uplo, char trans, int n, int k, double alpha, int lda, double a[], int ldb, double b[], double _beta, int ldc, double c[]);
extern void _dtrmm(char _side, char uplo, char transa, char diag, int m, int n, double alpha, int lda, double a[], int ldb, double b[]);
extern void _dtrsm(char _side, char uplo, char transa, char diag, int m, int n, double alpha, int lda, double a[], int ldb, double b[]);
#endif

/*
 * D1b-2. Elementary matrix operations (BLAS 3) (Complex)
 */
#if defined(_VLARRAY)
extern void _zgemm(char transa, char transb, int m, int n, int k, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex _beta, int ldc, doublecomplex c[][ldc]);
extern void _zsymm(char _side, char uplo, int m, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex _beta, int ldc, doublecomplex c[][ldc]);
extern void _zhemm(char _side, char uplo, int m, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex _beta, int ldc, doublecomplex c[][ldc]);
extern void _zsyrk(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[][lda], doublecomplex _beta, int ldc, doublecomplex c[][ldc]);
extern void _zherk(char uplo, char trans, int n, int k, double alpha, int lda, doublecomplex a[][lda], double _beta, int ldc, doublecomplex c[][ldc]);
extern void _zsyr2k(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex _beta, int ldc, doublecomplex c[][ldc]);
extern void _zher2k(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double _beta, int ldc, doublecomplex c[][ldc]);
extern void _ztrmm(char _side, char uplo, char transa, char diag, int m, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb]);
extern void _ztrsm(char _side, char uplo, char transa, char diag, int m, int n, doublecomplex alpha, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb]);
#else
extern void _zgemm(char transa, char transb, int m, int n, int k, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex _beta, int ldc, doublecomplex c[]);
extern void _zsymm(char _side, char uplo, int m, int n, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex _beta, int ldc, doublecomplex c[]);
extern void _zhemm(char _side, char uplo, int m, int n, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex _beta, int ldc, doublecomplex c[]);
extern void _zsyrk(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[], doublecomplex _beta, int ldc, doublecomplex c[]);
extern void _zherk(char uplo, char trans, int n, int k, double alpha, int lda, doublecomplex a[], double _beta, int ldc, doublecomplex c[]);
extern void _zsyr2k(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex _beta, int ldc, doublecomplex c[]);
void _zher2k(char uplo, char trans, int n, int k, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[], double _beta, int ldc, doublecomplex c[]);
extern void _ztrmm(char _side, char uplo, char transa, char diag, int m, int n, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[]);
extern void _ztrsm(char _side, char uplo, char transa, char diag, int m, int n, doublecomplex alpha, int lda, doublecomplex a[], int ldb, doublecomplex b[]);
#endif

/*
 * D1b-3. Elementary matrix operations (LAPACK auxiliary routines)
 */
#if defined(_VLARRAY)
extern double _dlange(char norm, int m, int n, int lda, double a[][lda], double work[]);
extern double _dlangb(char norm, int n, int kl, int ku, int ldab, double ab[][ldab], double work[]);
extern double _dlangt(char norm, int n, double dl[], double d[], double du[]);
extern double _dlansy(char norm, char uplo, int n, int lda, double a[][lda], double work[]);
extern double _dlansb(char norm, char uplo, int n, int k, int ldab, double ab[][ldab], double work[]);
extern double _dlansp(char norm, char uplo, int n, double ap[], double work[]);
extern double _dlanst(char norm, int n, double d[], double e[]);
extern double _dlantr(char norm, char uplo, char diag, int m, int n, int lda, double a[][lda], double work[]);
#else
extern double _dlange(char norm, int m, int n, int lda, double a[], double work[]);
extern double _dlangb(char norm, int n, int kl, int ku, int ldab, double ab[], double work[]);
extern double _dlangt(char norm, int n, double dl[], double d[], double du[]);
extern double _dlansy(char norm, char uplo, int n, int lda, double a[], double work[]);
extern double _dlansb(char norm, char uplo, int n, int k, int ldab, double ab[], double work[]);
extern double _dlansp(char norm, char uplo, int n, double ap[], double work[]);
extern double _dlanst(char norm, int n, double d[], double e[]);
extern double _dlantr(char norm, char uplo, char diag, int m, int n, int lda, double a[], double work[]);
#endif

/*
 * D1b-4. Elementary matrix operations (LAPACK auxiliary routines) (Complex)
 */
#if defined(_VLARRAY)
extern double _zlange(char norm, int m, int n, int lda, doublecomplex a[][lda], double work[]);
extern double _zlangb(char norm, int n, int kl, int ku, int ldab, doublecomplex ab[][ldab], double work[]);
extern double _zlangt(char norm, int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[]);
extern double _zlansy(char norm, char uplo, int n, int lda, doublecomplex a[][lda], double work[]);
extern double _zlansb(char norm, char uplo, int n, int k, int ldab, doublecomplex ab[][ldab], double work[]);
extern double _zlansp(char norm, char uplo, int n, doublecomplex ap[], double work[]);
extern double _zlanhe(char norm, char uplo, int n, int lda, doublecomplex a[][lda], double work[]);
extern double _zlanhb(char norm, char uplo, int n, int k, int ldab, doublecomplex ab[][ldab], double work[]);
extern double _zlanhp(char norm, char uplo, int n, doublecomplex ap[], double work[]);
extern double _zlanht(char norm, int n, double d[], doublecomplex e[]);
extern double _zlantr(char norm, char uplo, char diag, int m, int n, int lda, doublecomplex a[][lda], double work[]);
#else
extern double _zlange(char norm, int m, int n, int lda, doublecomplex a[], double work[]);
extern double _zlangb(char norm, int n, int kl, int ku, int ldab, doublecomplex ab[], double work[]);
extern double _zlangt(char norm, int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[]);
extern double _zlansy(char norm, char uplo, int n, int lda, doublecomplex a[], double work[]);
extern double _zlansb(char norm, char uplo, int n, int k, int ldab, doublecomplex ab[], double work[]);
extern double _zlansp(char norm, char uplo, int n, doublecomplex ap[], double work[]);
extern double _zlanhe(char norm, char uplo, int n, int lda, doublecomplex a[], double work[]);
extern double _zlanhb(char norm, char uplo, int n, int k, int ldab, doublecomplex ab[], double work[]);
extern double _zlanhp(char norm, char uplo, int n, doublecomplex ap[], double work[]);
extern double _zlanht(char norm, int n, double d[], doublecomplex e[]);
extern double _zlantr(char norm, char uplo, char diag, int m, int n, int lda, doublecomplex a[], double work[]);
#endif

/*
 * D2a. Solution of systems of _linear equations (real nonsymmetric matrices)
 */
#if defined(_VLARRAY)
extern void _dgesv(int n, int nrhs, int lda, double a[][lda], int ipiv[], int ldb, double b[][ldb], int *info);
extern void _dgetrf(int m, int n, int lda, double a[][lda], int ipiv[], int *info);
extern void _dgetrs(char trans, int n, int nrhs, int lda, double a[][lda], int ipiv[], int ldb, double b[][ldb], int *info);
extern void _dgetri(int n, int lda, double a[][lda], int ipiv[], double work[], int lwork, int *info);
extern void _dgesvx(char fact, char trans, int n, int nrhs, int lda, double a[][lda], int ldaf, double af[][ldaf], int ipiv[], char *equed, double r[], double c[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dsgesv(int n, int nrhs, int lda, double a[][lda], int ipiv[], int ldb, double b[][ldb], int ldx, double x[][ldx], double work[], float swork[], int *iter, int *info);
extern void _dgecon(char norm, int n, int lda, double a[][lda], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void _dgbsv(int n, int kl, int ku, int nrhs, int ldab, double ab[][ldab], int ipiv[], int ldb, double b[][ldb], int *info);
extern void _dgbtrf(int m, int n, int kl, int ku, int ldab, double ab[][ldab], int ipiv[], int *info);
extern void _dgbtrs(char trans, int n, int kl, int ku, int nrhs, int ldab, double ab[][ldab], int ipiv[], int ldb, double b[][ldb], int *info);
extern void _dgbsvx(char fact, char trans, int n, int kl, int ku, int nrhs, int ldab, double ab[][ldab], int ldafb, double afb[][ldafb], int ipiv[], char *equed, double r[], double c[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dgbcon(char norm, int n, int kl, int ku, int ldab, double ab[][ldab], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void _dgtsv(int n, int nrhs, double dl[], double d[], double du[], int ldb, double b[][ldb], int *info);
extern void _dgttrf(int n, double dl[], double d[], double du[], double du2[], int ipiv[], int *info);
extern void _dgttrs(char trans, int n, int nrhs, double dl[], double d[], double du[], double du2[], int ipiv[], int ldb, double b[][ldb], int *info);
extern void _dgtsvx(char fact, char trans, int n, int nrhs, double dl[], double d[], double du[], double dlf[], double df[], double duf[], double du2[], int ipiv[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dgtcon(char norm, int n, double dl[], double d[], double du[], double du2[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
#else
extern void _dgesv(int n, int nrhs, int lda, double a[], int ipiv[], int ldb, double b[], int *info);
extern void _dgetrf(int m, int n, int lda, double a[], int ipiv[], int *info);
extern void _dgetrs(char trans, int n, int nrhs, int lda, double a[], int ipiv[], int ldb, double b[], int *info);
extern void _dgetri(int n, int lda, double a[], int ipiv[], double work[], int lwork, int *info);
extern void _dgesvx(char fact, char trans, int n, int nrhs, int lda, double a[], int ldaf, double af[], int ipiv[], char *equed, double r[], double c[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dsgesv(int n, int nrhs, int lda, double a[], int ipiv[], int ldb, double b[], int ldx, double x[], double work[], float swork[], int *iter, int *info);
extern void _dgecon(char norm, int n, int lda, double a[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void _dgbsv(int n, int kl, int ku, int nrhs, int ldab, double ab[], int ipiv[], int ldb, double b[], int *info);
extern void _dgbtrf(int m, int n, int kl, int ku, int ldab, double ab[], int ipiv[], int *info);
extern void _dgbtrs(char trans, int n, int kl, int ku, int nrhs, int ldab, double ab[], int ipiv[], int ldb, double b[], int *info);
extern void _dgbsvx(char fact, char trans, int n, int kl, int ku, int nrhs, int ldab, double ab[], int ldafb, double afb[], int ipiv[], char *equed, double r[], double c[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dgbcon(char norm, int n, int kl, int ku, int ldab, double ab[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void _dgtsv(int n, int nrhs, double dl[], double d[], double du[], int ldb, double b[], int *info);
extern void _dgttrf(int n, double dl[], double d[], double du[], double du2[], int ipiv[], int *info);
extern void _dgttrs(char trans, int n, int nrhs, double dl[], double d[], double du[], double du2[], int ipiv[], int ldb, double b[], int *info);
extern void _dgtsvx(char fact, char trans, int n, int nrhs, double dl[], double d[], double du[], double dlf[], double df[], double duf[], double du2[], int ipiv[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dgtcon(char norm, int n, double dl[], double d[], double du[], double du2[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
#endif

/*
 * D2a3. Solution of systems of _linear equations (real triangular matrices)
 */
#if defined(_VLARRAY)
extern void _dtrtrs(char uplo, char trans, char diag, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], int *info);
extern void _dtrtri(char uplo, char diag, int n, int lda, double a[][lda], int *info);
extern void _dtrcon(char norm, char uplo, char diag, int n, int lda, double a[][lda], double *rcond, double work[], int iwork[], int *info);
extern void _dtptrs(char uplo, char trans, char diag, int n, int nrhs, double ap[], int ldb, double b[][ldb], int *info);
extern void _dtptri(char uplo, char diag, int n, double ap[], int *info);
extern void _dtpcon(char norm, char uplo, char diag, int n, double ap[], double *rcond, double work[], int iwork[], int *info);
extern void _dtbtrs(char uplo, char trans, char diag, int n, int kd, int nrhs, int ldab, double ab[][ldab], int ldb, double b[][ldb], int *info);
extern void _dtbcon(char norm, char uplo, char diag, int n, int kd, int ldab, double ab[][ldab], double *rcond, double work[], int iwork[], int *info);
#else
extern void _dtrtrs(char uplo, char trans, char diag, int n, int nrhs, int lda, double a[], int ldb, double b[], int *info);
extern void _dtrtri(char uplo, char diag, int n, int lda, double a[], int *info);
extern void _dtrcon(char norm, char uplo, char diag, int n, int lda, double a[], double *rcond, double work[], int iwork[], int *info);
extern void _dtptrs(char uplo, char trans, char diag, int n, int nrhs, double ap[], int ldb, double b[], int *info);
extern void _dtptri(char uplo, char diag, int n, double ap[], int *info);
extern void _dtpcon(char norm, char uplo, char diag, int n, double ap[], double *rcond, double work[], int iwork[], int *info);
extern void _dtbtrs(char uplo, char trans, char diag, int n, int kd, int nrhs, int ldab, double ab[], int ldb, double b[], int *info);
extern void _dtbcon(char norm, char uplo, char diag, int n, int kd, int ldab, double ab[], double *rcond, double work[], int iwork[], int *info);
#endif

/*
 * D2b1a. Solution of systems of _linear equations (real symmetric indefinite matrices)
 */
#if defined(_VLARRAY)
extern void _dsysv(char uplo, int n, int nrhs, int lda, double a[][lda], int ipiv[], int ldb, double b[][ldb], double work[], int lwork, int *info);
extern void _dsytrf(char uplo, int n, int lda, double a[][lda], int ipiv[], double work[], int lwork, int *info);
extern void _dsytrs(char uplo, int n, int nrhs, int lda, double a[][lda], int ipiv[], int ldb, double b[][ldb], int *info);
extern void _dsytri(char uplo, int n, int lda, double a[][lda], int ipiv[], double work[], int *info);
extern void _dsysvx(char fact, char uplo, int n, int nrhs, int lda, double a[][lda], int ldaf, double af[][ldaf], int ipiv[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int lwork, int iwork[], int *info);
extern void _dsycon(char uplo, int n, int lda, double a[][lda], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void _dspsv(char uplo, int n, int nrhs, double ap[], int ipiv[], int ldb, double b[][ldb], int *info);
extern void _dsptrf(char uplo, int n, double ap[], int ipiv[], int *info);
extern void _dsptrs(char uplo, int n, int nrhs, double ap[], int ipiv[], int ldb, double b[][ldb], int *info);
extern void _dsptri(char uplo, int n, double ap[], int ipiv[], double work[], int *info);
extern void _dspsvx(char fact, char uplo, int n, int nrhs, double ap[], double afp[], int ipiv[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dspcon(char uplo, int n, double ap[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
#else
extern void _dsysv(char uplo, int n, int nrhs, int lda, double a[], int ipiv[], int ldb, double b[], double work[], int lwork, int *info);
extern void _dsytrf(char uplo, int n, int lda, double a[], int ipiv[], double work[], int lwork, int *info);
extern void _dsytrs(char uplo, int n, int nrhs, int lda, double a[], int ipiv[], int ldb, double b[], int *info);
extern void _dsytri(char uplo, int n, int lda, double a[], int ipiv[], double work[], int *info);
extern void _dsysvx(char fact, char uplo, int n, int nrhs, int lda, double a[], int ldaf, double af[], int ipiv[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int lwork, int iwork[], int *info);
extern void _dsycon(char uplo, int n, int lda, double a[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void _dspsv(char uplo, int n, int nrhs, double ap[], int ipiv[], int ldb, double b[], int *info);
extern void _dsptrf(char uplo, int n, double ap[], int ipiv[], int *info);
extern void _dsptrs(char uplo, int n, int nrhs, double ap[], int ipiv[], int ldb, double b[], int *info);
extern void _dsptri(char uplo, int n, double ap[], int ipiv[], double work[], int *info);
extern void _dspsvx(char fact, char uplo, int n, int nrhs, double ap[], double afp[], int ipiv[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dspcon(char uplo, int n, double ap[], int ipiv[], double anorm, double *rcond, double work[], int iwork[], int *info);
#endif

/*
 * D2b1b. Solution of systems of _linear equations (real symmetric positive definite matrices)
 */
#if defined(_VLARRAY)
extern void _dposv(char uplo, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], int *info);
extern void _dpotrf(char uplo, int n, int lda, double a[][lda], int *info);
extern void _dpotrs(char uplo, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], int *info);
extern void _dpotri(char uplo, int n, int lda, double a[][lda], int *info);
extern void _dposvx(char fact, char uplo, int n, int nrhs, int lda, double a[][lda], int ldaf, double af[][ldaf], char *equed, double s[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dsposv(char uplo, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], int ldx, double x[][ldx], double work[], float swork[], int *iter, int *info);
extern void _dpocon(char uplo, int n, int lda, double a[][lda], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void _dppsv(char uplo, int n, int nrhs, double ap[], int ldb, double b[][ldb], int *info);
extern void _dpptrf(char uplo, int n, double ap[], int *info);
extern void _dpptrs(char uplo, int n, int nrhs, double ap[], int ldb, double b[][ldb], int *info);
extern void _dpptri(char uplo, int n, double ap[], int *info);
extern void _dppsvx(char fact, char uplo, int n, int nrhs, double ap[], double afp[], char *equed, double s[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dppcon(char uplo, int n, double ap[], double anorm, double *rcond, double work[], int iwork[], int *info);
#else
extern void _dposv(char uplo, int n, int nrhs, int lda, double a[], int ldb, double b[], int *info);
extern void _dpotrf(char uplo, int n, int lda, double a[], int *info);
extern void _dpotrs(char uplo, int n, int nrhs, int lda, double a[], int ldb, double b[], int *info);
extern void _dpotri(char uplo, int n, int lda, double a[], int *info);
extern void _dposvx(char fact, char uplo, int n, int nrhs, int lda, double a[], int ldaf, double af[], char *equed, double s[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dsposv(char uplo, int n, int nrhs, int lda, double a[], int ldb, double b[], int ldx, double x[], double work[], float swork[], int *iter, int *info);
extern void _dpocon(char uplo, int n, int lda, double a[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void _dppsv(char uplo, int n, int nrhs, double ap[], int ldb, double b[], int *info);
extern void _dpptrf(char uplo, int n, double ap[], int *info);
extern void _dpptrs(char uplo, int n, int nrhs, double ap[], int ldb, double b[], int *info);
extern void _dpptri(char uplo, int n, double ap[], int *info);
extern void _dppsvx(char fact, char uplo, int n, int nrhs, double ap[], double afp[], char *equed, double s[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int iwork[], int *info);
extern void _dppcon(char uplo, int n, double ap[], double anorm, double *rcond, double work[], int iwork[], int *info);
#endif

/*
 * D2b2. Solution of systems of _linear equations (real symmetric positive definite band matrices)
 */
#if defined(_VLARRAY)
extern void _dpbsv(char uplo, int n, int kd, int nrhs, int ldab, double ab[][ldab], int ldb, double b[][ldb], int *info);
extern void _dpbtrf(char uplo, int n, int kd, int ldab, double ab[][ldab], int *info);
extern void _dpbtrs(char uplo, int n, int kd, int nrhs, int ldab, double ab[][ldab], int ldb, double b[][ldb], int *info);
extern void _dpbsvx(char fact, char uplo, int n, int kd, int nrhs, int ldab, double ab[][ldab], int ldafb, double afb[][ldafb], char *equed, double s[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[],  double berr[], double work[], int iwork[], int *info);
extern void _dpbcon(char uplo, int n, int kd, int ldab, double ab[][ldab], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void _dptsv(int n, int nrhs, double d[], double e[], int ldb, double b[][ldb], int *info);
extern void _dpttrf(int n, double d[], double e[], int *info);
extern void _dpttrs(int n, int nrhs, double d[], double e[], int ldb, double b[][ldb], int *info);
extern void _dptsvx(char fact, int n, int nrhs, double d[], double e[], double df[], double ef[], int ldb, double b[][ldb], int ldx, double x[][ldx], double *rcond, double ferr[], double berr[], double work[], int *info);
extern void _dptcon(int n, double d[], double e[], double anorm, double *rcond, double work[], int *info);
#else
extern void _dpbsv(char uplo, int n, int kd, int nrhs, int ldab, double ab[], int ldb, double b[], int *info);
extern void _dpbtrf(char uplo, int n, int kd, int ldab, double ab[], int *info);
extern void _dpbtrs(char uplo, int n, int kd, int nrhs, int ldab, double ab[], int ldb, double b[], int *info);
extern void _dpbsvx(char fact, char uplo, int n, int kd, int nrhs, int ldab, double ab[], int ldafb, double afb[], char *equed, double s[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[],  double berr[], double work[], int iwork[], int *info);
extern void _dpbcon(char uplo, int n, int kd, int ldab, double ab[], double anorm, double *rcond, double work[], int iwork[], int *info);
extern void _dptsv(int n, int nrhs, double d[], double e[], int ldb, double b[], int *info);
extern void _dpttrf(int n, double d[], double e[], int *info);
extern void _dpttrs(int n, int nrhs, double d[], double e[], int ldb, double b[], int *info);
extern void _dptsvx(char fact, int n, int nrhs, double d[], double e[], double df[], double ef[], int ldb, double b[], int ldx, double x[], double *rcond, double ferr[], double berr[], double work[], int *info);
extern void _dptcon(int n, double d[], double e[], double anorm, double *rcond, double work[], int *info);
#endif

/*
 * D2c. Solution of systems of _linear equations (complex non-Hermitian matrices)
 */
#if defined(_VLARRAY)
extern void _zgesv(int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zgetrf(int m, int n, int lda, doublecomplex a[][lda], int ipiv[], int *info);
extern void _zgetrs(char trans, int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zgetri(int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void _zgesvx(char fact, char trans, int n, int nrhs, int lda, doublecomplex a[][lda], int ldaf, doublecomplex af[][ldaf], int ipiv[], char *equed, double r[], double c[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zcgesv(int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], doublecomplex work[], floatcomplex swork[], double rwork[], int *iter, int *info);
extern void _zgecon(char norm, int n, int lda, doublecomplex a[][lda], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _zgbsv(int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[][ldab], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zgbtrf(int m, int n, int kl, int ku, int ldab, doublecomplex ab[][ldab], int ipiv[], int *info);
extern void _zgbtrs(char trans, int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[][ldab], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zgbsvx(char fact, char trans, int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[][ldab], int ldafb, doublecomplex afb[][ldafb], int ipiv[], char *equed, double r[], double c[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zgbcon(char norm, int n, int kl, int ku, int ldab, doublecomplex ab[][ldab], int ipiv[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _zgtsv(int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zgttrf(int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], int *info);
extern void _zgttrs(char trans, int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zgtsvx(char fact, char trans, int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex dlf[], doublecomplex df[], doublecomplex duf[], doublecomplex du2[], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zgtcon(char norm, int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void _zsysv(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], doublecomplex work[], int lwork, int *info);
extern void _zsytrf(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void _zsytrs(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zsytri(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex work[], int *info);
extern void _zsysvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldaf, doublecomplex af[][ldaf], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zsycon(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void _zspsv(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zsptrf(char uplo, int n, doublecomplex ap[], int ipiv[], int *info);
extern void _zsptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zsptri(char uplo, int n, doublecomplex ap[], int ipiv[], doublecomplex work[], int *info);
extern void _zspsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zspcon(char uplo, int n, doublecomplex ap[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
#else
extern void _zgesv(int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void _zgetrf(int m, int n, int lda, doublecomplex a[], int ipiv[], int *info);
extern void _zgetrs(char trans, int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void _zgetri(int n, int lda, doublecomplex a[], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void _zgesvx(char fact, char trans, int n, int nrhs, int lda, doublecomplex a[], int ldaf, doublecomplex af[], int ipiv[], char *equed, double r[], double c[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zcgesv(int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], doublecomplex work[], floatcomplex swork[], double rwork[], int *iter, int *info);
extern void _zgecon(char norm, int n, int lda, doublecomplex a[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _zgbsv(int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void _zgbtrf(int m, int n, int kl, int ku, int ldab, doublecomplex ab[], int ipiv[], int *info);
extern void _zgbtrs(char trans, int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void _zgbsvx(char fact, char trans, int n, int kl, int ku, int nrhs, int ldab, doublecomplex ab[], int ldafb, doublecomplex afb[], int ipiv[], char *equed, double r[], double c[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zgbcon(char norm, int n, int kl, int ku, int ldab, doublecomplex ab[], int ipiv[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _zgtsv(int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], int ldb, doublecomplex b[], int *info);
extern void _zgttrf(int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], int *info);
extern void _zgttrs(char trans, int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void _zgtsvx(char fact, char trans, int n, int nrhs, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex dlf[], doublecomplex df[], doublecomplex duf[], doublecomplex du2[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zgtcon(char norm, int n, doublecomplex dl[], doublecomplex d[], doublecomplex du[], doublecomplex du2[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void _zsysv(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], doublecomplex work[], int lwork, int *info);
extern void _zsytrf(char uplo, int n, int lda, doublecomplex a[], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void _zsytrs(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void _zsytri(char uplo, int n, int lda, doublecomplex a[], int ipiv[], doublecomplex work[], int *info);
extern void _zsysvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldaf, doublecomplex af[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zsycon(char uplo, int n, int lda, doublecomplex a[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void _zspsv(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void _zsptrf(char uplo, int n, doublecomplex ap[], int ipiv[], int *info);
extern void _zsptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void _zsptri(char uplo, int n, doublecomplex ap[], int ipiv[], doublecomplex work[], int *info);
extern void _zspsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zspcon(char uplo, int n, doublecomplex ap[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
#endif

/*
 * D2c3. Solution of systems of _linear equations (complex triangular matrices)
 */
#if defined(_VLARRAY)
extern void _ztrtrs(char uplo, char trans, char diag, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *info);
extern void _ztrtri(char uplo, char diag, int n, int lda, doublecomplex a[][lda], int *info);
extern void _ztrcon(char norm, char uplo, char diag, int n, int lda, doublecomplex a[][lda], double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _ztptrs(char uplo, char trans, char diag, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[][ldb], int *info);
extern void _ztptri(char uplo, char diag, int n, doublecomplex ap[], int *info);
extern void _ztpcon(char norm, char uplo, char diag, int n, doublecomplex ap[], double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _ztbtrs(char uplo, char trans, char diag, int n, int kd, int nrhs, int ldab, doublecomplex ab[][ldab], int ldb, doublecomplex b[][ldb], int *info);
extern void _ztbcon(char norm, char uplo, char diag, int n, int kd, int ldab, doublecomplex ab[][ldab], double *rcond, doublecomplex work[], double rwork[], int *info);
#else
extern void _ztrtrs(char uplo, char trans, char diag, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *info);
extern void _ztrtri(char uplo, char diag, int n, int lda, doublecomplex a[], int *info);
extern void _ztrcon(char norm, char uplo, char diag, int n, int lda, doublecomplex a[], double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _ztptrs(char uplo, char trans, char diag, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[], int *info);
extern void _ztptri(char uplo, char diag, int n, doublecomplex ap[], int *info);
extern void _ztpcon(char norm, char uplo, char diag, int n, doublecomplex ap[], double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _ztbtrs(char uplo, char trans, char diag, int n, int kd, int nrhs, int ldab, doublecomplex ab[], int ldb, doublecomplex b[], int *info);
extern void _ztbcon(char norm, char uplo, char diag, int n, int kd, int ldab, doublecomplex ab[], double *rcond, doublecomplex work[], double rwork[], int *info);
#endif

/*
 * D2d1a. Solution of systems of _linear equations (Hermitian matrices)
 */
#if defined(_VLARRAY)
extern void _zhesv(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], doublecomplex work[], int lwork, int *info);
extern void _zhetrf(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void _zhetrs(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zhetri(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex work[], int *info);
extern void _zhesvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldaf, doublecomplex af[][ldaf], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zhecon(char uplo, int n, int lda, doublecomplex a[][lda], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void _zhpsv(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zhptrf(char uplo, int n, doublecomplex ap[], int ipiv[], int *info);
extern void _zhptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zhptri(char uplo, int n, doublecomplex ap[], int ipiv[], doublecomplex work[], int *info);
extern void _zhpsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], int ipiv[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zhpcon(char uplo, int n, doublecomplex ap[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
#else
extern void _zhesv(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], doublecomplex work[], int lwork, int *info);
extern void _zhetrf(char uplo, int n, int lda, doublecomplex a[], int ipiv[], doublecomplex work[], int lwork, int *info);
extern void _zhetrs(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void _zhetri(char uplo, int n, int lda, doublecomplex a[], int ipiv[], doublecomplex work[], int *info);
extern void _zhesvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldaf, doublecomplex af[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zhecon(char uplo, int n, int lda, doublecomplex a[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
extern void _zhpsv(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void _zhptrf(char uplo, int n, doublecomplex ap[], int ipiv[], int *info);
extern void _zhptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ipiv[], int ldb, doublecomplex b[], int *info);
extern void _zhptri(char uplo, int n, doublecomplex ap[], int ipiv[], doublecomplex work[], int *info);
extern void _zhpsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], int ipiv[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zhpcon(char uplo, int n, doublecomplex ap[], int ipiv[], double anorm, double *rcond, doublecomplex work[], int *info);
#endif

/*
 * D2d1b. Solution of systems of _linear equations (Hermitian positive definite matrices)
 */
#if defined(_VLARRAY)
extern void _zposv(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *info);
extern void _zpotrf(char uplo, int n, int lda, doublecomplex a[][lda], int *info);
extern void _zpotrs(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *info);
extern void _zpotri(char uplo, int n, int lda, doublecomplex a[][lda], int *info);
extern void _zposvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldaf, doublecomplex af[][ldaf], char *equed, double s[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zcposv(char uplo, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], doublecomplex work[], floatcomplex swork[], double rwork[], int *iter, int *info);
extern void _zpocon(char uplo, int n, int lda, doublecomplex a[][lda], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _zppsv(char uplo, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zpptrf(char uplo, int n, doublecomplex ap[], int *info);
extern void _zpptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zpptri(char uplo, int n, doublecomplex ap[], int *info);
extern void _zppsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], char *equed, double s[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zppcon(char uplo, int n, doublecomplex ap[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
#else
extern void _zposv(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *info);
extern void _zpotrf(char uplo, int n, int lda, doublecomplex a[], int *info);
extern void _zpotrs(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *info);
extern void _zpotri(char uplo, int n, int lda, doublecomplex a[], int *info);
extern void _zposvx(char fact, char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldaf, doublecomplex af[], char *equed, double s[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zcposv(char uplo, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], doublecomplex work[], floatcomplex swork[], double rwork[], int *iter, int *info);
extern void _zpocon(char uplo, int n, int lda, doublecomplex a[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _zppsv(char uplo, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[], int *info);
extern void _zpptrf(char uplo, int n, doublecomplex ap[], int *info);
extern void _zpptrs(char uplo, int n, int nrhs, doublecomplex ap[], int ldb, doublecomplex b[], int *info);
extern void _zpptri(char uplo, int n, doublecomplex ap[], int *info);
extern void _zppsvx(char fact, char uplo, int n, int nrhs, doublecomplex ap[], doublecomplex afp[], char *equed, double s[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zppcon(char uplo, int n, doublecomplex ap[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
#endif

/*
 * D2d2. Solution of systems of _linear equations (Hermitian positive definite band matrices)
 */
#if defined(_VLARRAY)
extern void _zpbsv(char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[][ldab], int ldb, doublecomplex b[][ldb], int *info);
extern void _zpbtrf(char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], int *info);
extern void _zpbtrs(char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[][ldab], int ldb, doublecomplex b[][ldb], int *info);
extern void _zpbsvx(char fact, char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[][ldab], int ldafb, doublecomplex afb[][ldafb], char *equed, double s[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[],  double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zpbcon(char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _zptsv(int n, int nrhs, double d[], doublecomplex e[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zpttrf(int n, double d[], doublecomplex e[], int *info);
extern void _zpttrs(char uplo, int n, int nrhs, double d[], doublecomplex e[], int ldb, doublecomplex b[][ldb], int *info);
extern void _zptsvx(char fact, int n, int nrhs, double d[], doublecomplex e[], double df[], doublecomplex ef[], int ldb, doublecomplex b[][ldb], int ldx, doublecomplex x[][ldx], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zptcon(int n, double d[], doublecomplex e[], double anorm, double *rcond, double rwork[], int *info);
#else
extern void _zpbsv(char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[], int ldb, doublecomplex b[], int *info);
extern void _zpbtrf(char uplo, int n, int kd, int ldab, doublecomplex ab[], int *info);
extern void _zpbtrs(char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[], int ldb, doublecomplex b[], int *info);
extern void _zpbsvx(char fact, char uplo, int n, int kd, int nrhs, int ldab, doublecomplex ab[], int ldafb, doublecomplex afb[], char *equed, double s[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[],  double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zpbcon(char uplo, int n, int kd, int ldab, doublecomplex ab[], double anorm, double *rcond, doublecomplex work[], double rwork[], int *info);
extern void _zptsv(int n, int nrhs, double d[], doublecomplex e[], int ldb, doublecomplex b[], int *info);
extern void _zpttrf(int n, double d[], doublecomplex e[], int *info);
extern void _zpttrs(char uplo, int n, int nrhs, double d[], doublecomplex e[], int ldb, doublecomplex b[], int *info);
extern void _zptsvx(char fact, int n, int nrhs, double d[], doublecomplex e[], double df[], doublecomplex ef[], int ldb, doublecomplex b[], int ldx, doublecomplex x[], double *rcond, double ferr[], double berr[], doublecomplex work[], double rwork[], int *info);
extern void _zptcon(int n, double d[], doublecomplex e[], double anorm, double *rcond, double rwork[], int *info);
#endif

/*
 * D4a1. Real symmetric matrix _eigenvalue problems
 */
#if defined(_VLARRAY)
extern void _dsyev(char jobz, char uplo, int n,  int lda, double a[][lda], double w[], double work[], int lwork, int *info);
extern void _dsyevd(char jobz, char uplo, int n, int lda, double a[][lda], double w[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dsyevr(char jobz, char range, char uplo, int n, int lda, double a[][lda], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], int isuppz[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dsyevx(char jobz, char range, char uplo, int n, int lda, double a[][lda], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int ifail[], int *info);
extern void _dspev(char jobz, char uplo, int n, double ap[], double w[], int ldz, double z[][ldz], double work[], int *info);
extern void _dspevd(char jobz, char uplo, int n, double ap[], double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dspevx(char jobz, char range, char uplo, int n, double ap[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void _dsbev(char jobz, char uplo, int n, int kd, int ldab, double ab[][ldab], double w[], int ldz, double z[][ldz], double work[], int *info);
extern void _dsbevd(char jobz, char uplo, int n, int kd, int ldab, double ab[][ldab], double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dsbevx(char jobz, char range, char uplo, int n, int kd, int ldab, double ab[][ldab], int ldq, double q[][ldq], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void _dsbtrd(char vect, char uplo, int n, int kd, int ldab, double ab[][ldab], double d[], double e[], int ldq, double q[][ldq], double work[], int *info);
extern void _dstev(char jobz, int n, double d[], double e[], int ldz, double z[][ldz], double work[], int *info);
extern void _dstevd(char jobz, int n, double d[], double e[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dstevr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], int isuppz[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dstevx(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void _dsytrd(char uplo, int n, int lda, double a[][lda], double d[], double e[], double tau[], double work[], int lwork, int *info);
extern void _dorgtr(char uplo, int n, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void _dormtr(char _side, char uplo, char trans, int m, int n, int lda, double a[][lda], double tau[], int ldc, double c[][ldc], double work[], int lwork, int *info);
extern void _dopgtr(char uplo, int n, double ap[], double tau[], int ldq, double q[][ldq],  double work[], int *info);
extern void _dopmtr(char _side, char uplo, char trans, int m, int n, double ap[], double tau[], int ldc, double c[][ldc], double work[], int *info);
extern void _dsteqr(char compz, int n, double d[], double e[], int ldz, double z[][ldz], double work[], int *info);
extern void _dstedc(char compz, int n, double d[], double e[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dstemr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, int *m, double w[], int ldz, double z[][ldz], int nzc, int isuppz[], int *tryrac, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dstein(int n, double d[], double e[], int m, double w[], int iblock[], int isplit[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void _dpteqr(char compz, int n, double d[], double e[], int ldz, double z[][ldz], double work[], int *info);
#else
extern void _dsyev(char jobz, char uplo, int n,  int lda, double a[], double w[], double work[], int lwork, int *info);
extern void _dsyevd(char jobz, char uplo, int n, int lda, double a[], double w[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dsyevr(char jobz, char range, char uplo, int n, int lda, double a[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], int isuppz[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dsyevx(char jobz, char range, char uplo, int n, int lda, double a[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int lwork, int iwork[], int ifail[], int *info);
extern void _dspev(char jobz, char uplo, int n, double ap[], double w[], int ldz, double z[], double work[], int *info);
extern void _dspevd(char jobz, char uplo, int n, double ap[], double w[], int ldz, double z[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dspevx(char jobz, char range, char uplo, int n, double ap[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
extern void _dsbev(char jobz, char uplo, int n, int kd, int ldab, double ab[], double w[], int ldz, double z[], double work[], int *info);
extern void _dsbevd(char jobz, char uplo, int n, int kd, int ldab, double ab[], double w[], int ldz, double z[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dsbevx(char jobz, char range, char uplo, int n, int kd, int ldab, double ab[], int ldq, double q[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
extern void _dsbtrd(char vect, char uplo, int n, int kd, int ldab, double ab[], double d[], double e[], int ldq, double q[], double work[], int *info);
extern void _dstev(char jobz, int n, double d[], double e[], int ldz, double z[], double work[], int *info);
extern void _dstevd(char jobz, int n, double d[], double e[], int ldz, double z[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dstevr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], int isuppz[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dstevx(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
extern void _dsytrd(char uplo, int n, int lda, double a[], double d[], double e[], double tau[], double work[], int lwork, int *info);
extern void _dorgtr(char uplo, int n, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void _dormtr(char _side, char uplo, char trans, int m, int n, int lda, double a[], double tau[], int ldc, double c[], double work[], int lwork, int *info);
extern void _dopgtr(char uplo, int n, double ap[], double tau[], int ldq, double q[],  double work[], int *info);
extern void _dopmtr(char _side, char uplo, char trans, int m, int n, double ap[], double tau[], int ldc, double c[], double work[], int *info);
extern void _dsteqr(char compz, int n, double d[], double e[], int ldz, double z[], double work[], int *info);
extern void _dstedc(char compz, int n, double d[], double e[], int ldz, double z[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dstemr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, int *m, double w[], int ldz, double z[], int nzc, int isuppz[], int *tryrac, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dstein(int n, double d[], double e[], int m, double w[], int iblock[], int isplit[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
extern void _dpteqr(char compz, int n, double d[], double e[], int ldz, double z[], double work[], int *info);
#endif
extern void _dsterf(int n, double d[], double e[], int *info);
extern void _dstebz(char range, char order, int n, double vl, double vu, int il, int iu, double abstol, double d[], double e[], int *m, int *nsplit, double w[], int iblock[], int isplit[], double work[], int iwork[], int *info);
extern void _dsptrd(char uplo, int n, double ap[], double d[], double e[], double tau[], int *info);
extern void _ddisna(char job, int m, int n, double d[], double sep[], int *info);

/*
 * D4a2. Real nonsymmetric matrix _eigenvalue problems
 */
#if defined(_VLARRAY)
extern void _dgeev(char jobvl, char jobvr, int n, int lda, double a[][lda], double wr[], double wi[], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], double work[], int lwork, int *info);
extern void _dgeevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, double a[][lda], double wr[], double wi[], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], int *ilo, int *ihi, double scale[], double *abnrm, double _rconde[], double rcondv[], double work[], int lwork, int iwork[], int *info);
extern void _dgehrd(int n, int ilo, int ihi, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void _dgebal(char job, int n, int lda, double a[][lda], int *ilo, int *ihi, double scale[], int *info);
extern void _dgebak(char job, char _side, int n, int ilo, int ihi, double scale[], int m, int ldv, double v[][ldv], int *info);
extern void _dorghr(int n, int ilo, int ihi, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void _dormhr(char _side, char trans, int m, int n, int ilo, int ihi, int lda, double a[][lda], double tau[], int ldc, double c[][ldc], double work[], int lwork, int *info);
extern void _dhseqr(char job, char compz, int n, int ilo, int ihi, int ldh, double h[][ldh], double wr[], double wi[], int ldz, double z[][ldz], double work[], int lwork, int *info);
extern void _dhsein(char _side, char _eigsrc, char initv, int select[], int n, int ldh, double h[][ldh], double wr[], double wi[], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], int mm, int *m, double work[], int ifaill[], int ifailr[], int *info);
extern void _dtrevc3(char _side, char howmny, int select[], int n, int ldt, double t[][ldt], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], int mm, int *m, double work[], int lwork, int *info);
extern void _dtrexc(char compq, int n, int ldt, double t[][ldt], int ldq, double q[][ldq], int *ifst, int *ilst, double work[], int *info);
extern void _dtrsyl(char transa, char transb, int isgn, int m, int n, int lda, double a[][lda], int ldb, double b[][ldb], int ldc, double c[][ldc], double *scale, int *info);
extern void _dtrsna(char job, char howmny, int select[], int n, int ldt, double t[][ldt], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], double s[], double sep[], int mm, int *m, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dtrsen(char job, char compq, int select[], int n, int ldt, double t[][ldt], int ldq, double q[][ldq], double wr[], double wi[], int *m, double *s, double *sep, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dgees(char jobvs, char sort, int (*select)(double, double), int n, int lda, double a[][lda], int *sdim, double wr[], double wi[], int ldvs, double vs[][ldvs], double work[], int lwork, int bwork[], int *info);
extern void _dgeesx(char jobvs, char sort, int (*select)(double, double), char sense, int n, int lda, double a[][lda], int *sdim, double wr[], double wi[], int ldvs, double vs[][ldvs], double *rconde, double *rcondv, double work[], int lwork, int iwork[], int _liwork, int bwork[], int *info);
extern void _dgees_r(char jobvs, char sort, int n, int lda, double a[][lda], int *sdim, double wr[], double wi[], int ldvs, double vs[][ldvs], double work[], int lwork, int bwork[], int *info, int *irev);
extern void _dgeesx_r(char jobvs, char sort, char sense, int n, int lda, double a[][lda], int *sdim, double wr[], double wi[], int ldvs, double vs[][ldvs], double *rconde, double *rcondv, double work[], int lwork, int iwork[], int _liwork, int bwork[], int *info, int *irev);
#else
extern void _dgeev(char jobvl, char jobvr, int n, int lda, double a[], double wr[], double wi[], int ldvl, double vl[], int ldvr, double vr[], double work[], int lwork, int *info);
extern void _dgeevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, double a[], double wr[], double wi[], int ldvl, double vl[], int ldvr, double vr[], int *ilo, int *ihi, double scale[], double *abnrm, double _rconde[], double rcondv[], double work[], int lwork, int iwork[], int *info);
extern void _dgehrd(int n, int ilo, int ihi, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void _dgebal(char job, int n, int lda, double a[], int *ilo, int *ihi, double scale[], int *info);
extern void _dgebak(char job, char _side, int n, int ilo, int ihi, double scale[], int m, int ldv, double v[], int *info);
extern void _dorghr(int n, int ilo, int ihi, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void _dormhr(char _side, char trans, int m, int n, int ilo, int ihi, int lda, double a[], double tau[], int ldc, double c[], double work[], int lwork, int *info);
extern void _dhseqr(char job, char compz, int n, int ilo, int ihi, int ldh, double h[], double wr[], double wi[], int ldz, double z[], double work[], int lwork, int *info);
extern void _dhsein(char _side, char _eigsrc, char initv, int select[], int n, int ldh, double h[], double wr[], double wi[], int ldvl, double vl[], int ldvr, double vr[], int mm, int *m, double work[], int ifaill[], int ifailr[], int *info);
extern void _dtrevc3(char _side, char howmny, int select[], int n, int ldt, double t[], int ldvl, double vl[], int ldvr, double vr[], int mm, int *m, double work[], int lwork, int *info);
extern void _dtrexc(char compq, int n, int ldt, double t[], int ldq, double q[], int *ifst, int *ilst, double work[], int *info);
extern void _dtrsyl(char transa, char transb, int isgn, int m, int n, int lda, double a[], int ldb, double b[], int ldc, double c[], double *scale, int *info);
extern void _dtrsna(char job, char howmny, int select[], int n, int ldt, double t[], int ldvl, double vl[], int ldvr, double vr[], double s[], double sep[], int mm, int *m, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dtrsen(char job, char compq, int select[], int n, int ldt, double t[], int ldq, double q[], double wr[], double wi[], int *m, double *s, double *sep, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dgees(char jobvs, char sort, int (*select)(double, double), int n, int lda, double a[], int *sdim, double wr[], double wi[], int ldvs, double vs[], double work[], int lwork, int bwork[], int *info);
extern void _dgeesx(char jobvs, char sort, int (*select)(double, double), char sense, int n, int lda, double a[], int *sdim, double wr[], double wi[], int ldvs, double vs[], double *rconde, double *rcondv, double work[], int lwork, int iwork[], int _liwork, int bwork[], int *info);
extern void _dgees_r(char jobvs, char sort, int n, int lda, double a[], int *sdim, double wr[], double wi[], int ldvs, double vs[], double work[], int lwork, int bwork[], int *info, int *irev);
extern void _dgeesx_r(char jobvs, char sort, char sense, int n, int lda, double a[], int *sdim, double wr[], double wi[], int ldvs, double vs[], double *rconde, double *rcondv, double work[], int lwork, int iwork[], int _liwork, int bwork[], int *info, int *irev);
#endif

/*
 * D4a3. Hermitian matrix _eigenvalue problems
 */
#if defined(_VLARRAY)
extern void _zheev(char jobz, char uplo, int n,  int lda, doublecomplex a[][lda], double w[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zheevd(char jobz, char uplo, int n, int lda, doublecomplex a[][lda], double w[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zheevr(char jobz, char range, char uplo, int n, int lda, doublecomplex a[][lda], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], int isuppz[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zheevx(char jobz, char range, char uplo, int n, int lda, doublecomplex a[][lda], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int iwork[], int ifail[], int *info);
extern void _zhetrd(char uplo, int n, int lda, doublecomplex a[][lda], double d[], double e[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zungtr(char uplo, int n, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zunmtr(char _side, char uplo, char trans, int m, int n, int lda, doublecomplex a[][lda], doublecomplex tau[], int ldc, doublecomplex c[][ldc], doublecomplex work[], int lwork, int *info);
extern void _zupgtr(char uplo, int n, doublecomplex ap[], doublecomplex tau[], int ldq, doublecomplex q[][ldq], doublecomplex work[], int *info);
extern void _zupmtr(char _side, char uplo, char trans, int m, int n, doublecomplex ap[], doublecomplex tau[], int ldc, doublecomplex c[][ldc], doublecomplex work[], int *info);
extern void _zsteqr(char compz, int n, double d[], double e[], int ldz, doublecomplex z[][ldz], double work[], int *info);
extern void _zstedc(char compz, int n, double d[], double e[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zstemr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, int *m, double w[], int ldz, doublecomplex z[][ldz], int nzc, int isuppz[], int *tryrac, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _zstein(int n, double d[], double e[], int m, double w[], int iblock[], int isplit[], int ldz, doublecomplex z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void _zpteqr(char compz, int n, double d[], double e[], int ldz, doublecomplex z[][ldz], double work[], int *info);
extern void _zhpev(char jobz, char uplo, int n, doublecomplex ap[], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int *info);
extern void _zhpevd(char jobz, char uplo, int n, doublecomplex ap[], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zhpevx(char jobz, char range, char uplo, int n, doublecomplex ap[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void _zhbev(char jobz, char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int *info);
extern void _zhbevd(char jobz, char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zhbevx(char jobz, char range, char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], int ldq, doublecomplex q[][ldq], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void _zhbtrd(char vect, char uplo, int n, int kd, int ldab, doublecomplex ab[][ldab], double d[], double e[], int ldq, doublecomplex q[][ldq], doublecomplex work[], int *info);
#else
extern void _zheev(char jobz, char uplo, int n,  int lda, doublecomplex a[], double w[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zheevd(char jobz, char uplo, int n, int lda, doublecomplex a[], double w[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zheevr(char jobz, char range, char uplo, int n, int lda, doublecomplex a[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], int isuppz[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zheevx(char jobz, char range, char uplo, int n, int lda, doublecomplex a[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int iwork[], int ifail[], int *info);
extern void _zhetrd(char uplo, int n, int lda, doublecomplex a[], double d[], double e[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zungtr(char uplo, int n, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zunmtr(char _side, char uplo, char trans, int m, int n, int lda, doublecomplex a[], doublecomplex tau[], int ldc, doublecomplex c[], doublecomplex work[], int lwork, int *info);
extern void _zupgtr(char uplo, int n, doublecomplex ap[], doublecomplex tau[], int ldq, doublecomplex q[], doublecomplex work[], int *info);
extern void _zupmtr(char _side, char uplo, char trans, int m, int n, doublecomplex ap[], doublecomplex tau[], int ldc, doublecomplex c[], doublecomplex work[], int *info);
extern void _zsteqr(char compz, int n, double d[], double e[], int ldz, doublecomplex z[], double work[], int *info);
extern void _zstedc(char compz, int n, double d[], double e[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zstemr(char jobz, char range, int n, double d[], double e[], double vl, double vu, int il, int iu, int *m, double w[], int ldz, doublecomplex z[], int nzc, int isuppz[], int *tryrac, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _zstein(int n, double d[], double e[], int m, double w[], int iblock[], int isplit[], int ldz, doublecomplex z[], double work[], int iwork[], int ifail[], int *info);
extern void _zpteqr(char compz, int n, double d[], double e[], int ldz, doublecomplex z[], double work[], int *info);
extern void _zhpev(char jobz, char uplo, int n, doublecomplex ap[], double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int *info);
extern void _zhpevd(char jobz, char uplo, int n, doublecomplex ap[], double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zhpevx(char jobz, char range, char uplo, int n, doublecomplex ap[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void _zhbev(char jobz, char uplo, int n, int kd, int ldab, doublecomplex ab[], double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int *info);
extern void _zhbevd(char jobz, char uplo, int n, int kd, int ldab, doublecomplex ab[], double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zhbevx(char jobz, char range, char uplo, int n, int kd, int ldab, doublecomplex ab[], int ldq, doublecomplex q[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void _zhbtrd(char vect, char uplo, int n, int kd, int ldab, doublecomplex ab[], double d[], double e[], int ldq, doublecomplex q[], doublecomplex work[], int *info);
#endif
extern void _zhptrd(char uplo, int n, doublecomplex ap[], double d[], double e[], doublecomplex tau[], int *info);
/*extern void _ddisna(char job, int m, int n, double d[], double sep[], int *info);*/

/*
 * D4a4. Complex nonsymmetric matrix _eigenvalue problems
 */
#if defined(_VLARRAY)
extern void _zgeev(char jobvl, char jobvr, int n, int lda, doublecomplex a[][lda], doublecomplex w[], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgeevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, doublecomplex a[][lda], doublecomplex w[], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], int *ilo, int *ihi, double scale[], double *abnrm, double _rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgehrd(int n, int ilo, int ihi, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zgebal(char job, int n, int lda, doublecomplex a[][lda], int *ilo, int *ihi, double scale[], int *info);
extern void _zgebak(char job, char _side, int n, int ilo, int ihi, double scale[], int m, int ldv, doublecomplex v[][ldv], int *info);
extern void _zunghr(int n, int ilo, int ihi, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zunmhr(char _side, char trans, int m, int n, int ilo, int ihi, int lda, doublecomplex a[][lda], doublecomplex tau[], int ldc, doublecomplex c[][ldc], doublecomplex work[], int lwork, int *info);
extern void _zhseqr(char job, char compz, int n, int ilo, int ihi, int ldh, doublecomplex h[][ldh], doublecomplex w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, int *info);
extern void _zhsein(char _side, char _eigsrc, char initv, int select[], int n, int ldh, doublecomplex h[][ldh], doublecomplex w[], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], int mm, int *m, doublecomplex work[], double rwork[], int ifaill[], int ifailr[], int *info);
extern void _ztrevc3(char _side, char howmny, int select[], int n, int ldt, doublecomplex t[][ldt], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], int mm, int *m, doublecomplex work[], int lwork, double rwork[], int lrwork, int *info);
extern void _ztrexc(char compq, int n, int ldt, doublecomplex t[][ldt], int ldq, doublecomplex q[][ldq], int ifst, int ilst, int *info);
extern void _ztrsyl(char transa, char transb, int isgn, int m, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int ldc, doublecomplex c[][ldc], double *scale, int *info);
extern void _ztrsna(char job, char howmny, int select[], int n, int ldt, doublecomplex t[][ldt], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], double s[], double sep[], int mm, int *m, doublecomplex work[], int lwork, double rwork[], int *info);
extern void _ztrsen(char job, char compq, int select[], int n, int ldt, doublecomplex t[][ldt], int ldq, doublecomplex q[][ldq], doublecomplex w[], int *m, double *s, double *sep, doublecomplex work[], int lwork, int *info);
extern void _zgees(char jobvs, char sort, int (*select)(doublecomplex), int n, int lda, doublecomplex a[][lda], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[][ldvs], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void _zgeesx(char jobvs, char sort, int (*select)(doublecomplex), char sense, int n, int lda, doublecomplex a[][lda], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[][ldvs], double *rconde, double *rcondv, doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void _zgees_r(char jobvs, char sort, int n, int lda, doublecomplex a[][lda], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[][ldvs], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
extern void _zgeesx_r(char jobvs, char sort, char sense, int n, int lda, doublecomplex a[][lda], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[][ldvs], double *rconde, double *rcondv, doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
#else
extern void _zgeev(char jobvl, char jobvr, int n, int lda, doublecomplex a[], doublecomplex w[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgeevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, doublecomplex a[], doublecomplex w[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], int *ilo, int *ihi, double scale[], double *abnrm, double _rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgehrd(int n, int ilo, int ihi, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zgebal(char job, int n, int lda, doublecomplex a[], int *ilo, int *ihi, double scale[], int *info);
extern void _zgebak(char job, char _side, int n, int ilo, int ihi, double scale[], int m, int ldv, doublecomplex v[], int *info);
extern void _zunghr(int n, int ilo, int ihi, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zunmhr(char _side, char trans, int m, int n, int ilo, int ihi, int lda, doublecomplex a[], doublecomplex tau[], int ldc, doublecomplex c[], doublecomplex work[], int lwork, int *info);
extern void _zhseqr(char job, char compz, int n, int ilo, int ihi, int ldh, doublecomplex h[], doublecomplex w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, int *info);
extern void _zhsein(char _side, char _eigsrc, char initv, int select[], int n, int ldh, doublecomplex h[], doublecomplex w[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], int mm, int *m, doublecomplex work[], double rwork[], int ifaill[], int ifailr[], int *info);

extern void _ztrevc3(char _side, char howmny, int select[], int n, int ldt, doublecomplex t[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], int mm, int *m, doublecomplex work[], int lwork, double rwork[], int lrwork, int *info);
extern void _ztrexc(char compq, int n, int ldt, doublecomplex t[], int ldq, doublecomplex q[], int ifst, int ilst, int *info);
extern void _ztrsyl(char transa, char transb, int isgn, int m, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], int ldc, doublecomplex c[], double *scale, int *info);

extern void _ztrsna(char job, char howmny, int select[], int n, int ldt, doublecomplex t[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], double s[], double sep[], int mm, int *m, doublecomplex work[], int lwork, double rwork[], int *info);
extern void _ztrsen(char job, char compq, int select[], int n, int ldt, doublecomplex t[], int ldq, doublecomplex q[], doublecomplex w[], int *m, double *s, double *sep, doublecomplex work[], int lwork, int *info);

extern void _zgees(char jobvs, char sort, int (*select)(doublecomplex), int n, int lda, doublecomplex a[], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void _zgeesx(char jobvs, char sort, int (*select)(doublecomplex), char sense, int n, int lda, doublecomplex a[], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[], double *rconde, double *rcondv, doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void _zgees_r(char jobvs, char sort, int n, int lda, doublecomplex a[], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
extern void _zgeesx_r(char jobvs, char sort, char sense, int n, int lda, doublecomplex a[], int *sdim, doublecomplex w[], int ldvs, doublecomplex vs[], double *rconde, double *rcondv, doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
#endif

/*
 * D4b1. Real symmetric generalized matrix _eigenvalue problems
 */
#if defined(_VLARRAY)
extern void _dsygv(int itype, char jobz, char uplo, int n, int lda, double a[][lda], int ldb, double b[][ldb], double w[], double work[], int lwork, int *info);
extern void _dsygvd(int itype, char jobz, char uplo, int n, int lda, double a[][lda], int ldb, double b[][ldb], double w[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dsygvx(int itype, char jobz, char range, char uplo, int n, int lda, double a[][lda], int ldb, double b[][ldb], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int ifail[], int *info);
extern void _dspgv(int itype, char jobz, char uplo, int n, double ap[], double bp[], double w[], int ldz, double z[][ldz], double work[], int *info);
extern void _dspgvd(int itype, char jobz, char uplo, int n, double ap[], double bp[], double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dspgvx(int itype, char jobz, char range, char uplo, int n, double ap[], double bp[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
extern void _dsbgv(char jobz, char uplo, int n, int ka, int kb, int ldab, double ab[][ldab], int ldbb, double bb[][ldbb], double w[], int ldz, double z[][ldz], double work[], int *info);
extern void _dsbgvd(char jobz, char uplo, int n, int ka, int kb, int ldab, double ab[][ldab], int ldbb, double bb[][ldbb], double w[], int ldz, double z[][ldz], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dsbgvx(char jobz, char range, char uplo, int n, int ka, int kb, int ldab, double ab[][ldab], int ldbb, double bb[][ldbb], int ldq, double q[][ldq], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[][ldz], double work[], int iwork[], int ifail[], int *info);
#else
extern void _dsygv(int itype, char jobz, char uplo, int n, int lda, double a[], int ldb, double b[], double w[], double work[], int lwork, int *info);
extern void _dsygvd(int itype, char jobz, char uplo, int n, int lda, double a[], int ldb, double b[], double w[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dsygvx(int itype, char jobz, char range, char uplo, int n, int lda, double a[], int ldb, double b[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int lwork, int iwork[], int ifail[], int *info);
extern void _dspgv(int itype, char jobz, char uplo, int n, double ap[], double bp[], double w[], int ldz, double z[], double work[], int *info);
extern void _dspgvd(int itype, char jobz, char uplo, int n, double ap[], double bp[], double w[], int ldz, double z[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dspgvx(int itype, char jobz, char range, char uplo, int n, double ap[], double bp[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
extern void _dsbgv(char jobz, char uplo, int n, int ka, int kb, int ldab, double ab[], int ldbb, double bb[], double w[], int ldz, double z[], double work[], int *info);
extern void _dsbgvd(char jobz, char uplo, int n, int ka, int kb, int ldab, double ab[], int ldbb, double bb[], double w[], int ldz, double z[], double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _dsbgvx(char jobz, char range, char uplo, int n, int ka, int kb, int ldab, double ab[], int ldbb, double bb[], int ldq, double q[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, double z[], double work[], int iwork[], int ifail[], int *info);
#endif

/*
 * D4b2. Real generalized matrix _eigenvalue problems
 */
#if defined(_VLARRAY)
extern void _dggev(char jobvl, char jobvr, int n, int lda, double a[][lda], int ldb, double b[][ldb], double alphar[], double alphai[], double _beta[], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], double work[], int lwork, int *info);
extern void _dggevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, double a[][lda], int ldb, double b[][ldb], double alphar[], double alphai[], double _beta[], int ldvl, double vl[][ldvl], int ldvr, double vr[][ldvr], int *ilo, int *ihi, double lscale[], double rscale[], double *abnrm, double *bbnrm, double _rconde[], double rcondv[], double work[], int lwork, int iwork[], int bwork[], int *info);
extern void _dgges(char jobvsl, char jobvsr, char sort, int (*selctg)(double, double, double), int n, int lda, double a[][lda], int ldb, double b[][ldb], int *sdim, double alphar[], double alphai[], double _beta[], int ldvsl, double vsl[][ldvsl], int ldvsr, double vsr[][ldvsr], double work[], int lwork, int bwork[], int *info);
extern void _dggesx(char jobvsl, char jobvsr, char sort, int (*selctg)(double, double, double), char sense, int n, int lda, double a[][lda], int ldb, double b[][ldb], int *sdim, double alphar[], double alphai[], double _beta[], int ldvsl, double vsl[][ldvsl], int ldvsr, double vsr[][ldvsr], double _rconde[], double rcondv[], double work[], int lwork, int iwork[], int _liwork, int bwork[], int *info);
extern void _dgges_r(char jobvsl, char jobvsr, char sort, int n, int lda, double a[][lda], int ldb, double b[][ldb], int *sdim, double alphar[], double alphai[], double _beta[], int ldvsl, double vsl[][ldvsl], int ldvsr, double vsr[][ldvsr], double work[], int lwork, int bwork[], int *info, int *irev);
extern void _dggesx_r(char jobvsl, char jobvsr, char sort, char sense, int n, int lda, double a[][lda], int ldb, double b[][ldb], int *sdim, double alphar[], double alphai[], double _beta[], int ldvsl, double vsl[][ldvsl], int ldvsr, double vsr[][ldvsr], double _rconde[], double rcondv[], double work[], int lwork, int iwork[], int _liwork, int bwork[], int *info, int *irev);
#else
extern void _dggev(char jobvl, char jobvr, int n, int lda, double a[], int ldb, double b[], double alphar[], double alphai[], double _beta[], int ldvl, double vl[], int ldvr, double vr[], double work[], int lwork, int *info);
extern void _dggevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, double a[], int ldb, double b[], double alphar[], double alphai[], double _beta[], int ldvl, double vl[], int ldvr, double vr[], int *ilo, int *ihi, double lscale[], double rscale[], double *abnrm, double *bbnrm, double _rconde[], double rcondv[], double work[], int lwork, int iwork[], int bwork[], int *info);
extern void _dgges(char jobvsl, char jobvsr, char sort, int (*selctg)(double, double, double), int n, int lda, double a[], int ldb, double b[], int *sdim, double alphar[], double alphai[], double _beta[], int ldvsl, double vsl[], int ldvsr, double vsr[], double work[], int lwork, int bwork[], int *info);
extern void _dggesx(char jobvsl, char jobvsr, char sort, int (*selctg)(double, double, double), char sense, int n, int lda, double a[], int ldb, double b[], int *sdim, double alphar[], double alphai[], double _beta[], int ldvsl, double vsl[], int ldvsr, double vsr[], double _rconde[], double rcondv[], double work[], int lwork, int iwork[], int _liwork, int bwork[], int *info);
extern void _dgges_r(char jobvsl, char jobvsr, char sort, int n, int lda, double a[], int ldb, double b[], int *sdim, double alphar[], double alphai[], double _beta[], int ldvsl, double vsl[], int ldvsr, double vsr[], double work[], int lwork, int bwork[], int *info, int *irev);
extern void _dggesx_r(char jobvsl, char jobvsr, char sort, char sense, int n, int lda, double a[], int ldb, double b[], int *sdim, double alphar[], double alphai[], double _beta[], int ldvsl, double vsl[], int ldvsr, double vsr[], double _rconde[], double rcondv[], double work[], int lwork, int iwork[], int _liwork, int bwork[], int *info, int *irev);
#endif

/*
 * D4b3. Complex Hermitian generalized matrix _eigenvalue problems
 */
#if defined(_VLARRAY)
extern void _zhegv(int itype, char jobz, char uplo, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double w[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zhegvd(int itype, char jobz, char uplo, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double w[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zhegvx(int itype, char jobz, char range, char uplo, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int iwork[], int ifail[], int *info);
extern void _zhpgv(int itype, char jobz, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int *info);
extern void _zhpgvd(int itype, char jobz, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zhpgvx(int itype, char jobz, char range, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void _zhbgv(char jobz, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[][ldab], int ldbb, doublecomplex bb[][ldbb], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int *info);
extern void _zhbgvd(char jobz, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[][ldab], int ldbb, doublecomplex bb[][ldbb], double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zhbgvx(char jobz, char range, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[][ldab], int ldbb, doublecomplex bb[][ldbb], int ldq, doublecomplex q[][ldq], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[][ldz], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
#else
extern void _zhegv(int itype, char jobz, char uplo, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], double w[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zhegvd(int itype, char jobz, char uplo, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], double w[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zhegvx(int itype, char jobz, char range, char uplo, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int iwork[], int ifail[], int *info);
extern void _zhpgv(int itype, char jobz, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int *info);
extern void _zhpgvd(int itype, char jobz, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zhpgvx(int itype, char jobz, char range, char uplo, int n, doublecomplex ap[], doublecomplex bp[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
extern void _zhbgv(char jobz, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[], int ldbb, doublecomplex bb[], double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int *info);
extern void _zhbgvd(char jobz, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[], int ldbb, doublecomplex bb[], double w[], int ldz, doublecomplex z[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zhbgvx(char jobz, char range, char uplo, int n, int ka, int kb, int ldab, doublecomplex ab[], int ldbb, doublecomplex bb[], int ldq, doublecomplex q[], double vl, double vu, int il, int iu, double abstol, int *m, double w[], int ldz, doublecomplex z[], doublecomplex work[], double rwork[], int iwork[], int ifail[], int *info);
#endif

/*
 * D4b4. Complex generalized matrix _eigenvalue problems
 */
#if defined(_VLARRAY)
extern void _zggev(char jobvl, char jobvr, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex alpha[], doublecomplex _beta[], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zggevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex alpha[], doublecomplex _beta[], int ldvl, doublecomplex vl[][ldvl], int ldvr, doublecomplex vr[][ldvr], int *ilo, int *ihi, double lscale[], double rscale[], double *abnrm, double *bbnrm, double _rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int bwork[], int *info);
extern void _zgges(char jobvsl, char jobvsr, char sort, int (*selctg)(doublecomplex, doublecomplex), int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *sdim, doublecomplex alpha[], doublecomplex _beta[], int ldvsl, doublecomplex vsl[][ldvsl], int ldvsr, doublecomplex vsr[][ldvsr], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void _zggesx(char jobvsl, char jobvsr, char sort, int (*selctg)(doublecomplex, doublecomplex), char sense, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *sdim, doublecomplex alpha[], doublecomplex _beta[], int ldvsl, doublecomplex vsl[][ldvsl], int ldvsr, doublecomplex vsr[][ldvsr], double _rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int _liwork, int bwork[], int *info);
extern void _zgges_r(char jobvsl, char jobvsr, char sort, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *sdim, doublecomplex alpha[], doublecomplex _beta[], int ldvsl, doublecomplex vsl[][ldvsl], int ldvsr, doublecomplex vsr[][ldvsr], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
extern void _zggesx_r(char jobvsl, char jobvsr, char sort, char sense, int n, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int *sdim, doublecomplex alpha[], doublecomplex _beta[], int ldvsl, doublecomplex vsl[][ldvsl], int ldvsr, doublecomplex vsr[][ldvsr], double _rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int _liwork, int bwork[], int *info, int *irev);
#else
extern void _zggev(char jobvl, char jobvr, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex alpha[], doublecomplex _beta[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zggevx(char balanc, char jobvl, char jobvr, char sense, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex alpha[], doublecomplex _beta[], int ldvl, doublecomplex vl[], int ldvr, doublecomplex vr[], int *ilo, int *ihi, double lscale[], double rscale[], double *abnrm, double *bbnrm, double _rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int bwork[], int *info);
extern void _zgges(char jobvsl, char jobvsr, char sort, int (*selctg)(doublecomplex, doublecomplex), int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *sdim, doublecomplex alpha[], doublecomplex _beta[], int ldvsl, doublecomplex vsl[], int ldvsr, doublecomplex vsr[], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info);
extern void _zggesx(char jobvsl, char jobvsr, char sort, int (*selctg)(doublecomplex, doublecomplex), char sense, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *sdim, doublecomplex alpha[], doublecomplex _beta[], int ldvsl, doublecomplex vsl[], int ldvsr, doublecomplex vsr[], double _rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int _liwork, int bwork[], int *info);
extern void _zgges_r(char jobvsl, char jobvsr, char sort, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *sdim, doublecomplex alpha[], doublecomplex _beta[], int ldvsl, doublecomplex vsl[], int ldvsr, doublecomplex vsr[], doublecomplex work[], int lwork, double rwork[], int bwork[], int *info, int *irev);
extern void _zggesx_r(char jobvsl, char jobvsr, char sort, char sense, int n, int lda, doublecomplex a[], int ldb, doublecomplex b[], int *sdim, doublecomplex alpha[], doublecomplex _beta[], int ldvsl, doublecomplex vsl[], int ldvsr, doublecomplex vsr[], double _rconde[], double rcondv[], doublecomplex work[], int lwork, double rwork[], int iwork[], int _liwork, int bwork[], int *info, int *irev);
#endif

/*
 * D5. QR factorization
 */
#if defined(_VLARRAY)
extern void _dgeqp3(int m, int n, int lda, double a[][lda], int jpvt[], double tau[], double work[], int lwork, int *info);
extern void _dgeqrf(int m, int n, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void _dorgqr(int m, int n, int k, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void _dormqr(char _side, char trans, int m, int n, int k, int lda, double a[][lda], double tau[], int ldc, double c[][ldc], double work[], int lwork, int *info);
extern void _dgelqf(int m, int n, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void _dorglq(int m, int n, int k, int lda, double a[][lda], double tau[], double work[], int lwork, int *info);
extern void _dormlq(char _side, char trans, int m, int n, int k, int lda, double a[][lda], double tau[], int ldc, double c[][ldc], double work[], int lwork, int *info);
#else
extern void _dgeqp3(int m, int n, int lda, double a[], int jpvt[], double tau[], double work[], int lwork, int *info);
extern void _dgeqrf(int m, int n, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void _dorgqr(int m, int n, int k, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void _dormqr(char _side, char trans, int m, int n, int k, int lda, double a[], double tau[], int ldc, double c[], double work[], int lwork, int *info);
extern void _dgelqf(int m, int n, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void _dorglq(int m, int n, int k, int lda, double a[], double tau[], double work[], int lwork, int *info);
extern void _dormlq(char _side, char trans, int m, int n, int k, int lda, double a[], double tau[], int ldc, double c[], double work[], int lwork, int *info);
#endif

/*
 * D5-2. QR factorization (complex matrices)
 */
#if defined(_VLARRAY)
extern void _zgeqp3(int m, int n, int lda, doublecomplex a[][lda], int jpvt[], doublecomplex tau[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgeqrf(int m, int n, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zungqr(int m, int n, int k, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zunmqr(char _side, char trans, int m, int n, int k, int lda, doublecomplex a[][lda], doublecomplex tau[], int ldc, doublecomplex c[][ldc], doublecomplex work[], int lwork, int *info);
extern void _zgelqf(int m, int n, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zunglq(int m, int n, int k, int lda, doublecomplex a[][lda], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zunmlq(char _side, char trans, int m, int n, int k, int lda, doublecomplex a[][lda], doublecomplex tau[], int ldc, doublecomplex c[][ldc], doublecomplex work[], int lwork, int *info);
#else
extern void _zgeqp3(int m, int n, int lda, doublecomplex a[], int jpvt[], doublecomplex tau[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgeqrf(int m, int n, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zungqr(int m, int n, int k, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zunmqr(char _side, char trans, int m, int n, int k, int lda, doublecomplex a[], doublecomplex tau[], int ldc, doublecomplex c[], doublecomplex work[], int lwork, int *info);
extern void _zgelqf(int m, int n, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zunglq(int m, int n, int k, int lda, doublecomplex a[], doublecomplex tau[], doublecomplex work[], int lwork, int *info);
extern void _zunmlq(char _side, char trans, int m, int n, int k, int lda, doublecomplex a[], doublecomplex tau[], int ldc, doublecomplex c[], doublecomplex work[], int lwork, int *info);
#endif

/*
 * D6. Singular value decomposition (SVD)
 */
#if defined(_VLARRAY)
extern void _dgesvd(char jobu, char jobvt, int m, int n, int lda, double a[][lda], double s[], int ldu, double u[][ldu], int ldvt, double vt[][ldvt], double work[], int lwork, int *info);
extern void _dgesvdx(char jobu, char jobvt, char range, int m, int n, int lda, double a[][lda], double vl, double vu, int il, int iu, int *ns, double s[], int ldu, double u[][ldu], int ldvt, double vt[][ldvt], double work[], int lwork, int iwork[], int *info);
extern void _dgesdd(char jobz, int m, int n, int lda, double a[][lda], double s[], int ldu, double u[][ldu], int ldvt, double vt[][ldvt], double work[], int lwork, int iwork[], int *info);
extern void _dgesvdq(char joba, char jobp, char jobr, char jobu, char jobv, int m, int n, int lda, double a[][lda], double s[], int ldu, double u[][ldu], int ldv, double v[][ldv], int *numrank, double work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _dgejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, int m, int n, int lda, double a[][lda], double sva[], int ldu, double u[][ldu], int ldv, double v[][ldv], double work[], int lwork, int iwork[], int *info);
extern void _dggsvd3(char jobu, char jobv, char jobq, int m, int n, int p, int *k, int *l, int lda, double a[][lda], int ldb, double b[][ldb], double alpha[], double _beta[], int ldu, double u[][ldu], int ldv, double v[][ldv], int ldq, double q[][ldq], double work[], int lwork, int iwork[], int *info);
#else
extern void _dgesvd(char jobu, char jobvt, int m, int n, int lda, double a[], double s[], int ldu, double u[], int ldvt, double vt[], double work[], int lwork, int *info);
extern void _dgesvdx(char jobu, char jobvt, char range, int m, int n, int lda, double a[], double vl, double vu, int il, int iu, int *ns, double s[], int ldu, double u[], int ldvt, double vt[], double work[], int lwork, int iwork[], int *info);
extern void _dgesdd(char jobz, int m, int n, int lda, double a[], double s[], int ldu, double u[], int ldvt, double vt[], double work[], int lwork, int iwork[], int *info);
extern void _dgesvdq(char joba, char jobp, char jobr, char jobu, char jobv, int m, int n, int lda, double a[], double s[], int ldu, double u[], int ldv, double v[], int *numrank, double work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _dgejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, int m, int n, int lda, double a[], double sva[], int ldu, double u[], int ldv, double v[], double work[], int lwork, int iwork[], int *info);
extern void _dggsvd3(char jobu, char jobv, char jobq, int m, int n, int p, int *k, int *l, int lda, double a[], int ldb, double b[], double alpha[], double _beta[], int ldu, double u[], int ldv, double v[], int ldq, double q[], double work[], int lwork, int iwork[], int *info);
#endif

/*
 * D6-2. Singular value decomposition (SVD) (complex matrices)
 */
#if defined(_VLARRAY)
extern void _zgesvd(char jobu, char jobvt, int m, int n, int lda, doublecomplex a[][lda], double s[], int ldu, doublecomplex u[][ldu], int ldvt, doublecomplex vt[][ldvt], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgesvdx(char jobu, char jobvt, char range, int m, int n, int lda, doublecomplex a[][lda], double vl, double vu, int il, int iu, int *ns, double s[], int ldu, doublecomplex u[][ldu], int ldvt, doublecomplex vt[][ldvt], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void _zgesdd(char jobz, int m, int n, int lda, doublecomplex a[][lda], double s[], int ldu, doublecomplex u[][ldu], int ldvt, doublecomplex vt[][ldvt], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void _zgesvdq(char joba, char jobp, char jobr, char jobu, char jobv, int m, int n, int lda, doublecomplex a[][lda], double s[], int ldu, doublecomplex u[][ldu], int ldv, doublecomplex v[][ldv], int *numrank, doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zgejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, int m, int n, int lda, doublecomplex a[][lda], double sva[], int ldu, doublecomplex u[][ldu], int ldv, doublecomplex v[][ldv], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int *info);
extern void _zggsvd3(char jobu, char jobv, char jobq, int m, int n, int p, int *k, int *l, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double alpha[], double _beta[], int ldu, doublecomplex u[][ldu], int ldv, doublecomplex v[][ldv], int ldq, doublecomplex q[][ldq], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
#else
extern void _zgesvd(char jobu, char jobvt, int m, int n, int lda, doublecomplex a[], double s[], int ldu, doublecomplex u[], int ldvt, doublecomplex vt[], doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgesvdx(char jobu, char jobvt, char range, int m, int n, int lda, doublecomplex a[], double vl, double vu, int il, int iu, int *ns, double s[], int ldu, doublecomplex u[], int ldvt, doublecomplex vt[], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void _zgesdd(char jobz, int m, int n, int lda, doublecomplex a[], double s[], int ldu, doublecomplex u[], int ldvt, doublecomplex vt[], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void _zgesvdq(char joba, char jobp, char jobr, char jobu, char jobv, int m, int n, int lda, doublecomplex a[], double s[], int ldu, doublecomplex u[], int ldv, doublecomplex v[], int *numrank, doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int _liwork, int *info);
extern void _zgejsv(char joba, char jobu, char jobv, char jobr, char jobt, char jobp, int m, int n, int lda, doublecomplex a[], double sva[], int ldu, doublecomplex u[], int ldv, doublecomplex v[], doublecomplex work[], int lwork, double rwork[], int lrwork, int iwork[], int *info);
extern void _zggsvd3(char jobu, char jobv, char jobq, int m, int n, int p, int *k, int *l, int lda, doublecomplex a[], int ldb, doublecomplex b[], double alpha[], double _beta[], int ldu, doublecomplex u[], int ldv, doublecomplex v[], int ldq, doublecomplex q[], doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
#endif

/*
 * D9a. Singular, overdetermined or underdetermined systems of _linear equations without constraints
 */
#if defined(_VLARRAY)
extern void _dgels(char trans, int m, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], double work[], int lwork, int *info);
extern void _dgelsy(int m, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], int jpvt[], double _rcond, int *rank, double work[], int lwork, int *info);
extern void _dgelss(int m, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], double s[], double _rcond, int *rank, double work[], int lwork, int *info);
extern void _dgetsls(char trans, int m, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], double work[], int lwork, int *info);
extern void _dgelsd(int m, int n, int nrhs, int lda, double a[][lda], int ldb, double b[][ldb], double s[], double _rcond, int *rank, double work[], int lwork, int iwork[], int *info);
extern void _dgecov(int job, int n, int lda, double a[][lda], double _ci[], int *info);
extern void _dgecovy(int job, int n, int lda, double a[][lda], int ipiv[], double _ci[], int iwork[], int *info);
extern void _dgecovs(int job, int n, int lda, double a[][lda], double s[], double _ci[], double work[], int *info);
#else
extern void _dgels(char trans, int m, int n, int nrhs, int lda, double a[], int ldb, double b[], double work[], int lwork, int *info);
extern void _dgelsy(int m, int n, int nrhs, int lda, double a[], int ldb, double b[], int jpvt[], double _rcond, int *rank, double work[], int lwork, int *info);
extern void _dgelss(int m, int n, int nrhs, int lda, double a[], int ldb, double b[], double s[], double _rcond, int *rank, double work[], int lwork, int *info);
extern void _dgetsls(char trans, int m, int n, int nrhs, int lda, double a[], int ldb, double b[], double work[], int lwork, int *info);
extern void _dgelsd(int m, int n, int nrhs, int lda, double a[], int ldb, double b[], double s[], double _rcond, int *rank, double work[], int lwork, int iwork[], int *info);
extern void _dgecov(int job, int n, int lda, double a[], double _ci[], int *info);
extern void _dgecovy(int job, int n, int lda, double a[], int ipiv[], double _ci[], int iwork[], int *info);
extern void _dgecovs(int job, int n, int lda, double a[], double s[], double _ci[], double work[], int *info);
#endif

/*
 * D9a-2. Singular, overdetermined or underdetermined systems of _linear equations without constraints (complex matrices)
 */
#if defined(_VLARRAY)
extern void _zgels(char trans, int m, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex work[], int lwork, int *info);
extern void _zgelsy(int m, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], int jpvt[], double _rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgelss(int m, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double s[], double _rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgetsls(char trans, int m, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex work[], int lwork, int *info);
extern void _zgelsd(int m, int n, int nrhs, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], double s[], double _rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void _zgecov(int job, int n, int lda, doublecomplex a[][lda], doublecomplex _ci[], int *info);
extern void _zgecovy(int job, int n, int lda, doublecomplex a[][lda], int ipiv[], doublecomplex _ci[], int iwork[], int *info);
extern void _zgecovs(int job, int n, int lda, doublecomplex a[][lda], double s[], doublecomplex _ci[], doublecomplex work[], int *info);
#else
extern void _zgels(char trans, int m, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex work[], int lwork, int *info);
extern void _zgelsy(int m, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], int jpvt[], double _rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgelss(int m, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], double s[], double _rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int *info);
extern void _zgetsls(char trans, int m, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex work[], int lwork, int *info);
extern void _zgelsd(int m, int n, int nrhs, int lda, doublecomplex a[], int ldb, doublecomplex b[], double s[], double _rcond, int *rank, doublecomplex work[], int lwork, double rwork[], int iwork[], int *info);
extern void _zgecov(int job, int n, int lda, doublecomplex a[], doublecomplex _ci[], int *info);
extern void _zgecovy(int job, int n, int lda, doublecomplex a[], int ipiv[], doublecomplex _ci[], int iwork[], int *info);
extern void _zgecovs(int job, int n, int lda, doublecomplex a[], double s[], doublecomplex _ci[], doublecomplex work[], int *info);
#endif

/*
 * D9b. Singular, overdetermined or underdetermined systems of _linear equations with constraints
 */
#if defined(_VLARRAY)
extern void _dgglse(int m, int n, int p, int lda, double a[][lda], int ldb, double b[][ldb], double c[], double d[], double x[], double work[], int lwork, int *info);
extern void _dggglm(int n, int m, int p, int lda, double a[][lda], int ldb, double b[][ldb], double d[], double x[], double y[], double work[], int lwork, int *info);
#else
extern void _dgglse(int m, int n, int p, int lda, double a[], int ldb, double b[], double c[], double d[], double x[], double work[], int lwork, int *info);
extern void _dggglm(int n, int m, int p, int lda, double a[], int ldb, double b[], double d[], double x[], double y[], double work[], int lwork, int *info);
#endif

/*
 * D9b-2. Singular, overdetermined or underdetermined systems of complex _linear equations with constraints (complex matrices)
 */
#if defined(_VLARRAY)
extern void _zgglse(int m, int n, int p, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex c[], doublecomplex d[], doublecomplex x[], doublecomplex work[], int lwork, int *info);
extern void _zggglm(int n, int m, int p, int lda, doublecomplex a[][lda], int ldb, doublecomplex b[][ldb], doublecomplex d[], doublecomplex x[], doublecomplex y[], doublecomplex work[], int lwork, int *info);
#else
extern void _zgglse(int m, int n, int p, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex c[], doublecomplex d[], doublecomplex x[], doublecomplex work[], int lwork, int *info);
extern void _zggglm(int n, int m, int p, int lda, doublecomplex a[], int ldb, doublecomplex b[], doublecomplex d[], doublecomplex x[], doublecomplex y[], doublecomplex work[], int lwork, int *info);
#endif

/*
 * E. Interpolation
 */
extern void _polint(int n, double x[], double y[], double c[], int *info);
extern void _polyvl(int nder, double xx, double *yfit, double yp[], int n, double x[], double c[], double work[], int *info);
extern void _polcof(double xx, int n, double x[], double c[], double d[], double work[]);
extern double _fitlag(double x, int n, double a[], double f[], int m, double eps, int *info);
extern void _pchim(int n, double x[], double f[], double d[], int incfd, int *info);
extern void _pchic(int ic[], double vc[], double sw, int n, double x[], double f[], double d[], int incfd, double work[], int lwork, int *info);
extern void _pchse(int n, double x[], double f[], double d[], int incfd, double work[], int lwork, int *info);
extern void _pchsp(int ic[], double vc[], int n, double x[], double f[], double d[], int incfd, double work[], int lwork, int *info);
extern void _pchfe(int n, double x[], double f[], double d[], int incfd, int skip, int ne, double xe[], double fe[], int *info);
extern void _pchfd(int n, double x[], double f[], double d[], int incfd, int skip, int ne, double xe[], double fe[], double de[], int *info);
extern void _chfev(double x1, double x2, double f1, double f2, double d1, double d2, int ne, double xe[], double fe[], int next[], int *info);
extern void _chfdv(double x1, double x2, double f1, double f2, double d1, double d2, int ne, double xe[], double fe[], double de[], int next[], int *info);
extern void _pchbs(int n, double x[], double f[], double d[], int incfd, int knotyp, int *nknots, double t[], double bcoef[], int *ndim, int *kord, int *info);
extern void _pchcm(int n, double x[], double f[], double d[], int incfd, int skip, int ismon[], int *info);

extern void _bint4(double x[], double y[], int ndata, int ibcl, int ibcr, double fbcl, double fbcr, int kntopt, double t[], double bcoef[], int *n, int *k, double work[], int *info);
extern void _bintk(double x[], double y[], double t[], int n, int k, double bcoef[], double q[], double work[], int *info);
extern double _bvalue(double t[], double a[], int n, int k, int ideriv, double x, int *inbv, double work[], int *info);
#if defined(_VLARRAY)
extern double _ppvalu(int ldc, double c[][ldc], double xi[], int lxi, int k, int ideriv, double x, int *inppv, int *info);
extern void _bsplpp(double t[], double a[], int n, int k, int ldc, double c[][ldc], double xi[], int *lxi, double work[], int *info);
extern void _bsplvd(double t[], int k, int nderiv, double x, int ileft, int ldvnikx, double vnikx[][ldvnikx], double work[], int *info);
extern void _banfac(int n, int kl, int ku, int ldab, double ab[][ldab], int *info);
extern void _banslv(int n, int kl, int ku, int ldab, double ab[][ldab], double b[], int *info);
#else
extern double _ppvalu(int ldc, double c[], double xi[], int lxi, int k, int ideriv, double x, int *inppv, int *info);
extern void _bsplpp(double t[], double a[], int n, int k, int ldc, double c[], double xi[], int *lxi, double work[], int *info);
extern void _bsplvd(double t[], int k, int nderiv, double x, int ileft, int ldvnikx, double vnikx[], double work[], int *info);
extern void _banfac(int n, int kl, int ku, int ldab, double ab[], int *info);
extern void _banslv(int n, int kl, int ku, int ldab, double ab[], double b[], int *info);
#endif
extern void _bsplvn(double t[], int jhigh, int k, int index, double x, int ileft, double vnikx[], double work[], int *iwork, int *info);
extern void _bspldr(double t[], double a[], int n, int k, int nderiv, double ad[], int *info);
extern void _bsplev(double t[], double ad[], int n, int k, int nderiv, double x, int *inev, double svalue[], double work[], int *info);
extern void _interv(double xt[], int lxt, double x, int *ilo, int *ileft, int *info);

extern double _pchia(int n, double x[], double f[], double d[], int incfd, int skip, double a, double b, int *info);
extern double _pchid(int n, double x[], double f[], double d[], int incfd, int skip, int ia, int ib, int *info);
extern void _bsqad(double t[], double bcoef[], int n, int k, double x1, double x2, double *bquad, double work[], int *info);
extern void _bfqad(double (*f)(double), double t[], double bcoef[], int n, int k, int id, double x1, double x2, double tol, double *quad, double work[], int *info);
extern void _bfqad_r(double t[], double bcoef[], int n, int k, int id, double x1, double x2, double tol, double *quad, double work[], int *info, double *xx, double yy, int *irev);
#if defined(_VLARRAY)
extern void _ppqad(int ldc, double c[][ldc], double xi[], int lxi, int k, double x1, double x2, double *pquad, int *info);
extern void _pfqad(double (*f)(double), int ldc, double c[][ldc], double xi[], int lxi, int k, int id, double x1, double x2, double tol, double *quad, int *info);
extern void _pfqad_r(int ldc, double c[][ldc], double xi[], int lxi, int k, int id, double x1, double x2, double tol, double *quad, int *info, double *xx, double yy, int *irev);
#else
extern void _ppqad(int ldc, double c[], double xi[], int lxi, int k, double x1, double x2, double *pquad, int *info);
extern void _pfqad(double (*f)(double), int ldc, double c[], double xi[], int lxi, int k, int id, double x1, double x2, double tol, double *quad, int *info);
extern void _pfqad_r(int ldc, double c[], double xi[], int lxi, int k, int id, double x1, double x2, double tol, double *quad, int *info, double *xx, double yy, int *irev);
#endif

/*
 * F1a. Roots of polynomials
 */
extern void _cpzero(int n, doublecomplex a[], doublecomplex r[], int iflag, int maxiter, int *iter, double s[], doublecomplex work[], int *info);
extern void _rpzero(int n, double a[], doublecomplex r[], int iflag, int maxiter, int *iter, double s[], doublecomplex work[], int *info);
extern void _rpzero2(int n, double a[], double rr[], double ri[], int iflag, int maxiter, int *iter, double s[], double work[], int *info);
extern void _cpqr79(int n, doublecomplex a[], doublecomplex r[], doublecomplex work[], int lwork, int *info);
extern void _rpqr79(int n, double a[], doublecomplex r[], double work[], int lwork, int *info);
extern void _dka(int n, doublecomplex a[], doublecomplex r[], int maxiter, int *iter, doublecomplex work[], double rwork[], int *info);

/*
 * F1b. Solution of _single general nonlinear equation
 */
extern void _dfzero(double (*f)(double), double *b, double *c, double r, double re, double ae, int *info);
extern void _dfzero_r(double *b, double *c, double r, double re, double ae, int *info, double *xx, double yy, int *irev);

/*
 * F2. Solution of a system of nonlinear equations
 */
#if defined(_VLARRAY)
extern void _hybrj(void (*fcn)(int, double *, double *, int, double (*)[*], int *), int n, double x[], double fvec[], double xtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, double work[], int lwork, int *info);
extern void _hybrj1(void (*fcn)(int, double *, double *, int, double (*)[*], int *), int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info);
#else
extern void _hybrj(void (*fcn)(int, double *, double *, int, double *, int *), int n, double x[], double fvec[], double xtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, double work[], int lwork, int *info);
extern void _hybrj1(void (*fcn)(int, double *, double *, int, double *, int *), int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info);
#endif
extern void _hybrd(void (*fcn)(int, double *, double *, int *), int n, double x[], double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, double work[], int lwork, int *info);
extern void _hybrd1(void (*fcn)(int, double *, double *, int *), int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info);
#if defined(_VLARRAY)
extern void _chkder(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double xp[], double fvecp[], int mode, double err[], int *info);
#else
extern void _chkder(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double xp[], double fvecp[], int mode, double err[], int *info);
#endif
extern void _sos(void (*fnc)(int, double *, int, double *), int neq, double x[], double rtolx, double atolx, double tolf, int nprint, int maxiter, int *iter, double work[], int lwork, int iwork[], int _liwork, int *info);

#if defined(_VLARRAY)
extern void _hybrj_r(int n, double x[], double fvec[], double xtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, double work[], int lwork, int *info, double xx[], double yy[], int ldyypd, double yypd[][ldyypd], int *irev);
extern void _hybrj1_r(int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info, double xx[], double yy[], int ldyypd, double yypd[][ldyypd], int *irev);
#else
extern void _hybrj_r(int n, double x[], double fvec[], double xtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, double work[], int lwork, int *info, double xx[], double yy[], int ldyypd, double yypd[], int *irev);
extern void _hybrj1_r(int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info, double xx[], double yy[], int ldyypd, double yypd[], int *irev);
#endif
extern void _hybrd_r(int n, double x[], double fvec[], double xtol, int maxfev, int ml, int mu, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void _hybrd1_r(int n, double x[], double fvec[], double xtol, double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void _sos_r(int neq, double x[], double rtolx, double atolx, double tolf, int nprint, int maxiter, int *iter, double work[], int lwork, int iwork[], int _liwork, int *info, int *k, double xx[], double *yy, int *irev);

/*
 * G1a. Unconstrained optimization of a general univariate function
 */
extern double _dfmin(double ax, double bx, double (*f)(double), double tol);
extern void _dfmin_r(double ax, double bx, double tol, double *xx, double yy, int *irev);

/*
 * G1b. Unconstrained optimization of a general multivariate function
 */
#if defined(_VLARRAY)
extern void _optif9(int n, double x[], void (*fcn)(int, double *, double *), void (*d1fcn)(int, double *, double *), void (*d2fcn)(int, double *, int, double [][*]), double typsiz[], double fscale, int method, int iexp, int ndigit, int maxiter, int iagflg, int iahflg, double dlt, double gradtl, double stepmx, double steptl, void (*result)(int, double *, double, double *, int, double [][*], double *, int, int), int iresult, double xpls[], double *fpls, double gpls[], int *iter, double work[], int lwork, int *info);
#else
extern void _optif9(int n, double x[], void (*fcn)(int, double *, double *), void (*d1fcn)(int, double *, double *), void (*d2fcn)(int, double *, int, double *), double typsiz[], double fscale, int method, int iexp, int ndigit, int maxiter, int iagflg, int iahflg, double dlt, double gradtl, double stepmx, double steptl, void (*result)(int, double *, double, double *, int, double *, double *, int, int), int iresult, double xpls[], double *fpls, double gpls[], int *iter, double work[], int lwork, int *info);
#endif
extern void _optif0(int n, double x[], void (*fcn)(int, double *, double *), double xpls[], double *fpls, double work[], int lwork, int *info);
extern void _mng(int n, double x[], void (*calcf)(int, double *, int *, double *), void (*calcg)(int, double *, int *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);
extern void _mnf(int n, double x[], void (*calcf)(int, double *, int *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);
extern void _mnh(int n, double x[], void (*calcf)(int, double *, int *, double *), void (*calcgh)(int, double *, int *, double *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);
extern void _ivset(int alg, double v[], int lv, int iv[], int _liv);
extern void _subplex(int n, double x[], void (*f)(int, double *, double *), double tol, int maxeval, double scale[], double *fx, int *neval, int nsmin, int nsmax, int icont, int nfstop, double fstop, int minf, double alpha, double _beta, double gamma, double delta, double psi, double omega, int irepl, int ifxsw, double bonus, double work[], int lwork, int iwork[], int _liwork, int *info);

#if defined(_VLARRAY)
extern void _optif9_r(int n, double x[], double typsiz[], double fscale, int method, int iexp, int ndigit, int maxiter, int iagflg, int iahflg, double dlt, double gradtl, double stepmx, double steptl, int iresult, double xpls[], double *fpls, double gpls[], int *iter, double work[], int lwork, int *info, double xx[], double *yy, double yyp[], int ldyyp2, double yyp2[][ldyyp2], int *irev);
#else
extern void _optif9_r(int n, double x[], double typsiz[], double fscale, int method, int iexp, int ndigit, int maxiter, int iagflg, int iahflg, double dlt, double gradtl, double stepmx, double steptl, int iresult, double xpls[], double *fpls, double gpls[], int *iter, double work[], int lwork, int *info, double xx[], double *yy, double yyp[], int ldyyp2, double yyp2[], int *irev);
#endif
extern void _optif0_r(int n, double x[], double xpls[], double *fpls, double work[], int lwork, int *info, double xx[], double yy, int *irev);
extern void _mng_r(int n, double x[], double d[], double v[], int lv, int iv[], int _liv, int *info, double *yy, double yyp[], int *irev);
extern void _mnf_r(int n, double x[], double d[], double v[], int lv, int iv[], int _liv, int *info, double *yy, int *irev);
extern void _mnh_r(int n, double x[], double d[], double v[], int lv, int iv[], int _liv, int *info, double *yy, double yyp[], double yypd[], int *irev);
extern void _subplex_r(int n, double x[], double tol, int maxeval, double scale[], double *fx, int *neval, int nsmin, int nsmax, int icont, int nfstop, double fstop, int minf, double alpha, double _beta, double gamma, double delta, double psi, double omega, int irepl, int ifxsw, double bonus, double work[], int lwork, int iwork[], int _liwork, int *info, double yy, int *irev);

/*
 * G2. Constrained optimization
 */
extern void _mngb(int n, double x[], double b[][2], void (*calcf)(int, double *, int *, double *), void (*calcg)(int, double *, int *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);
extern void _mnfb(int n, double x[], double b[][2], void (*calcf)(int, double *, int *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);
extern void _mnhb(int n, double x[], double b[][2], void (*calcf)(int, double *, int *, double *), void (*calcgh)(int, double *, int *, double *, double *), double d[], void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);

extern void _mngb_r(int n, double x[], double b[][2], double d[], double v[], int lv, int iv[], int _liv, int *info, double *yy, double yyp[], int *irev);
extern void _mnfb_r(int n, double x[], double b[][2], double d[], double v[], int lv, int iv[], int _liv, int *info, double *yy, int *irev);
extern void _mnhb_r(int n, double x[], double b[][2], double d[], double v[], int lv, int iv[], int _liv, int *info, double *yy, double yyp[], double yypd[], int *irev);

/*
 * H2a1a. 1-D finite _interval quadrature, integrand available via user-defined procedure
 */
extern void _qk15(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void _qk21(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void _qk31(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void _qk41(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void _qk51(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void _qk61(double (*f)(double), double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void _qng(double (*f)(double), double a, double b, double epsabs, double epsrel, double *result, double *abserr, int *neval, int *info);
extern void _qag(double (*f)(double), double a, double b, double epsabs, double epsrel, int key, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info);
extern void _qags(double (*f)(double), double a, double b, double epsabs, double epsrel, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info);
extern void _defint(double (*f)(double), double a, double b, double eps, int l, double *result, int *neval, int *info);

extern void _qk15_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void _qk21_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void _qk31_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void _qk41_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void _qk51_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void _qk61_r(double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void _qng_r(double a, double b, double epsabs, double epsrel, double *result, double *abserr, int *neval, int *info, double *xx, double yy, int *irev);
extern void _qag_r(double a, double b, double epsabs, double epsrel, int key, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info, double *xx, double yy, int *irev);
extern void _qags_r(double a, double b, double epsabs, double epsrel, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info, double *xx, double yy, int *irev);
extern void _defint_r(double a, double b, double eps, int l, double *result, int *neval, int *info, double *xx, double yy, int *irev);

/*
 * H2a1b. 1-D finite _interval quadrature, integrand available on a grid
 */
extern void _avint(int n, double x[], double y[], double a, double b, double *result, int *info);

/*
 * H2a2a. 1-D finite _interval quadrature, special integrand, integrand available via user-defined procedure
 */
extern void _qagp(double (*f)(double), double a, double b, int npts, double points[], double epsabs, double epsrel, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _qawc(double (*f)(double), double a, double b, double c, double epsabs, double epsrel, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info);
extern void _qaws(double (*f)(double), double a, double b, double alfa, double _beta, int integr, double epsabs, double epsrel, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info);
extern void _qawo(double (*f)(double), double a, double b, double omega, int integr, double epsabs, double epsrel, int _limit, int maxp1, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int liwork, int *info);

extern void _qagp_r(double a, double b, int npts, double points[], double epsabs, double epsrel, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int liwork, int *info, double *xx, double yy, int *irev);
extern void _qawc_r(double a, double b, double c, double epsabs, double epsrel, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info, double *xx, double yy, int *irev);
extern void _qaws_r(double a, double b, double alfa, double _beta, int integr, double epsabs, double epsrel, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info, double *xx, double yy, int *irev);
extern void _qawo_r(double a, double b, double omega, int integr, double epsabs, double epsrel, int _limit, int maxp1, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int liwork, int *info, double *xx, double yy, int *irev);

/*
 * H2a3a. 1-D semi-infinite _interval quadrature, integrand available via
user-defined procedure
 */
extern void _qk15i(double (*f)(double), double bound, int inf, double a, double b, double *result, double *abserr, double *resabs, double *resasc);
extern void _qagi(double (*f)(double), double bound, int inf, double epsabs, double epsrel, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info);
extern void _qawf(double (*f)(double), double a, double omega, int integr, double epsabs, int _limlst, int limit, int maxp1, double *result, double *abserr, int *neval, int *lst, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _dehint(double (*f)(double), double a, double eps, double *result, int *neval, int *l, int *info);
extern void _deoint(double (*f)(double), double a, double omega, int iw, double eps, double *result, int *neval, double *err, int *info);

extern void _qk15i_r(double bound, int inf, double a, double b, double *result, double *abserr, double *resabs, double *resasc, double *xx, double yy, int *irev);
extern void _qagi_r(double bound, int inf, double epsabs, double epsrel, int _limit, double *result, double *abserr, int *neval, int *last, double work[], int lwork, int iwork[], int *info, double *xx, double yy, int *irev);
extern void _qawf_r(double a, double omega, int integr, double epsabs, int _limlst, int limit, int maxp1, double *result, double *abserr, int *neval, int *lst, double work[], int lwork, int iwork[], int liwork, int *info, double *xx, double yy, int *irev);
extern void _dehint_r(double a, double eps, double *result, int *neval, int *l, int *info, double *xx, double yy, int *irev);
extern void _deoint_r(double a, double omega, int iw, double eps, double *result, int *neval, double *err, int *info, double *xx, double yy, int *irev);

/*
 * H2a4. 1-D infinite _interval quadrature
 */
extern void _deiint(double (*f)(double), double eps, double *result, int *neval, int *info);
extern void _deiint_r(double eps, double *result, int *neval, int *info, double *xx, double yy, int *irev);

/*
 * I1a1. Nonstiff initial value problem for ordinary differential equations (ODEs)
 */
extern void _derkf(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _derkf_int(int n, double t, double y[], double work[]);
extern void _dopri5(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info);
extern double _contd5(int ii, double x, double _rcont[], int icont[]);
extern void _dverk(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double tol, double c[], int lc, double work[], int lwork, int *info);
extern void _dverk_int(int n, double t, double y[], double c[], double work[]);
extern void _dop853(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info);
extern double _contd8(int ii, double x, double _rcont[], int icont[]);
extern void _deabm(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, int itstop, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _odex(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info);
extern double _contx1(int ii, double x, double _rcont[], int icont[]);
extern void _doprin(int n, void (*f)(int, double, double *, double *), double *t, double y[], double yp[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double *, double *, int, int *), int iout, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _odex2(int n, void (*f)(int, double, double *, double *), double *t, double y[], double yp[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double, double *, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info);
extern double _contx2(int ii, double x, double _rcont[], int icont[]);
extern void _retard(int n, void (*f)(int, double, double *, double *, double *, int *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*solout)(int, double, double, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info);
extern double _ylag(int ii, double x, double (*phi)(int, double), double _rcont[], int icont[]);

extern void _derkfa(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
void _dopri5a(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
void _dop853a(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _dverka(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double tend, double tol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _odexa(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _dopn43(int n, void (*f2)(int, double, double *, double *), double *t, double y[], double yp[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _dopn64(int n, void (*f2)(int, double, double *, double *), double *t, double y[], double yp[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _dopn86(int n, void (*f2)(int, double, double *, double *), double *t, double y[], double yp[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _dopn1210(int n, void (*f2)(int, double, double *, double *), double *t, double y[], double yp[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _odex2a(int n, void (*f2)(int, double, double *, double *), double *t, double y[], double yp[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _retarda(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _ylaga(int i, int n, double t,  double y[], void (*phi)(int, int, double, double *, int *), double work[], int iwork[], int *info);

extern void _derkf_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int _liwork, int *info, double *xx, double yy[], double yyp[], int *irev);
extern void _dopri5_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int *irtrn, int *irev);
extern void _dverk_r(int n, double *t, double y[], double tout, double tol, double c[], int lc, double work[], int lwork, int *info, double *tt, double yy[], double yyp[], int *irev);
extern void _dop853_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int *irtrn, int *irev);
extern void _deabm_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, int itstop, double work[], int lwork, int iwork[], int _liwork, int *info, double *tt, double yy[], double yyp[], int *irev);
extern void _odex_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int *irtrn, int *irev);
extern void _doprin_r(int n, double *t, double y[], double yp[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int _liwork, int *info, double *tt, double yy[], double yypp[], int *irtrn, int *irev);
extern void _odex2_r(int n, double *t, double y[], double yp[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int *irtrn, int *irev);
extern void _retard_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int *irtrn, int *irev);
extern void _ylag_r(int ii, double x, double _rcont[], int icont[], double *yy, int *irev);

extern void _derkfa_r(int n, double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yyp[], int *irev);
void _dopri5a_r(int n, double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yyp[], int *irev);
void _dop853a_r(int n, double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yyp[], int *irev);
extern void _dverka_r(int n, double *t, double y[], double tout, double tend, double tol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yyp[], int *irev);
extern void _odexa_r(int n, double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yyp[], int *irev);
extern void _dopn43_r(int n, double *t, double y[], double yp[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yypp[], int *irev);
extern void _dopn64_r(int n, double *t, double y[], double yp[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yypp[], int *irev);
extern void _dopn86_r(int n, double *t, double y[], double yp[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yypp[], int *irev);
extern void _dopn1210_r(int n, double *t, double y[], double yp[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yypp[], int *irev);
extern void _odex2a_r(int n, double *t, double y[], double yp[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yypp[], int *irev);
extern void _retarda_r(int n, double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yyp[], int *irev);

/*
 * I1a2. Stiff initial value problem for ordinary differential equations (ODEs)
 */
extern void _debdf(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, int itstop, void (*djac)(int, double, double *, int, double *), int idjac, double work[], int lwork, int iwork[], int _liwork, int *info);
extern void _radau5(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*jac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*mas)(int, int, double *), int imas, int mlmas, int mumas, void (*solout)(int, double, double, double *, int, double *, int *), int iout, double work[], int lwork, int iwork[], int _liwork, double cont[], int lcont, int *info);
extern double _contr5(int ii, double x, double cont[]);
extern void _radaup(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*jac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*mas)(int, int, double *), int imas, int mlmas, int mumas, void (*solout)(int, double, double, double *, int, double *, int *), int iout, double work[], int lwork, int iwork[], int _liwork, double cont[], int lcont, int *info);
extern double _contrp(int ii, double x, double cont[]);
extern void _radau(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*jac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*mas)(int, int, double *), int imas, int mlmas, int mumas, void (*solout)(int, double, double, double *, int, double *, int *), int iout, double work[], int lwork, int iwork[], int _liwork, double cont[], int lcont, int *info);
extern double _contra(int ii, double x, double cont[]);
extern void _rodas(int n, void (*f)(int, double, double *, double *), int ifcn, double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*jac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*dfx)(int, double, double *, double *), int idfx, void (*mas)(int, int, double *), int imas, int mlmas, int mumas, void (*solout)(int, double, double, double *, int, double *, int *), int iout, double work[], int lwork, int iwork[], int _liwork, double cont[], int lcont, int *info);
extern double _contro(int ii, double x, double cont[]);
extern void _seulex(int n, void (*f)(int, double, double *, double *), int ifcn, double *t, double y[], double tout, double *rtol, double *atol, int itol, void (*jac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*mas)(int, int, double *), int imas, int mlmas, int mumas, void (*solout)(int, double, double, double *, int, double *, int *, int *), int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info);
extern double _contex(int ii, double x, double _rcont[], int icont[]);
extern void _dassl(int n, void (*res)(int, double, double *, double *, double *, int *), double *t, double y[], double yp[], double tout, int iopt[], double *rtol, double *atol, void (*jac)(int, double, double *, double *, int, double *, double), double work[], int lwork, int iwork[], int _liwork, int *info);

extern void _debdf_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int mode, int itstop, int idjac, double work[], int lwork, int iwork[], int _liwork, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irev);
extern void _radau5_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int imas, int mlmas, int mumas, int iout, double work[], int lwork, int iwork[], int _liwork, double cont[], int lcont, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irtrn, int *irev);
extern void _radaup_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int imas, int mlmas, int mumas, int iout, double work[], int lwork, int iwork[], int _liwork, double cont[], int lcont, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irtrn, int *irev);
extern void _radau_r(int n, double *t, double y[], double tout, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int imas, int mlmas, int mumas, int iout, double work[], int lwork, int iwork[], int _liwork, double cont[], int lcont, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irtrn, int *irev);
extern void _rodas_r(int n, int ifcn, double *t, double y[], double tout, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int idfx, int imas, int mlmas, int mumas, int iout, double work[], int lwork, int iwork[], int _liwork, double cont[], int lcont, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irtrn, int *irev);
extern void _seulex_r(int n, int ifcn, double *t, double y[], double tout, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int imas, int mlmas, int mumas, int iout, double work[], int lwork, int iwork[], int _liwork, double _rcont[], int lrcont, int icont[], int licont, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irtrn, int *irev);
extern void _dassl_r(int n, double *t, double y[], double yp[], double tout, int iopt[], double *rtol, double *atol, double work[], int lwork, int iwork[], int _liwork, int *info, double *tt, double yyp[], int ldyypd, double yypd[], double *cj, int *ires, int *irev);

extern void _radaua(int n, void (*f)(int, double, double *, double *), double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, void (*fjac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*fmas)(int, int, double *), int imas, int mlmas, int mumas, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _rodasa(int n, void (*f)(int, double, double *, double *), int ifcn, double *t, double y[], double tout, double tend, double rtol[], double atol[], int itol, void (*fjac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*fdfx)(int, double, double *, double *), int idfx, void (*fmas)(int, int, double *), int imas, int mlmas, int mumas, int mode, double work[], int lwork, int iwork[], int liwork, int *info);
extern void _seulexa(int n, void (*f)(int, double, double *, double *), int ifcn, double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, void (*fjac)(int, double, double *, int, double *), int ijac, int mljac, int mujac, void (*fmas)(int, int, double *), int imas, int mlmas, int mumas, int mode, double work[], int lwork, int iwork[], int liwork, int *info);

extern void _radaua_r(int n, double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int imas, int mlmas, int mumas, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irev);
extern void _rodasa_r(int n, int ifcn, double *t, double y[], double tout, double tend, double rtol[], double atol[], int itol, int ijac, int mljac, int mujac, int idfx, int imas, int mlmas, int mumas, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irev);
extern void _seulexa_r(int n, int ifcn, double *t, double y[], double tout, double tend, double *rtol, double *atol, int itol, int ijac, int mljac, int mujac, int imas, int mlmas, int mumas, int mode, double work[], int lwork, int iwork[], int liwork, int *info, double *tt, double yy[], double yyp[], int ldyypd, double yypd[], int *irev);

/*
 * J1a1. One-dimensional real fast Fourier transforms
 */
extern void _rfft1f(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _rfft1b(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _rfft1i(int n, double wsave[], int lwsave, int *info);
extern void _rfftmf(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _rfftmb(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _rfftmi(int n, double wsave[], int lwsave, int *info);

/*
 * J1a2. One-dimensional complex fast Fourier transforms
 */
extern void _cfft1f(int n, int inc, doublecomplex c[], int lc, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cfft1b(int n, int inc, doublecomplex c[], int lc, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cfft1i(int n, double wsave[], int lwsave, int *info);
extern void _cfftmf(int lot, int jump, int n, int inc, doublecomplex c[], int lc, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cfftmb(int lot, int jump, int n, int inc, doublecomplex c[], int lc, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cfftmi(int n, double wsave[], int lwsave, int *info);

/*
 * J1a3. One-dimensional trigonometric fast Fourier transforms
 */
extern void _sint1f(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _sint1b(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _sint1i(int n, double wsave[], int lwsave, int *info);
extern void _sintmf(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _sintmb(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _sintmi(int n, double wsave[], int lwsave, int *info);
extern void _cost1f(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cost1b(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cost1i(int n, double wsave[], int lwsave, int *info);
extern void _costmf(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _costmb(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _costmi(int n, double wsave[], int lwsave, int *info);

extern void _sinq1f(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _sinq1b(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _sinq1i(int n, double wsave[], int lwsave, int *info);
extern void _sinqmf(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _sinqmb(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _sinqmi(int n, double wsave[], int lwsave, int *info);
extern void _cosq1f(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cosq1b(int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cosq1i(int n, double wsave[], int lwsave, int *info);
extern void _cosqmf(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cosqmb(int lot, int jump, int n, int inc, double r[], int lr, double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cosqmi(int n, double wsave[], int lwsave, int *info);

/*
 * J1b. Multidimensional fast Fourier transforms
 */
#if defined(_VLARRAY)
extern void _rfft2f(int l, int m, int ldr, double r[][ldr], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _rfft2b(int l, int m, int ldr, double r[][ldr], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _rfft2c(int l, int m, int ldr, double r[][ldr], int ldc, doublecomplex c[][ldc], int *info);
#else
extern void _rfft2f(int l, int m, int ldr, double r[], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _rfft2b(int l, int m, int ldr, double r[], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _rfft2c(int l, int m, int ldr, double r[], int ldc, doublecomplex c[], int *info);
#endif
extern void _rfft2i(int l, int m, double wsave[], int lwsave, int *info);
#if defined(_VLARRAY)
extern void _cfft2f(int l, int m, int ldc, doublecomplex c[][ldc], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cfft2b(int l, int m, int ldc, doublecomplex c[][ldc], double wsave[], int lwsave, double work[], int lwork, int *info);
#else
extern void _cfft2f(int l, int m, int ldc, doublecomplex c[], double wsave[], int lwsave, double work[], int lwork, int *info);
extern void _cfft2b(int l, int m, int ldc, doublecomplex c[], double wsave[], int lwsave, double work[], int lwork, int *info);
#endif
extern void _cfft2i(int l, int m, double wsave[], int lwsave, int *info);

/*
 * K1b1. Unconstrained nonlinear least squares approximation
 */
#if defined(_VLARRAY)
extern void _lmder(void (*fcn)(int, int, double *, double *, int, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info);
extern void _lmder1(void (*fcn)(int, int, double *, double *, int, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double tol, int ipvt[], double work[], int lwork, int *info);
extern void _lmstr(void (*fcn)(int, int, double *, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info);
extern void _lmstr1(void (*fcn)(int, int, double *, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double tol, int ipvt[], double work[], int lwork, int *info);
extern void _lmdif(void (*fcn)(int, int, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, int ipvt[], double work[], int lwork, int *info);
extern void _lmdif1(void (*fcn)(int, int, double *, double *, int *), int m, int n, double x[], double fvec[], double tol, double work[], int lwork, int iwork[], int *info);
extern void _covar(int n, int ldr, double r[][ldr], int ipvt[], double tol, double work[], int *info);
#else
extern void _lmder(void (*fcn)(int, int, double *, double *, int, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info);
extern void _lmder1(void (*fcn)(int, int, double *, double *, int, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double tol, int ipvt[], double work[], int lwork, int *info);
extern void _lmstr(void (*fcn)(int, int, double *, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info);
extern void _lmstr1(void (*fcn)(int, int, double *, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double tol, int ipvt[], double work[], int lwork, int *info);
extern void _lmdif(void (*fcn)(int, int, double *, double *, int *), int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, int ipvt[], double work[], int lwork, int *info);
extern void _lmdif1(void (*fcn)(int, int, double *, double *, int *), int m, int n, double x[], double fvec[], double tol, double work[], int lwork, int iwork[], int *info);
extern void _covar(int n, int ldr, double r[], int ipvt[], double tol, double work[], int *info);
#endif
extern void _n2g(int n, int p, double x[], void (*calcr)(int, int, double *, int *, double *), void (*calcj)(int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);
extern void _n2f(int n, int p, double x[], void (*calcr)(int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);
extern void _n2p(int n, int md, int p, double x[], void (*calcr)(int, int, int, int *, int, double *, int *, double *), void (*calcj)(int, int, int, int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);

#if defined(_VLARRAY)
extern void _lmder_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void _lmder1_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double tol, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void _lmstr_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], double yypr[], int *irev);
extern void _lmstr1_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double tol, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], double yypr[], int *irev);
extern void _lmdif_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[][ldfjac], double ftol, double xtol, double gtol, int maxfev, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
#else
extern void _lmder_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void _lmder1_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double tol, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
extern void _lmstr_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double diag[], int mode, double factor, int nprint, int *nfev, int *njev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], double yypr[], int *irev);
extern void _lmstr1_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double tol, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], double yypr[], int *irev);
extern void _lmdif_r(int m, int n, double x[], double fvec[], int ldfjac, double fjac[], double ftol, double xtol, double gtol, int maxfev, double epsfcn, double diag[], int mode, double factor, int nprint, int *nfev, int ipvt[], double work[], int lwork, int *info, double xx[], double yy[], int *irev);
#endif
extern void _lmdif1_r(int m, int n, double x[], double fvec[], double tol, double work[], int lwork, int iwork[], int *info, double xx[], double yy[], int *irev);
#if defined(_VLARRAY)
extern void _n2g_r(int n, int p, double x[], double v[], int lv, int iv[], int _liv, int *info, double yy[], int ldyyp, double yyp[][ldyyp], int *irev);
extern void _n2p_r(int n, int md, int p, double x[], double v[], int lv, int iv[], int _liv, int *info, int *m1, int *m2, double yy[], int ldyyp, double yyp[][ldyyp], int *irev);
#else
extern void _n2g_r(int n, int p, double x[], double v[], int lv, int iv[], int _liv, int *info, double yy[], int ldyyp, double yyp[], int *irev);
extern void _n2p_r(int n, int md, int p, double x[], double v[], int lv, int iv[], int _liv, int *info, int *m1, int *m2, double yy[], int ldyyp, double yyp[], int *irev);
#endif
extern void _n2f_r(int n, int p, double x[], double v[], int lv, int iv[], int _liv, int *info, double yy[], int *irev);

/*
 * K1b2. Constrained nonlinear least squares approximation
 */
extern void _n2gb(int n, int p, double x[], double b[][2], void (*calcr)(int, int, double *, int *, double *), void (*calcj)(int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);
extern void _n2fb(int n, int p, double x[], double b[][2], void (*calcr)(int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);
extern void _n2pb(int n, int nd, int p, double x[], double b[][2], void (*calcr)(int, int, int, int *, int, double *, int *, double *), void (*calcj)(int, int, int, int, int, double *, int *, double *), void (*itsum)(int, double *, int, int, int, double), double v[], int lv, int iv[], int _liv, int *info);

#if defined(_VLARRAY)
extern void _n2gb_r(int n, int p, double x[], double b[][2], double v[], int lv, int iv[], int _liv, int *info, double yy[], int ldyyp, double yyp[][ldyyp], int *irev);
extern void _n2pb_r(int n, int nd, int p, double x[], double b[][2], double v[], int lv, int iv[], int _liv, int *info, int *n1, int *n2, double yy[], int ldyyp, double yyp[][ldyyp], int *irev);
#else
extern void _n2gb_r(int n, int p, double x[], double b[][2], double v[], int lv, int iv[], int _liv, int *info, double yy[], int ldyyp, double yyp[], int *irev);
extern void _n2pb_r(int n, int nd, int p, double x[], double b[][2], double v[], int lv, int iv[], int _liv, int *info, int *n1, int *n2, double yy[], int ldyyp, double yyp[], int *irev);
#endif
extern void _n2fb_r(int n, int p, double x[], double b[][2], double v[], int lv, int iv[], int _liv, int *info, double yy[], int *irev);

/*
 * L6a5. Exponential random numbers
 */
extern double _zigexp(long (*iurand)(void), double (*durand)(void), double theta);
extern void _init_zigexp(int n);
extern void _zigexp_r(double theta, double *r, int ii, double yy, int *irev);
extern void _init_zigexp_r(int n);

/*
 * L6a7. Gamma random numbers
 */
extern double _rgama(double (*durand)(void), double (*dnrand)(void), double alpha, double _beta);
extern void _rgama_r(double alpha, double _beta, double *r, double yy, int *irev);

/*
 * L6a14. Normal random numbers
 */
extern double _zignorm(long (*iurand)(void), double (*durand)(void), double mu, double _sigma);
extern void _init_zignorm(int n);
extern void _zignorm_r(double mu, double _sigma, double *r, int ii, double yy, int *irev);
extern void _init_zignorm_r(int n);

/*
 * L6a21. Uniform random numbers
 */
extern void _init_genrand(unsigned long s);
extern void _init_by_array(unsigned long init_key[], int key_length);
extern unsigned long _genrand_int32(void);
extern long _genrand_int31(void);
extern double _genrand_real1(void);
extern double _genrand_real2(void);
extern double _genrand_real3(void);
extern double _genrand_res53(void);

extern void _init_genrand64(unsigned long long s);
extern void _init_by_array64(unsigned long long initkey[], unsigned long long keylength);
extern unsigned long long _genrand64_int64(void);
extern long long _genrand64_int63(void);
extern double _genrand64_real1(void);
extern double _genrand64_real2(void);
extern double _genrand64_real3(void);

extern void _ran_start(long seed);
extern void _ran_array(long aa[], int n);
extern long _ran_arr_next(void);
extern void _ranf_start(long seed);
extern void _ranf_array(double aa[], int n);
extern double _ranf_arr_next(void);

extern void _srand48(long seed);
extern unsigned short *seed48(unsigned short xseed[]);
extern void _lcong48(unsigned short p[]);
extern double _drand48(void);
extern double _erand48(unsigned short xseed[]);
extern long _lrand48(void);
extern long _nrand48(unsigned short xseed[]);
extern long _mrand48(void);
extern long _jrand48(unsigned short xseed[]);

/*
 * R. Service routines
 */
extern double _dlamch(char cmach);
extern float _slamch(char cmach);
extern double _d1mach(int i);
extern float _r1mach(int i);
extern int _i1mach(int i);

/*
 * Z. Others
 */
#if defined(_VLARRAY)
extern void _dlatme(int n, char dist, int iseed[], double d[], int mode, double cond, double dmax, char _ei[], char rsign, char upper, char _sim, double ds[], int modes, double conds, int kl, int ku, double anorm, int lda, double a[][lda], double work[], int *info);
extern void _dlatmr(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, char rsign, char grade, double dl[], int model, double condl, double dr[], int moder, double condr, char pivtng, int ipivot[], int kl, int ku, double sparse, double anorm, char pack, int lda, double a[][lda], int iwork[], int *info);
extern void _dlatms(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int kl, int ku, char pack, int lda, double a[][lda], double work[], int *info);
extern void _dlatmt(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int rank, int kl, int ku, char pack, int lda, double a[][lda], double work[], int *info);
extern void _zlatme(int n, char dist, int iseed[], doublecomplex d[], int mode, double cond, doublecomplex dmax, char rsign, char upper, char _sim, double ds[], int modes, double conds, int kl, int ku, double anorm, int lda, doublecomplex a[][lda], doublecomplex work[], int *info);
extern void _zlatmr(int m, int n, char dist, int iseed[], char sym, doublecomplex d[], int mode, double cond, doublecomplex dmax, char rsign, char grade, doublecomplex dl[], int model, double condl, doublecomplex dr[], int moder, double condr, char pivtng, int ipivot[], int kl, int ku, double sparse, double anorm, char pack, int lda, doublecomplex a[][lda], int iwork[], int *info);
extern void _zlatms(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int kl, int ku, char pack, int lda, doublecomplex a[][lda], doublecomplex work[], int *info);
extern void _zlatmt(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int rank, int kl, int ku, char pack, int lda, doublecomplex a[][lda], doublecomplex work[], int *info);
#else
extern void _dlatme(int n, char dist, int iseed[], double d[], int mode, double cond, double dmax, char _ei[], char rsign, char upper, char _sim, double ds[], int modes, double conds, int kl, int ku, double anorm, int lda, double a[], double work[], int *info);
extern void _dlatmr(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, char rsign, char grade, double dl[], int model, double condl, double dr[], int moder, double condr, char pivtng, int ipivot[], int kl, int ku, double sparse, double anorm, char pack, int lda, double a[], int iwork[], int *info);
extern void _dlatms(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int kl, int ku, char pack, int lda, double a[], double work[], int *info);
extern void _dlatmt(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int rank, int kl, int ku, char pack, int lda, double a[], double work[], int *info);
extern void _zlatme(int n, char dist, int iseed[], doublecomplex d[], int mode, double cond, doublecomplex dmax, char rsign, char upper, char _sim, double ds[], int modes, double conds, int kl, int ku, double anorm, int lda, doublecomplex a[], doublecomplex work[], int *info);
extern void _zlatmr(int m, int n, char dist, int iseed[], char sym, doublecomplex d[], int mode, double cond, doublecomplex dmax, char rsign, char grade, doublecomplex dl[], int model, double condl, doublecomplex dr[], int moder, double condr, char pivtng, int ipivot[], int kl, int ku, double sparse, double anorm, char pack, int lda, doublecomplex a[], int iwork[], int *info);
extern void _zlatms(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int kl, int ku, char pack, int lda, doublecomplex a[], doublecomplex work[], int *info);
extern void _zlatmt(int m, int n, char dist, int iseed[], char sym, double d[], int mode, double cond, double dmax, int rank, int kl, int ku, char pack, int lda, doublecomplex a[], doublecomplex work[], int *info);
#endif

#if defined(__cplusplus)
}
#endif
