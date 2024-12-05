/****************************************
 *                                      *
 *  XLPack Numerical Library            *
 *  Version 7.0 (January 31, 2023)      *
 *  (C) 2014-2023  K Technologies       *
 *                                      *
 ****************************************/
#pragma once

/* M0 */

#define d1num	_d1num

#define cadd	_cadd
#define cadd_sub	_cadd_sub
#define cdadd	_cdadd
#define cdadd_sub	_cdadd_sub
#define cddiv	_cddiv
#define cddiv_sub	_cddiv_sub
#define cdiv	_cdiv
#define cdiv_sub	_cdiv_sub
#define cdmul	_cdmul
#define cdmul_sub	_cdmul_sub
#define cdpow	_cdpow
#define cdpow_sub	_cdpow_sub
#define cdsub	_cdsub
#define cdsub_sub	_cdsub_sub
#define cipow	_cipow
#define cipow_sub	_cipow_sub
#define cminus	_cminus
#define cminus_sub	_cminus_sub
#define cmplx	_cmplx
#define cmplx_sub	_cmplx_sub
#define cmul	_cmul
#define cmul_sub	_cmul_sub
#define conj_sub	_conj_sub
#define cpolar	_cpolar
#define cpolar_sub	_cpolar_sub
#define cpow_sub	_cpow_sub
#define cproj_sub	_cproj_sub
#define csub	_csub
#define csub_sub	_csub_sub
#define dcdiv	_dcdiv
#define dcdiv_sub	_dcdiv_sub
#define dcsub	_dcsub
#define dcsub_sub	_dcsub_sub

#define factorial	_factorial

#define ccbrt	_ccbrt
#define ccbrt_sub	_ccbrt_sub
#define csqrt_sub	_csqrt_sub

#define cacosh_sub	_cacosh_sub
#define cacos_sub	_cacos_sub
#define casinh_sub	_casinh_sub
#define casin_sub	_casin_sub
#define catanh_sub	_catanh_sub
#define catan_sub	_catan_sub
#define ccosh_sub	_ccosh_sub
#define ccos_sub	_ccos_sub
#define ccot	_ccot
#define ccot_sub	_ccot_sub
#define cexpm1	_cexpm1
#define cexpm1_sub	_cexpm1_sub
#define cexp_sub	_cexp_sub
#define clog1p	_clog1p
#define clog1p_sub	_clog1p_sub
#define clog_sub	_clog_sub
#define cos_pi	_cos_pi
#define csinh_sub	_csinh_sub
#define csin_sub	_csin_sub
#define ctanh_sub	_ctanh_sub
#define ctan_sub	_ctan_sub
#define powm1	_powm1
#define sin_pi	_sin_pi
#define sqrt1pm1	_sqrt1pm1

#define ccot	_ccot
#define cexpm1	_cexpm1
#define clog1p	_clog1p

#define dconst	_dconst

#define slamch	_slamch
#define dlamch	_dlamch
#define r1mach	_r1mach
#define d1mach	_d1mach
#define i1mach	_i1mach

/* M1 */

#define dasum	_dasum
#define daxpy	_daxpy
#define dcopy	_dcopy
#define ddot	_ddot
#define dnrm2	_dnrm2
#define drot	_drot
#define drotg	_drotg
#define drotm	_drotm
#define drotmg	_drotmg
#define dscal	_dscal
#define dswap	_dswap
#define idamax	_idamax
#define lsame	_lsame

#define dgbmv	_dgbmv
#define dgemv	_dgemv
#define dger	_dger
#define dsbmv	_dsbmv
#define dspmv	_dspmv
#define dspr	_dspr
#define dspr2	_dspr2
#define dsymv	_dsymv
#define dsyr	_dsyr
#define dsyr2	_dsyr2
#define dtbmv	_dtbmv
#define dtbsv	_dtbsv
#define dtpmv	_dtpmv
#define dtpsv	_dtpsv
#define dtrmv	_dtrmv
#define dtrsv	_dtrsv

#define dgemm	_dgemm
#define dsymm	_dsymm
#define dsyr2k	_dsyr2k
#define dsyrk	_dsyrk
#define dtrmm	_dtrmm
#define dtrsm	_dtrsm

#define dlangb	_dlangb
#define dlange	_dlange
#define dlangt	_dlangt
#define dlansb	_dlansb
#define dlansp	_dlansp
#define dlanst	_dlanst
#define dlansy	_dlansy
#define dlantr	_dlantr

#define dgbcon	_dgbcon
#define dgbsv	_dgbsv
#define dgbsvx	_dgbsvx
#define dgbtrf	_dgbtrf
#define dgbtrs	_dgbtrs
#define dgecon	_dgecon
#define dgesv	_dgesv
#define dgesvx	_dgesvx
#define dgetrf	_dgetrf
#define dgetri	_dgetri
#define dgetrs	_dgetrs
#define dgtcon	_dgtcon
#define dgtsv	_dgtsv
#define dgtsvx	_dgtsvx
#define dgttrf	_dgttrf
#define dgttrs	_dgttrs
#define dsgesv	_dsgesv

#define dtbcon	_dtbcon
#define dtbtrs	_dtbtrs
#define dtpcon	_dtpcon
#define dtptri	_dtptri
#define dtptrs	_dtptrs
#define dtrcon	_dtrcon
#define dtrtri	_dtrtri
#define dtrtrs	_dtrtrs

#define dspcon	_dspcon
#define dspsv	_dspsv
#define dspsvx	_dspsvx
#define dsptrf	_dsptrf
#define dsptri	_dsptri
#define dsptrs	_dsptrs
#define dsycon	_dsycon
#define dsysv	_dsysv
#define dsysvx	_dsysvx
#define dsytrf	_dsytrf
#define dsytri	_dsytri
#define dsytrs	_dsytrs

#define dpocon	_dpocon
#define dposv	_dposv
#define dposvx	_dposvx
#define dpotrf	_dpotrf
#define dpotri	_dpotri
#define dpotrs	_dpotrs
#define dppcon	_dppcon
#define dppsv	_dppsv
#define dppsvx	_dppsvx
#define dpptrf	_dpptrf
#define dpptri	_dpptri
#define dpptrs	_dpptrs
#define dsposv	_dsposv

#define dpbcon	_dpbcon
#define dpbsv	_dpbsv
#define dpbsvx	_dpbsvx
#define dpbtrf	_dpbtrf
#define dpbtrs	_dpbtrs
#define dptcon	_dptcon
#define dptsv	_dptsv
#define dptsvx	_dptsvx
#define dpttrf	_dpttrf
#define dpttrs	_dpttrs

#define ddisna	_ddisna
#define dopgtr	_dopgtr
#define dopmtr	_dopmtr
#define dorgtr	_dorgtr
#define dormtr	_dormtr
#define dpteqr	_dpteqr
#define dsbev	_dsbev
#define dsbevd	_dsbevd
#define dsbevx	_dsbevx
#define dsbtrd	_dsbtrd
#define dspev	_dspev
#define dspevd	_dspevd
#define dspevx	_dspevx
#define dsptrd	_dsptrd
#define dstebz	_dstebz
#define dstedc	_dstedc
#define dstein	_dstein
#define dstemr	_dstemr
#define dsteqr	_dsteqr
#define dsterf	_dsterf
#define dstev	_dstev
#define dstevd	_dstevd
#define dstevr	_dstevr
#define dstevx	_dstevx
#define dsyev	_dsyev
#define dsyevd	_dsyevd
#define dsyevr	_dsyevr
#define dsyevx	_dsyevx
#define dsytrd	_dsytrd

#define dgebak	_dgebak
#define dgebal	_dgebal
#define dgees	_dgees
#define dgeesx	_dgeesx
#define dgeesx_r	_dgeesx_r
#define dgees_r	_dgees_r
#define dgeev	_dgeev
#define dgeevx	_dgeevx
#define dgehrd	_dgehrd
#define dhsein	_dhsein
#define dhseqr	_dhseqr
#define dorghr	_dorghr
#define dormhr	_dormhr
#define dtrevc3	_dtrevc3
#define dtrexc	_dtrexc
#define dtrsen	_dtrsen
#define dtrsna	_dtrsna
#define dtrsyl	_dtrsyl

#define dsbgv	_dsbgv
#define dsbgvd	_dsbgvd
#define dsbgvx	_dsbgvx
#define dspgv	_dspgv
#define dspgvd	_dspgvd
#define dspgvx	_dspgvx
#define dsygv	_dsygv
#define dsygvd	_dsygvd
#define dsygvx	_dsygvx

#define dgges	_dgges
#define dggesx	_dggesx
#define dggesx_r	_dggesx_r
#define dgges_r	_dgges_r
#define dggev	_dggev
#define dggevx	_dggevx

#define dgelqf	_dgelqf
#define dgeqp3	_dgeqp3
#define dgeqrf	_dgeqrf
#define dorglq	_dorglq
#define dorgqr	_dorgqr
#define dormlq	_dormlq
#define dormqr	_dormqr

#define dgejsv	_dgejsv
#define dgesdd	_dgesdd
#define dgesvd	_dgesvd
#define dgesvdq	_dgesvdq
#define dgesvdx	_dgesvdx
#define dggsvd3	_dggsvd3

#define dgels	_dgels
#define dgelsd	_dgelsd
#define dgelss	_dgelss
#define dgelsy	_dgelsy
#define dgetsls	_dgetsls

#define dggglm	_dggglm
#define dgglse	_dgglse

#define dlatme	_dlatme
#define dlatmr	_dlatmr
#define dlatms	_dlatms
#define dlatmt	_dlatmt

#define dgecov	_dgecov
#define dgecovs	_dgecovs
#define dgecovy	_dgecovy

/* M2 */

#define dcabs1	_dcabs1
#define dzasum	_dzasum
#define dznrm2	_dznrm2
#define izamax	_izamax
#define zaxpy	_zaxpy
#define zcopy	_zcopy
#define zdotc	_zdotc
#define zdotc_sub	_zdotc_sub
#define zdotu	_zdotu
#define zdotu_sub	_zdotu_sub
#define zdrot	_zdrot
#define zdscal	_zdscal
#define zrot	_zrot
#define zrotg	_zrotg
#define zscal	_zscal
#define zswap	_zswap

#define zgbmv	_zgbmv
#define zgemv	_zgemv
#define zgerc	_zgerc
#define zgeru	_zgeru
#define zhbmv	_zhbmv
#define zhemv	_zhemv
#define zher	_zher
#define zher2	_zher2
#define zhpmv	_zhpmv
#define zhpr	_zhpr
#define zhpr2	_zhpr2
#define zsbmv	_zsbmv
#define zspmv	_zspmv
#define zspr	_zspr
#define zspr2	_zspr2
#define zsymv	_zsymv
#define zsyr	_zsyr
#define zsyr2	_zsyr2
#define ztbmv	_ztbmv
#define ztbsv	_ztbsv
#define ztpmv	_ztpmv
#define ztpsv	_ztpsv
#define ztrmv	_ztrmv
#define ztrsv	_ztrsv

#define zgemm	_zgemm
#define zhemm	_zhemm
#define zher2k	_zher2k
#define zherk	_zherk
#define zsymm	_zsymm
#define zsyr2k	_zsyr2k
#define zsyrk	_zsyrk
#define ztrmm	_ztrmm
#define ztrsm	_ztrsm

#define zlangb	_zlangb
#define zlange	_zlange
#define zlangt	_zlangt
#define zlanhb	_zlanhb
#define zlanhe	_zlanhe
#define zlanhp	_zlanhp
#define zlanht	_zlanht
#define zlansb	_zlansb
#define zlansp	_zlansp
#define zlansy	_zlansy
#define zlantr	_zlantr

#define zcgesv	_zcgesv
#define zgbcon	_zgbcon
#define zgbsv	_zgbsv
#define zgbsvx	_zgbsvx
#define zgbtrf	_zgbtrf
#define zgbtrs	_zgbtrs
#define zgecon	_zgecon
#define zgesv	_zgesv
#define zgesvx	_zgesvx
#define zgetrf	_zgetrf
#define zgetri	_zgetri
#define zgetrs	_zgetrs
#define zgtcon	_zgtcon
#define zgtsv	_zgtsv
#define zgtsvx	_zgtsvx
#define zgttrf	_zgttrf
#define zgttrs	_zgttrs
#define zspcon	_zspcon
#define zspsv	_zspsv
#define zspsvx	_zspsvx
#define zsptrf	_zsptrf
#define zsptri	_zsptri
#define zsptrs	_zsptrs
#define zsycon	_zsycon
#define zsysv	_zsysv
#define zsysvx	_zsysvx
#define zsytrf	_zsytrf
#define zsytri	_zsytri
#define zsytrs	_zsytrs

#define ztbcon	_ztbcon
#define ztbtrs	_ztbtrs
#define ztpcon	_ztpcon
#define ztptri	_ztptri
#define ztptrs	_ztptrs
#define ztrcon	_ztrcon
#define ztrtri	_ztrtri
#define ztrtrs	_ztrtrs

#define zhecon	_zhecon
#define zhesv	_zhesv
#define zhesvx	_zhesvx
#define zhetrf	_zhetrf
#define zhetri	_zhetri
#define zhetrs	_zhetrs
#define zhpcon	_zhpcon
#define zhpsv	_zhpsv
#define zhpsvx	_zhpsvx
#define zhptrf	_zhptrf
#define zhptri	_zhptri
#define zhptrs	_zhptrs

#define zcposv	_zcposv
#define zpocon	_zpocon
#define zposv	_zposv
#define zposvx	_zposvx
#define zpotrf	_zpotrf
#define zpotri	_zpotri
#define zpotrs	_zpotrs
#define zppcon	_zppcon
#define zppsv	_zppsv
#define zppsvx	_zppsvx
#define zpptrf	_zpptrf
#define zpptri	_zpptri
#define zpptrs	_zpptrs

#define zpbcon	_zpbcon
#define zpbsv	_zpbsv
#define zpbsvx	_zpbsvx
#define zpbtrf	_zpbtrf
#define zpbtrs	_zpbtrs
#define zptcon	_zptcon
#define zptsv	_zptsv
#define zptsvx	_zptsvx
#define zpttrf	_zpttrf
#define zpttrs	_zpttrs

#define zhbev	_zhbev
#define zhbevd	_zhbevd
#define zhbevx	_zhbevx
#define zhbtrd	_zhbtrd
#define zheev	_zheev
#define zheevd	_zheevd
#define zheevr	_zheevr
#define zheevx	_zheevx
#define zhetrd	_zhetrd
#define zhpev	_zhpev
#define zhpevd	_zhpevd
#define zhpevx	_zhpevx
#define zhptrd	_zhptrd
#define zpteqr	_zpteqr
#define zstedc	_zstedc
#define zstein	_zstein
#define zstemr	_zstemr
#define zsteqr	_zsteqr
#define zungtr	_zungtr
#define zunmtr	_zunmtr
#define zupgtr	_zupgtr
#define zupmtr	_zupmtr

#define zgebak	_zgebak
#define zgebal	_zgebal
#define zgees	_zgees
#define zgeesx	_zgeesx
#define zgeesx_r	_zgeesx_r
#define zgees_r	_zgees_r
#define zgeev	_zgeev
#define zgeevx	_zgeevx
#define zgehrd	_zgehrd
#define zhsein	_zhsein
#define zhseqr	_zhseqr
#define ztrevc3	_ztrevc3
#define ztrexc	_ztrexc
#define ztrsen	_ztrsen
#define ztrsna	_ztrsna
#define ztrsyl	_ztrsyl
#define zunghr	_zunghr
#define zunmhr	_zunmhr

#define zhbgv	_zhbgv
#define zhbgvd	_zhbgvd
#define zhbgvx	_zhbgvx
#define zhegv	_zhegv
#define zhegvd	_zhegvd
#define zhegvx	_zhegvx
#define zhpgv	_zhpgv
#define zhpgvd	_zhpgvd
#define zhpgvx	_zhpgvx

#define zgges	_zgges
#define zggesx	_zggesx
#define zggesx_r	_zggesx_r
#define zgges_r	_zgges_r
#define zggev	_zggev
#define zggevx	_zggevx

#define zgelqf	_zgelqf
#define zgeqp3	_zgeqp3
#define zgeqrf	_zgeqrf
#define zunglq	_zunglq
#define zungqr	_zungqr
#define zunmlq	_zunmlq
#define zunmqr	_zunmqr

#define zgejsv	_zgejsv
#define zgesdd	_zgesdd
#define zgesvd	_zgesvd
#define zgesvdq	_zgesvdq
#define zgesvdx	_zgesvdx
#define zggsvd3	_zggsvd3

#define zgels	_zgels
#define zgelsd	_zgelsd
#define zgelss	_zgelss
#define zgelsy	_zgelsy
#define zgetsls	_zgetsls

#define zggglm	_zggglm
#define zgglse	_zgglse

#define zlatme	_zlatme
#define zlatmr	_zlatmr
#define zlatms	_zlatms
#define zlatmt	_zlatmt

#define zgecov	_zgecov
#define zgecovs	_zgecovs
#define zgecovy	_zgecovy

/* M3 */

#define laguerre	_laguerre
#define alaguerre	_alaguerre
#define legendre	_legendre
#define legendred	_legendred
#define alegendre	_alegendre
#define sharmonic	_sharmonic
#define sharmonic_sub	_sharmonic_sub
#define sharmonicr	_sharmonicr
#define sharmonici	_sharmonici
#define hermite	_hermite
#define chebt	_chebt
#define chebtd	_chebtd
#define chebu	_chebu
#define chebs	_chebs

#define gegenbauer	_gegenbauer
#define gegenbauerd1	_gegenbauerd1
#define gegenbauerd	_gegenbauerd
#define jacobi	_jacobi
#define jacobid1	_jacobid1
#define jacobid2	_jacobid2
#define jacobid	_jacobid

#define li	_li
#define ei	_ei
#define e1	_e1
#define en	_en
#define spence	_spence

#define ci	_ci
#define si	_si
#define chi	_chi
#define shi	_shi

#define tgamma1pm1	_tgamma1pm1
#define lgammas	_lgammas
#define rgamma	_rgamma
#define tgammaratio	_tgammaratio
#define tgammadratio	_tgammadratio
#define cgamma	_cgamma
#define clgamma	_clgamma
#define crgamma	_crgamma
#define cgamma_sub	_cgamma_sub
#define clgamma_sub	_clgamma_sub
#define crgamma_sub	_crgamma_sub
#define poch	_poch
#define poch1	_poch1

#define beta	_beta
#define lbeta	_lbeta
#define cbeta	_cbeta
#define cbeta_sub	_cbeta_sub
#define clbeta	_clbeta
#define clbeta_sub	_clbeta_sub

#define digamma	_digamma
#define trigamma	_trigamma
#define polygamma	_polygamma
#define cdigamma	_cdigamma
#define cdigamma_sub	_cdigamma_sub

#define gammai	_gammai
#define gammaic	_gammaic
#define gammait	_gammait
#define gammap	_gammap
#define gammaq	_gammaq
#define gammapi	_gammapi
#define gammaqi	_gammaqi
#define gammapia	_gammapia
#define gammaqia	_gammaqia
#define gammapd	_gammapd

#define betax	_betax
#define betaxc	_betaxc
#define ibeta	_ibeta
#define ibetac	_ibetac
#define ibetai	_ibetai
#define ibetaci	_ibetaci
#define ibetaia	_ibetaia
#define ibetacia	_ibetacia
#define ibetaib	_ibetaib
#define ibetacib	_ibetacib
#define ibetad	_ibetad

#define zeta	_zeta

#define erfi	_erfi
#define erfci	_erfci
#define dawson	_dawson
#define fresc	_fresc
#define fress	_fress

#define besj0	_besj0
#define besj1	_besj1
#define besjn	_besjn
#define besy0	_besy0
#define besy1	_besy1
#define besyn	_besyn
#define besjnu	_besjnu
#define besynu	_besynu
#define besjnd	_besjnd
#define besynd	_besynd
#define besjnud	_besjnud
#define besynud	_besynud
#define sbesjn	_sbesjn
#define sbesyn	_sbesyn
#define sbesjnu	_sbesjnu
#define sbesynu	_sbesynu
#define cbesh	_cbesh
#define cbesj	_cbesj
#define cbesy	_cbesy

#define besi0	_besi0
#define besi0e	_besi0e
#define besi1	_besi1
#define besi1e	_besi1e
#define besin	_besin
#define besk0	_besk0
#define besk0e	_besk0e
#define besk1	_besk1
#define besk1e	_besk1e
#define beskn	_beskn
#define besinu	_besinu
#define besknu	_besknu
#define besind	_besind
#define besknd	_besknd
#define besinud	_besinud
#define besknud	_besknud
#define sbesin	_sbesin
#define sbeskn	_sbeskn
#define sbesinu	_sbesinu
#define sbesknu	_sbesknu
#define cbesi	_cbesi
#define cbesk	_cbesk

#define airyai	_airyai
#define airybi	_airybi
#define airyaid	_airyaid
#define airybid	_airybid
#define cairy	_cairy
#define cbiry	_cbiry

#define chu	_chu
#define hyp1f1	_hyp1f1
#define hyp1f1r	_hyp1f1r
#define lhyp1f1	_lhyp1f1
#define hyp2f1	_hyp2f1
#define hyp0f1	_hyp0f1
#define hyp1f0	_hyp1f0
#define hyp2f0	_hyp2f0
#define hyppfq	_hyppfq

#define jelli	_jelli
#define jsn	_jsn
#define jcn	_jcn
#define jdn	_jdn
#define jns	_jns
#define jnc	_jnc
#define jnd	_jnd
#define jsc	_jsc
#define jsd	_jsd
#define jdc	_jdc
#define jds	_jds
#define jcs	_jcs
#define jcd	_jcd
#define jtheta1	_jtheta1
#define jtheta1t	_jtheta1t
#define jtheta2	_jtheta2
#define jtheta2t	_jtheta2t
#define jtheta3	_jtheta3
#define jtheta3t	_jtheta3t
#define jtheta3m1	_jtheta3m1
#define jtheta3m1t	_jtheta3m1t
#define jtheta4	_jtheta4
#define jtheta4t	_jtheta4t
#define jtheta4m1	_jtheta4m1
#define jtheta4m1t	_jtheta4m1t

#define celli1	_celli1
#define celli2	_celli2
#define celli3	_celli3
#define elli1	_elli1
#define elli2	_elli2
#define elli3	_elli3
#define rc	_rc
#define rd	_rd
#define rg	_rg
#define rf	_rf
#define rj	_rj
#define jzeta	_jzeta
#define hlambda	_hlambda

#define dfzero	_dfzero
#define dfzero_r	_dfzero_r

#define cpqr79	_cpqr79
#define rpqr79	_rpqr79
#define cpzero	_cpzero
#define rpzero	_rpzero
#define rpzero2	_rpzero2
#define dka	_dka

#define dfzero	_dfzero
#define dfzero_r_	_dfzero_r

#define chkder	_chkder
#define hybrd	_hybrd
#define hybrd1	_hybrd1
#define hybrd1_r	_hybrd1_r
#define hybrd_r	_hybrd_r
#define hybrj	_hybrj
#define hybrj1	_hybrj1
#define hybrj1_r	_hybrj1_r
#define hybrj_r	_hybrj_r
#define sos	_sos
#define sos_r	_sos_r

#define dfmin	_dfmin
#define dfmin_r	_dfmin_r

#define ivset	_ivset
#define mnf	_mnf
#define mnf_r	_mnf_r
#define mng	_mng
#define mng_r	_mng_r
#define mnh	_mnh
#define mnh_r	_mnh_r
#define optif0	_optif0
#define optif0_r	_optif0_r
#define optif9	_optif9
#define optif9_r	_optif9_r
#define subplex	_subplex
#define subplex_r	_subplex_r

#define mnfb	_mnfb
#define mnfb_r	_mnfb_r
#define mngb	_mngb
#define mngb_r	_mngb_r
#define mnhb	_mnhb
#define mnhb_r	_mnhb_r

#define rfft1b	_rfft1b
#define rfft1f	_rfft1f
#define rfft1i	_rfft1i
#define rfftmb	_rfftmb
#define rfftmf	_rfftmf
#define rfftmi	_rfftmi

#define cfft1b	_cfft1b
#define cfft1f	_cfft1f
#define cfft1i	_cfft1i
#define cfftmb	_cfftmb
#define cfftmf	_cfftmf
#define cfftmi	_cfftmi

#define cost1b	_cost1b
#define cost1f	_cost1f
#define cost1i	_cost1i
#define costmb	_costmb
#define costmf	_costmf
#define costmi	_costmi
#define sint1b	_sint1b
#define sint1f	_sint1f
#define sint1i	_sint1i
#define sintmb	_sintmb
#define sintmf	_sintmf
#define sintmi	_sintmi

#define cosq1b	_cosq1b
#define cosq1f	_cosq1f
#define cosq1i	_cosq1i
#define cosqmb	_cosqmb
#define cosqmf	_cosqmf
#define cosqmi	_cosqmi
#define sinq1b	_sinq1b
#define sinq1f	_sinq1f
#define sinq1i	_sinq1i
#define sinqmb	_sinqmb
#define sinqmf	_sinqmf
#define sinqmi	_sinqmi

#define cfft2b	_cfft2b
#define cfft2f	_cfft2f
#define cfft2i	_cfft2i
#define rfft2b	_rfft2b
#define rfft2c	_rfft2c
#define rfft2f	_rfft2f
#define rfft2i	_rfft2i

#define covar	_covar
#define lmder	_lmder
#define lmder1	_lmder1
#define lmder1_r	_lmder1_r
#define lmder_r	_lmder_r
#define lmdif	_lmdif
#define lmdif1	_lmdif1
#define lmdif1_r	_lmdif1_r
#define lmdif_r	_lmdif_r
#define lmstr	_lmstr
#define lmstr1	_lmstr1
#define lmstr1_r	_lmstr1_r
#define lmstr_r	_lmstr_r
#define n2f	_n2f
#define n2f_r	_n2f_r
#define n2g	_n2g
#define n2g_r	_n2g_r
#define n2p	_n2p
#define n2p_r	_n2p_r

#define n2fb	_n2fb
#define n2fb_r	_n2fb_r
#define n2gb	_n2gb
#define n2gb_r	_n2gb_r
#define n2pb	_n2pb
#define n2pb_r	_n2pb_r

/* M4 */

#define banfac	_banfac
#define banslv	_banslv
#define bint4	_bint4
#define bintk	_bintk
#define bspldr	_bspldr
#define bsplev	_bsplev
#define bsplpp	_bsplpp
#define bsplvd	_bsplvd
#define bsplvn	_bsplvn
#define bvalue	_bvalue
#define chfdv	_chfdv
#define chfev	_chfev
#define fitlag	_fitlag
#define interv	_interv
#define pchbs	_pchbs
#define pchcm	_pchcm
#define pchfd	_pchfd
#define pchfe	_pchfe
#define pchic	_pchic
#define pchim	_pchim
#define pchse	_pchse
#define pchsp	_pchsp
#define polcof	_polcof
#define polint	_polint
#define polyvl	_polyvl
#define ppvalu	_ppvalu

#define bfqad	_bfqad
#define bfqad_r	_bfqad_r
#define bsqad	_bsqad
#define pchia	_pchia
#define pchid	_pchid
#define pfqad	_pfqad
#define pfqad_r	_pfqad_r
#define ppqad	_ppqad

#define defint	_defint
#define defint_r	_defint_r
#define qag	_qag
#define qags	_qags
#define qags_r	_qags_r
#define qag_r	_qag_r
#define qk15	_qk15
#define qk15_r	_qk15_r
#define qk21	_qk21
#define qk21_r	_qk21_r
#define qk31	_qk31
#define qk31_r	_qk31_r
#define qk41	_qk41
#define qk41_r	_qk41_r
#define qk51	_qk51
#define qk51_r	_qk51_r
#define qk61	_qk61
#define qk61_r	_qk61_r
#define qng	_qng
#define qng_r	_qng_r

#define avint	_avint

#define qagp	_qagp
#define qagp_r	_qagp_r
#define qawc	_qawc
#define qawc_r	_qawc_r
#define qawo	_qawo
#define qawo_r	_qawo_r
#define qaws	_qaws
#define qaws_r	_qaws_r

#define dehint	_dehint
#define dehint_r	_dehint_r
#define deoint	_deoint
#define deoint_r	_deoint_r
#define qagi	_qagi
#define qagi_r	_qagi_r
#define qawf	_qawf
#define qawf_r	_qawf_r
#define qk15i	_qk15i
#define qk15i_r	_qk15i_r

#define deiint	_deiint
#define deiint_r	_deiint_r

#define contd5	_contd5
#define contd8	_contd8
#define contx1	_contx1
#define contx2	_contx2
#define deabm	_deabm
#define deabm_r	_deabm_r
#define derkf	_derkf
#define derkf_int	_derkf_int
#define derkf_r	_derkf_r
#define dop853	_dop853
#define dop853_r	_dop853_r
#define dopri5	_dopri5
#define dopri5_r	_dopri5_r
#define doprin	_doprin
#define doprin_r	_doprin_r
#define dverk	_dverk
#define dverk_int	_dverk_int
#define dverk_r	_dverk_r
#define odex	_odex
#define odex2	_odex2
#define odex2_r	_odex2_r
#define odex_r	_odex_r
#define retard	_retard
#define retard_r	_retard_r
#define ylag	_ylag
#define ylag_r	_ylag_r

#define contex	_contex
#define contr5	_contr5
#define contra	_contra
#define contro	_contro
#define contrp	_contrp
#define dassl	_dassl
#define dassl_r	_dassl_r
#define debdf	_debdf
#define debdf_r	_debdf_r
#define radau	_radau
#define radau5	_radau5
#define radau5_r	_radau5_r
#define radaup	_radaup
#define radaup_r	_radaup_r
#define radau_r	_radau_r
#define rodas	_rodas
#define rodas_r	_rodas_r
#define seulex	_seulex
#define seulex_r	_seulex_r

#define derkfa	_derkfa
#define derkfa_r	_derkfa_r
#define dopri5a	_dopri5a
#define dopri5a_r	_dopri5a_r
#define dverka	_dverka
#define dverka_r	_dverka_r
#define dop853a	_dop853a
#define dop853a_r	_dop853a_r
#define odexa	_odexa
#define odexa_r	_odexa_r
#define doprina	_doprina
#define doprina_r	_doprina_r
#define dopn43	_dopn43
#define dopn43_r	_dopn43_r
#define dopn64	_dopn64
#define dopn64_r	_dopn64_r
#define dopn86	_dopn86
#define dopn86_r	_dopn86_r
#define dopn1210	_dopn1210
#define dopn1210_r	_dopn1210_r
#define odex2a	_odex2a
#define odex2a_r	_odex2a_r
#define retarda	_retarda
#define retarda_r	_retarda_r
#define ylaga	_ylaga

#define radaua	_radaua
#define radaua_r	_radaua_r
#define rodasa	_rodasa
#define rodasa_r	_rodasa_r
#define seulexa	_seulexa
#define seulexa_r	_seulexa_r

#define init_zigexp	_init_zigexp
#define init_zigexp_r	_init_zigexp_r
#define zigexp	_zigexp
#define zigexp_r	_zigexp_r

#define rgama	_rgama
#define rgama_r	_rgama_r

#define init_zignorm	_init_zignorm
#define init_zignorm_r	_init_zignorm_r
#define zignorm	_zignorm
#define zignorm_r	_zignorm_r

#define seed48	_seed48
#define srand48	_srand48
#define lcong48	_lcong48
#define drand48	_drand48
#define erand48	_erand48
#define jrand48	_jrand48
#define lrand48	_lrand48
#define mrand48	_mrand48
#define nrand48	_nrand48

#define init_genrand	_init_genrand
#define init_by_array	_init_by_array
#define genrand_int31	_genrand_int31
#define genrand_int32	_genrand_int32
#define genrand_real1	_genrand_real1
#define genrand_real2	_genrand_real2
#define genrand_real3	_genrand_real3
#define genrand_res53	_genrand_res53

#define init_genrand64	_init_genrand64
#define init_by_array64	_init_by_array64
#define genrand64_int63	_genrand64_int63
#define genrand64_int64	_genrand64_int64
#define genrand64_real1	_genrand64_real1
#define genrand64_real2	_genrand64_real2
#define genrand64_real3	_genrand64_real3

#define ranf_array	_ranf_array
#define ranf_arr_next	_ranf_arr_next
#define ranf_start	_ranf_start
#define ran_array	_ran_array
#define ran_arr_next	_ran_arr_next
#define ran_start	_ran_start

#define zigexp	_zigexp
#define init_zigexp	_init_zigexp
#define zigexp_r	_zigexp_r
#define init_zigexp_r	_init_zigexp_r

#define rgama	_rgama
#define rgama_r	_rgama_r

#define zignorm	_zignorm
#define init_zignorm	_init_zignorm
#define zignorm_r	_zignorm_r
#define init_zignorm_r	_init_zignorm_r
