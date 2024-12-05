/****************************************
 *                                      *
 *  XLPack Numerical Library            *
 *  Version 7.0 (January 31, 2023)      *
 *  (C) 2014-2023  K Technologies       *
 *                                      *
 ****************************************/
#pragma once

/*
 * Solver
 */
#define cg	_cg
#define cg1	_cg1
#define cg_r	_cg_r
#define cr	_cr
#define cr_r	_cr_r

#define sor	_sor
#define sor_r	_sor_r

#define fom	_fom
#define fom_r	_fom_r
#define diom	_diom
#define diom_r	_diom_r
#define fgmres	_fgmres
#define fgmres_r	_fgmres_r
#define dqgmres	_dqgmres
#define dqgmres_r	_dqgmres_r
#define gcr	_gcr
#define gcr_r	_gcr_r
#define orthomin	_orthomin
#define orthomin_r	_orthomin_r
#define bicg	_bicg
#define bicg1	_bicg1
#define bicg_r	_bicg_r
#define cgs	_cgs
#define cgs_r	_cgs_r
#define gpbicg	_gpbicg
#define gpbicg_r	_gpbicg_r
#define qmr	_qmr
#define qmr_r	_qmr_r
#define tfqmr	_tfqmr
#define tfqmr_r	_tfqmr_r

#define z_cg	_z_cg
#define z_cg_r	_z_cg_r
#define z_cr	_z_cr
#define z_cr_r	_z_cr_r
#define z_cocg	_z_cocg
#define z_cocg_r	_z_cocg_r
#define z_cocr	_z_cocr
#define z_cocr_r	_z_cocr_r

#define z_sor	_z_sor
#define z_sor_r	_z_sor_r

#define z_fom	_z_fom
#define z_fom_r	_z_fom_r
#define z_diom	_z_diom
#define z_diom_r	_z_diom_r
#define z_fgmres	_z_fgmres
#define z_fgmres_r	_z_fgmres_r
#define z_dqgmres	_z_dqgmres
#define z_dqgmres_r	_z_dqgmres_r
#define z_gcr	_z_gcr
#define z_gcr_r	_z_gcr_r
#define z_orthomin	_z_orthomin
#define z_orthomin_r	_z_orthomin_r
#define z_bicg	_z_bicg
#define z_bicg_r	_z_bicg_r
#define z_cgs	_z_cgs
#define z_cgs_r	_z_cgs_r
#define z_gpbicg	_z_gpbicg
#define z_gpbicg_r	_z_gpbicg_r
#define z_qmr	_z_qmr
#define z_qmr_r	_z_qmr_r
#define z_tfqmr	_z_tfqmr
#define z_tfqmr_r	_z_tfqmr_r

/*
 * Precon
 */
#define csx_ds	_csx_ds
#define csx_ds_solve	_csx_ds_solve

#define csx_ssor	_csx_ssor
#define csc_ssor_solve	_csc_ssor_solve
#define csr_ssor_solve	_csr_ssor_solve
#define ssc_ssor_solve	_ssc_ssor_solve
#define ssr_ssor_solve	_ssr_ssor_solve

#define ssc_ic0	_ssc_ic0
#define ssr_ic0	_ssr_ic0
#define ssc_ic_solve	_ssc_ic_solve
#define ssr_ic_solve	_ssr_ic_solve

#define csc_ilu0	_csc_ilu0
#define csr_ilu0	_csr_ilu0
#define csc_ilu	_csc_ilu
#define csr_ilu	_csr_ilu
#define csc_ilu_solve	_csc_ilu_solve
#define csr_ilu_solve	_csr_ilu_solve

#define z_csx_ds	_z_csx_ds
#define z_csx_ds_solve	_z_csx_ds_solve

#define z_csx_ssor	_z_csx_ssor
#define z_csc_ssor_solve	_z_csc_ssor_solve
#define z_csr_ssor_solve	_z_csr_ssor_solve
#define z_ssc_ssor_solve	_z_ssc_ssor_solve
#define z_ssr_ssor_solve	_z_ssr_ssor_solve
#define z_hsc_ssor_solve	_z_hsc_ssor_solve
#define z_hsr_ssor_solve	_z_hsr_ssor_solve

#define z_hsc_ic0	_z_hsc_ic0
#define z_hsr_ic0	_z_hsr_ic0
#define z_hsc_ic_solve	_z_hsc_ic_solve
#define z_hsr_ic_solve	_z_hsr_ic_solve

#define z_csc_ilu0	_z_csc_ilu0
#define z_csr_ilu0	_z_csr_ilu0
#define z_csc_ilu	_z_csc_ilu
#define z_csr_ilu	_z_csr_ilu
#define z_csc_ilu_solve	_z_csc_ilu_solve
#define z_csr_ilu_solve	_z_csr_ilu_solve

/*
 * SpBlas
 */
#define csc_dusmv	_csc_dusmv
#define csr_dusmv	_csr_dusmv
#define ssc_dusmv	_ssc_dusmv
#define ssr_dusmv	_ssr_dusmv

#define csc_dussv	_csc_dussv
#define csr_dussv	_csr_dussv

#define csc_dusmm	_csc_dusmm
#define csr_dusmm	_csr_dusmm
#define csc_dussm	_csc_dussm
#define csr_dussm	_csr_dussm

#define csc_zusmv	_csc_zusmv
#define csr_zusmv	_csr_zusmv
#define hsc_zusmv	_hsc_zusmv
#define hsr_zusmv	_hsr_zusmv
#define ssc_zusmv	_ssc_zusmv
#define ssr_zusmv	_ssr_zusmv

#define csc_zussv	_csc_zussv
#define csr_zussv	_csr_zussv

#define csc_zusmm	_csc_zusmm
#define csr_zusmm	_csr_zusmm
#define csc_zussm	_csc_zussm
#define csr_zussm	_csr_zussm

/*
 * SpUtils
 */
#define coo_csc	_coo_csc
#define coo_csr	_coo_csr
#define csc_coo	_csc_coo
#define csr_coo	_csr_coo
#define csc_csr	_csc_csr
#define csr_csc	_csc_csr
#define ssc_csc	_ssc_csc
#define ssr_csr	_ssr_csr
#define csc_ssc	_csc_ssc
#define csr_ssr	_csr_ssr
#define csc_dense	_csc_dense
#define csr_dense	_csr_dense
#define coo_dense	_coo_dense
#define dense_csc	_dense_csc
#define dense_csr	_dense_csr
#define dense_coo	_dense_coo

#define csx_diag	_csx_diag
#define csx_sort	_csx_sort
#define dense_nnz	_dense_nnz

#define z_coo_csc	_z_coo_csc
#define z_coo_csr	_z_coo_csr
#define z_csc_coo	_z_csc_coo
#define z_csr_coo	_z_csr_coo
#define z_csc_csr	_z_csc_csr
#define z_csr_csc	_z_csc_csr
#define z_ssc_csc	_z_ssc_csc
#define z_ssr_csr	_z_ssr_csr
#define z_csc_ssc	_z_csc_ssc
#define z_csr_ssr	_z_csr_ssr
#define z_hsc_csc	_z_hsc_csc
#define z_hsr_csr	_z_hsr_csr
#define z_csc_hsc	_z_csc_hsc
#define z_csr_hsr	_z_csr_hsr
#define z_csc_dense	_z_csc_dense
#define z_csr_dense	_z_csr_dense
#define z_coo_dense	_z_coo_dense
#define z_dense_csc	_z_dense_csc
#define z_dense_csr	_z_dense_csr
#define z_dense_coo	_z_dense_coo

#define z_csx_diag	_z_csx_diag
#define z_csx_sort	_z_csx_sort
#define z_dense_nnz	_z_dense_nnz

#define csc_dussv_sor	_csc_dussv_sor
#define csr_dussv_sor	_csr_dussv_sor
#define csc_zussv_sor	_csc_zussv_sor
#define csr_zussv_sor	_csr_zussv_sor

#define csc_dusadd	_csc_dusadd
#define csr_dusadd	_csr_dusadd
#define csc_dzusadd	_csc_dzusadd
#define csr_dzusadd	_csr_dzusadd
#define csc_zusadd	_csc_zusadd
#define csr_zusadd	_csr_zusadd

/*
 * Input/output
 */
#define hb_read	_hb_read
#define hb_write	_hb_write
#define mm_read	_mm_read
#define mm_write	_mm_write

/*
 * Check
 */
#define csc_check	_csc_check
#define csr_check	_csr_check
#define csx_check_sym	_csx_check_sym
#define z_csc_check	_z_csc_check
#define z_csr_check	_z_csr_check
#define z_csx_check_sym	_z_csx_check_sym

/*
 * SuperLU
 */
#define dgssv	_dgssv
#define dgssvx	_dgssvx
#define zgssv	_zgssv
#define zgssvx	_zgssvx

/*
 * Arpack
 */
#define dssev	_dssev
#define dssevs	_dssevs
#define dgsev	_dgsev
#define dgsevs	_dgsevs
#define dssgv	_dssgv
#define dssgvs	_dssgvs
#define dgsgv	_dgsgv
#define dgsgvs	_dgsgvs
#define dgssvd	_dgssvd
#define zgsev	_zgsev
#define zgsevs	_zgsevs
#define zgsgv	_zgsgv
#define zgsgvs	_zgsgvs
#define zgssvd	_zgssvd
#define dsaupd	_dsaupd
#define dseupd	_dseupd
#define dnaupd	_dnaupd
#define dneupd	_dneupd
#define znaupd	_znaupd
#define zneupd	_zneupd
