# ****************************************
# *                                      *
# *  XLPack Numerical Library            *
# *  Version 7.0 (January 31, 2023)      *
# *  (C) 2014-2023  K Technologies       *
# *                                      *
# ****************************************

import platform
version = platform.python_version_tuple()
if platform.system() == 'Windows' and int(version[0]) >= 3 and int(version[1]) >= 8:
    # The following statements are used to run with python 3.8 or later
    import os
    # Set XLPack.dll install directory
    os.add_dll_directory(os.path.expanduser('~/AppData/Local/Programs/XLPack'))
import struct

import numpy as np
from math import *
from XLPack import *

# Test SPUtils

def Test_csr_dense():
    a = np.array([ \
        [0.2, 0.11, 0, -0.11, -0.8], \
        [0.11, 0.93, 0.81, 0.93, -0.92], \
        [0, -0.81, 0.37, 0.81, 0], \
        [0.8, 0.92, 0, 0.9, 0.86]])
    a_val = np.array([0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86])
    a_ptr = np.array([0, 3, 7, 9, 13, 16], dtype=np.int32)
    a_ind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3], dtype=np.int32)
    a_dense = np.empty((4, 5))
    print ('** csr_dense **')
    m = 5
    n = 4
    info = csr_dense(m, n, a_val, a_ptr, a_ind, 0, a_dense)
    print('info =', info)
    if info == 0:
        for i in range(m):
            for j in range(n):
                if a_dense[j, i] != a[j, i]:
                    print('error: unmatched', i, j)

def Test_dense_csr():
    a = np.array([ \
        [0.2, 0.11, 0, -0.11, -0.8], \
        [0.11, 0.93, 0.81, 0.93, -0.92], \
        [0, -0.81, 0.37, 0.81, 0], \
        [0.8, 0.92, 0, 0.9, 0.86]])
    a_val = np.array([0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86])
    a_ptr = np.array([0, 3, 7, 9, 13, 16], dtype=np.int32)
    a_ind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3], dtype=np.int32)
    print ('** dense_csr **')
    m = 5
    n = 4
    nnz = 16
    a_val_2 = np.empty(nnz)
    a_ptr_2 = np.empty(m + 1, np.int32)
    a_ind_2 = np.empty(nnz, np.int32)
    info = dense_csr(m, n, a, nnz, a_val_2, a_ptr_2, a_ind_2, 0)
    print('info =', info)
    if info == 0:
        for i in range(m + 1):
            if a_ptr_2[i] != a_ptr[i]:
                print('error: unmatched ptr', i)
        for i in range(nnz):
            if a_val_2[i] != a_val[i]:
                print('error: unmatched val', i)
            if a_ind_2[i] != a_ind[i]:
                print('error: unmatched ind', i)

def Test_csr_coo():
    a_val = np.array([0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86])
    a_rowptr = np.array([0, 3, 7, 9, 13, 16], dtype=np.int32)
    a_colind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3], dtype=np.int32)
    a_rowind = np.array([0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4], dtype=np.int32)
    print ('** csr_coo **')
    m = 5
    n = 4
    nnz = 16
    a2_val = np.empty(nnz)
    a2_rowind = np.empty(nnz, np.int32)
    a2_colind = np.empty(nnz, np.int32)
    info = csr_coo(m, n, a_val, a_rowptr, a_colind, 0, a2_val, a2_rowind, a2_colind, 0)
    print('info =', info)
    if info == 0:
        for i in range(nnz):
            if a2_val[i] != a_val[i]:
                print('error: unmatched val', i)
            if a2_rowind[i] != a_rowind[i]:
                print('error: unmatched rowind', i)
            if a2_colind[i] != a_colind[i]:
                print('error: unmatched colind', i)

def Test_coo_csr():
    a_val = np.array([0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86])
    a_rowptr = np.array([0, 3, 7, 9, 13, 16], dtype=np.int32)
    a_colind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3], dtype=np.int32)
    a_rowind = np.array([0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 3, 3, 4, 4, 4], dtype=np.int32)
    print ('** coo_csr **')
    m = 5
    n = 4
    nnz = 16
    a2_val = np.empty(nnz)
    a2_rowptr = np.empty(m + 1, np.int32)
    a2_colind = np.empty(nnz, np.int32)
    info = coo_csr(m, n, nnz, a_val, a_rowind, a_colind, 0, a2_val, a2_rowptr, a2_colind, 0)
    print('info =', info)
    if info == 0:
        for i in range(m + 1):
            if a2_rowptr[i] != a_rowptr[i]:
                print('error: unmatched rowptr', i)
        for i in range(nnz):
            if a2_val[i] != a_val[i]:
                print('error: unmatched val', i)
            if a2_colind[i] != a_colind[i]:
                print('error: unmatched colind', i)

def Test_csr_ssr():
    a_val = np.array([2.2, -0.11, -0.8, -0.11, 2.93, 0.81, -0.92, 0.81, 2.37, -0.8, -0.92, 2.86])
    a_ptr = np.array([0, 3, 7, 9, 12], dtype=np.int32)
    a_ind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 3], dtype=np.int32)
    s_val = np.array([2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86])
    s_ptr = np.array([0, 1, 3, 5, 8], dtype=np.int32)
    s_ind = np.array([0, 0, 1, 1, 2, 0, 1, 3], dtype=np.int32)
    print ('** csr_ssr **')
    n = 4
    nnz2 = 8
    a2_val = np.empty(nnz2)
    a2_ptr = np.empty(n + 1, np.int32)
    a2_ind = np.empty(nnz2, np.int32)
    info = csr_ssr('L', n, a_val, a_ptr, a_ind, 0, a2_val, a2_ptr, a2_ind, 0)
    print('info =', info)
    if info == 0:
        for i in range(n + 1):
            if a2_ptr[i] != s_ptr[i]:
                print('error: unmatched ptr', i)
        for i in range(nnz2):
            if a2_val[i] != s_val[i]:
                print('error: unmatched val', i)
            if a2_ind[i] != s_ind[i]:
                print('error: unmatched ind', i)

def Test_ssr_csr():
    a_val = np.array([2.2, -0.11, -0.8, -0.11, 2.93, 0.81, -0.92, 0.81, 2.37, -0.8, -0.92, 2.86])
    a_ptr = np.array([0, 3, 7, 9, 12], dtype=np.int32)
    a_ind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 3], dtype=np.int32)
    s_val = np.array([2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86])
    s_ptr = np.array([0, 1, 3, 5, 8], dtype=np.int32)
    s_ind = np.array([0, 0, 1, 1, 2, 0, 1, 3], dtype=np.int32)
    print ('** ssr_csr **')
    n = 4
    nnz = 12
    a2_val = np.empty(nnz)
    a2_ptr = np.empty(n + 1, np.int32)
    a2_ind = np.empty(nnz, np.int32)
    info = ssr_csr('L', n, s_val, s_ptr, s_ind, 0, a2_val, a2_ptr, a2_ind, 0)
    print('info =', info)
    if info == 0:
        for i in range(n + 1):
            if a2_ptr[i] != a_ptr[i]:
                print('error: unmatched ptr', i)
        for i in range(nnz):
            if a2_val[i] != a_val[i]:
                print('error: unmatched val', i)
            if a2_ind[i] != a_ind[i]:
                print('error: unmatched ind', i)

def Test_csr_trans():
    a_val = np.array([0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86])
    a_ptr = np.array([0, 3, 7, 9, 13, 16], dtype=np.int32)
    a_ind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3], dtype=np.int32)
    at_val = np.array([0.2, 0.11, -0.11, -0.8, 0.11, 0.93, 0.81, 0.93, -0.92, -0.81, 0.37, 0.81, 0.8, 0.92, 0.9, 0.86])
    at_ptr = np.array([0, 4, 9, 12, 16], dtype=np.int32)
    at_ind = np.array([0, 1, 3, 4, 0, 1, 2, 3, 4, 1, 2, 3, 0, 1, 3, 4], dtype=np.int32)
    print ('** csr_trans **')
    m = 5
    n = 4
    nnz = 16
    a2_val = np.empty(nnz)
    a2_ptr = np.empty(n + 1, np.int32)
    a2_ind = np.empty(nnz, np.int32)
    info = csr_trans(m, n, a_val, a_ptr, a_ind, 0, a2_val, a2_ptr, a2_ind, 0)
    print('info =', info)
    if info == 0:
        for i in range(n + 1):
            if a2_ptr[i] != at_ptr[i]:
                print('error: unmatched ptr', i)
        for i in range(nnz):
            if a2_val[i] != at_val[i]:
                print('error: unmatched val', i)
            if a2_ind[i] != at_ind[i]:
                print('error: unmatched ind', i)

# Test SPBLAS

def Test_csr_dusmv():
    a = np.array([ \
        [0.2, 0.11, 0, -0.11, -0.8], \
        [0.11, 0.93, 0.81, 0.93, -0.92], \
        [0, -0.81, 0.37, 0.81, 0], \
        [0.8, 0.92, 0, 0.9, 0.86]])
    a_val = np.array([0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86])
    a_ptr = np.array([0, 3, 7, 9, 13, 16], dtype=np.int32)
    a_ind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3], dtype=np.int32)
    x = np.array([0.73, 0.86, 0.44, 0.87])
    y = np.array([0.28, 0.30, 0.69, 0.50, 0.15])
    y_bup = np.copy(y)
    print ('** csr_dusmv **')
    m = 5
    n = 4
    alpha = 2.0
    beta = 3.0
    incx = 1
    incy = 1
    info = csr_dusmv('N', m, n, alpha, a_val, a_ptr, a_ind, 0, x, incx, beta, y, incy)
    print('info = ', info)
    if info == 0:
        err = 0.0
        for i in range(m):
            s = 0.0
            for j in range(n):
                s = s + a[j, i]*x[incx*j]
            err = err + abs(y[incy*i] - (alpha*s + beta*y_bup[incy*i]))
        print('err =', err/m)

def Test_csr_dussv():
    a = np.array([[2.2, -0.11, 0, -0.8], [0, 2.93, 0.81, -0.92], [0, 0, 2.37, 0], [0, 0, 0, 2.86]])
    a_ssr = np.array([2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86])
    a_ptr = np.array([0, 1, 3, 5, 8], dtype=np.int32)
    a_ind = np.array([0, 0, 1, 1, 2, 0, 1, 3], dtype=np.int32)
    b = np.array([ 0.73, 0.86, 0.44, 0.87])
    x = np.copy(b)
    print ('** csr_dussv **')
    n = 4
    incx = 1
    info = csr_dussv('L', 'N', 'N', n, a_ssr, a_ptr, a_ind, 0, x, incx)
    print('info = ', info)
    if info == 0:
        err = 0.0
        for i in range(n):
            s = 0.0
            for j in range(n):
                s = s + a[j, i]*x[j]
            err = err + abs(b[i] - s)
        print('err =', err/n)

def Test_csr_dusmm():
    a = np.array([ \
        [0.2, 0.11, 0, -0.11, -0.8], \
        [0.11, 0.93, 0.81, 0.93, -0.92], \
        [0, -0.81, 0.37, 0.81, 0], \
        [0.8, 0.92, 0, 0.9, 0.86]])
    a_val = np.array([0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86])
    a_ptr = np.array([0, 3, 7, 9, 13, 16], dtype=np.int32)
    a_ind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3], dtype=np.int32)
    b = np.array([ \
        [0.73, 0.86, 0.44, 0.87], \
        [0.0057, 0.89, 0.75, 0.18], \
        [0.71, 0.46, 0.69, 0.64]])
    c = np.array([ \
        [0.28, 0.30, 0.69, 0.50, 0.15], \
        [0.99, 0.24, 0.4, 0.11, 0.48], \
        [0.0029, 0.67, 0.38, 0.46, 0.014]])
    c_bup = np.copy(c)
    print ('** csr_dusmm **')
    m = 5
    n = 4
    l = 3
    alpha = 2.0
    beta = 3.0
    info = csr_dusmm('N', 'C', m, n, l, alpha, a_val, a_ptr, a_ind, 0, b, beta, c)
    print('info = ', info)
    if info == 0:
        err = 0.0
        for k in range(l):
            for i in range(m):
                s = 0.0
                for j in range(n):
                    s = s + a[j, i]*b[k, j]
                err = err + abs(c[k, i] - (alpha*s + beta*c_bup[k, i]))
        print('err =', err/(m*n))

def Test_csr_dussm():
    a = np.array([[2.2, -0.11, 0, -0.8], [0, 2.93, 0.81, -0.92], [0, 0, 2.37, 0], [0, 0, 0, 2.86]])
    a_ssr = np.array([2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86])
    a_ptr = np.array([0, 1, 3, 5, 8], dtype=np.int32)
    a_ind = np.array([0, 0, 1, 1, 2, 0, 1, 3], dtype=np.int32)
    b = np.array([ \
        [0.73, 0.86, 0.44, 0.87], \
        [0.0057, 0.89, 0.75, 0.18], \
        [0.71, 0.46, 0.69, 0.64]])
    x = np.copy(b)
    print ('** csr_dussm **')
    n = 4
    nrhs = 3
    info = csr_dussm('L', 'N', 'N', 'C', n, nrhs, a_ssr, a_ptr, a_ind, 0, x)
    print('info = ', info)
    if info == 0:
        err = 0.0
        for k in range(nrhs):
            for i in range(n):
                s = 0.0
                for j in range(n):
                    s = s + a[j, i]*x[k, j]
                err = err + abs(b[k, i] - s)
        print('err =', err/(n*nrhs))

def Test_ssr_dusmv():
    a = np.array([ \
        [2.2, -0.11, 0, -0.8], \
        [-0.11, 2.93, 0.81, -0.92], \
        [0, 0.81, 2.37, 0], \
        [-0.8, -0.92, 0, 2.86]])
    a_ssr = np.array([2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86])
    a_ptr = np.array([0, 1, 3, 5, 8], dtype=np.int32)
    a_ind = np.array([0, 0, 1, 1, 2, 0, 1, 3], dtype=np.int32)
    x = np.array([0.73, 0.86, 0.44, 0.87])
    y = np.array([0.28, 0.30, 0.69, 0.50])
    y_bup = np.copy(y)
    print ('** ssr_dusmv **')
    n = 4
    alpha = 2.0
    beta = 3.0
    incx = 1
    incy = 1
    info = ssr_dusmv('L', n, alpha, a_ssr, a_ptr, a_ind, 0, x, incx, beta, y, incy)
    print('info = ', info)
    if info == 0:
        err = 0.0
        for i in range(n):
            s = 0.0
            for j in range(n):
                s = s + a[j, i]*x[incx*j]
            err = err + abs(y[incy*i] - (alpha*s + beta*y_bup[incy*i]))
        print('err =', err/n)

# Test iterative solvers

def Test_cg1():
    a_val = np.array([2.2, -0.11, 2.93, 0.81, 2.37, -0.8, -0.92, 2.86])
    a_ptr = np.array([0, 1, 3, 5, 8], dtype=np.int32)
    a_ind = np.array([0, 0, 1, 1, 2, 0, 1, 3], dtype=np.int32)
    b = np.array([1.433, 1.3137, 2.3799, -1.4392])
    x = np.zeros(4)
    print ('** cg1 **')
    n = 4
    info, iter, res = cg1('L', n, a_val, a_ptr, a_ind, b, x)
    print('info =', info, ', iter =', iter, ', res =', res)
    print(x)

def Test_bicg1():
    a_val = np.array([2.2, -0.11, -0.8, 0.11, 2.93, 0.81, -0.92, 0.81, 2.37, 0.8, 0.92, 2.86])
    a_ptr = np.array([0, 3, 7, 9, 12], dtype=np.int32)
    a_ind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 3], dtype=np.int32)
    b = np.array([1.433, 1.4325, 2.3799, -0.3912])
    x = np.zeros(4)
    print ('** bicg1 **')
    n = 4
    info, iter, res = bicg1(n, a_val, a_ptr, a_ind, b, x)
    print('info =', info, ', iter =', iter, ', res =', res)
    print(x)

# Test I/O routines

def Test_mmio():
    fname = 'test.mtx'
    a_val = np.array([0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86])
    a_ptr = np.array([0, 3, 7, 9, 13, 16], dtype=np.int32)
    a_ind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3], dtype=np.int32)
    b_val = np.empty(16)
    b_ptr = np.empty(6, dtype=np.int32)
    b_ind = np.empty(16, dtype=np.int32)
    print ('** mmio **')
    m = 5
    n = 4
    nnz = 16
    # Write a (CSR) to MM file
    info = mm_write1(fname, m, n, a_val, a_ptr, a_ind)
    print('mm_write1: info =', info)
    if info != 0:
        return
    # Read MM file info
    info, m1, n1, nnz1, dtype, mshape, mform = mm_read_info(fname)
    print('mm_read_info: info =', info)
    if info != 0:
        return
    print(m1, n1, nnz1, dtype, mshape, mform)
    # Read MM file
    info = mm_read1(fname, b_val, b_ptr, b_ind)
    print('mm_read1: info =', info)
    if info[0] != 0:
        return
    # Compare data
    for i in range(m + 1):
        if b_ptr[i] != a_ptr[i]:
            print('error: unmatched ptr', i)
    for i in range(nnz):
        if b_val[i] != a_val[i]:
            print('error: unmatched val', i)
        if b_ind[i] != a_ind[i]:
            print('error: unmatched ind', i)

# Test check routines

def Test_csr_check():
    a_val = np.array([0.2, 0.11, 0.8, 0.11, 0.93, -0.81, 0.92,  0.81, 0.37, -0.11, 0.93, 0.81, 0.9, -0.8, -0.92, 0.86])
    a_ptr = np.array([0, 3, 7, 9, 13, 16], dtype=np.int32)
    a_ind = np.array([0, 1, 3, 0, 1, 2, 3, 1, 2, 0, 1, 2, 3, 0, 1, 3], dtype=np.int32)
    print ('** csr_check **')
    m = 5
    n = 4
    result = np.zeros(11, dtype=np.int32)
    print('m = ', m, ', n = ', n)
    info = csr_check(m, n, a_val, a_ptr, a_ind, result)
    print('csr_check: info = ', info)
    print('  base = ', result[0])
    print('  nnz = ', result[1])
    print('  nnz_lower = ', result[2])
    print('  nnz_upper = ', result[3])
    print('  nnz_diag = ', result[4])
    print('  n_zero = ', result[5])
    print('  n_zero_diag = ', result[6])
    print('  empty_rows = ', result[7])
    print('  unsorted_rows = ', result[8])
    print('  invalid_inds = ', result[9])
    print('  symmetric = ', result[10])

# -----------------

Test_csr_dense()
Test_dense_csr()
Test_csr_coo()
Test_coo_csr()
Test_csr_ssr()
Test_ssr_csr()
Test_csr_trans()

Test_csr_dusmv()
Test_csr_dussv()
Test_csr_dusmm()
Test_csr_dussm()
Test_ssr_dusmv()

Test_cg1()
Test_bicg1()

Test_mmio()

Test_csr_check()
