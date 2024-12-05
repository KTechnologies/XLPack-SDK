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

# Test fem2p

# Solve Laplace equation
#   d2u/dx2 + d2u/dy2 = 0
#   u(x, 0) = 0, u(0, y) = 0, u(x, 1) = x, u(1, y) = y (BC 1)

def Test_fem2p():
    print ('** fem2p **')
    # Generate mesh data
    nx = 8
    ny = 8
    n = (nx + 1)*(ny + 1)
    ne = 2*nx*ny
    nb = 2*(nx + ny)
    x = np.empty(n)
    y = np.empty(n)
    knc = np.empty((ne, 4), dtype=np.int32)
    ks = np.empty((nb, 3), dtype=np.int32)
    lb = np.empty(nb, dtype=np.int32)
    info = mesh23(nx, ny, x, y, knc, ks, lb)
    print('mesh23 end:', info[0], info[1], info[2], info[3])
    if info[0] != 0:
        return
    # Set BC 1
    nns = 2
    ib1 = np.empty(nb, dtype=np.int32)
    bv1 = np.empty(nb)
    id = np.zeros(n, dtype=np.int32)
    for i in range(nb):
        for j in range(nns):
            id[ks[i, j + 1] - 1] = lb[i]
    k = 0
    for i in range(n):
        if id[i] != 0:
            ib1[k] = i + 1
            if id[i] == 1 or id[i] == 4:
                bv1[k] = 0.0
            elif id[i] == 2:
                bv1[k] = y[i]
            elif id[i] == 3:
                bv1[k] = x[i]
            k = k + 1
    nb1 = k
    # Set BC 2
    nb2 = 0
    ks2 = np.empty((nb2, nns + 1), dtype=np.int32)
    alpha = np.empty((nb2, nns))
    beta = np.empty((nb2, nns))
    # Assemble FEM matrix
    p = np.empty(n)
    q = np.empty(n)
    f = np.empty(n)
    for i in range(n):
        p[i] = 1.0
        q[i] = 0.0
        f[i] = 0.0
    nnz = n*n
    val = np.empty(nnz)
    rowptr = np.empty(n + 1, dtype=np.int32)
    colind = np.empty(nnz, dtype=np.int32)
    b = np.empty(n)
    info = fem2p(n, ne, x, y, knc, p, q, f, nb1, ib1, bv1, nb2, ks2, alpha, beta, val, rowptr, colind, b)
    print('fem2p end: info =', info)
    if info != 0:
        return
    # Solve FEM equation by cg
    u = np.zeros(n)
    info, iter, res = cg1('F', n, val, rowptr, colind, b, u)
    print('cg1 end: info =', info, ', iter =', iter, ', res =', res)
    print(u)
    # Write to VTK file
    fname = 'test.vtk'
    z = np.zeros(n)
    info = writevtkug(fname, n, x, y, z, ne, knc, u)
    print('writevtkug end: info =', info)

# Test gmsh22 I/O

def Test_gmsh22():
    print ('** gmsh22 **')
    # Generate mesh data
    nx = 8
    ny = 8
    n = (nx + 1)*(ny + 1)
    ne = 2*nx*ny
    nb = 2*(nx + ny)
    x = np.empty(n)
    y = np.empty(n)
    ldknc = 4
    knc = np.empty((ne, ldknc), dtype=np.int32)
    ldks = 3
    ks = np.empty((nb, ldks), dtype=np.int32)
    lb = np.empty(nb, dtype=np.int32)
    info = mesh23(nx, ny, x, y, knc, ks, lb)
    print('mesh23: info =', info[0], info[1], info[2], info[3])
    if info[0] != 0:
        return
    # Set gmsh data
    ne2 = ne + nb;
    ldkc2 = 4
    kc2 = np.empty((ne2, ldkc2), dtype=np.int32)
    ldlb2 = 2
    lb2 = np.empty((ne2, ldlb2), dtype=np.int32)
    for j in range(ne):
        for i in range(ldknc):
            kc2[j, i] = knc[j, i]
        lb2[j, 0] = 1
        lb2[j, 1] = 11
    for j in range(nb):
        for i in range(ldks):
            kc2[j + ne, i] = ks[j, i]
        lb2[j + ne, 0] = 1
        lb2[j + ne, 1] = lb[j]
    z = np.zeros(n)
    # Write to file
    fname = 'test.mesh'
    info = writegmsh22(fname, n, x, y, z, ne2, kc2, lb2)
    print('writegmsh22: info =', info)
    if info != 0:
        return
    # Read a file
    x3 = np.empty(n)
    y3 = np.empty(n)
    z3 = np.empty(n)
    ldkc3 = 4
    kc3 = np.empty((ne2, ldkc3), dtype=np.int32)
    ldlb3 = 2
    lb3 = np.empty((ne2, ldlb3), dtype=np.int32)
    info = readgmsh22(fname, x3, y3, z3, kc3, lb3)
    print('readgmsh22: info =', info[0], info[1], info[2])
    if info[0] != 0:
        return
    # Compare data
    if info[1] != n or info[2] != ne2:
        print('unmatched n or ne:', info[1], info[2])
        return
    for i in range(n):
        if (x3[i] != x[i]):
            print('unmatched x data:', i)
        if (y3[i] != y[i]):
            print('unmatched y data:', i)
        if (z3[i] != z[i]):
            print('unmatched z data:', i)
    for j in range(ne2):
        ns = 3
        if (kc3[j, 0] == 2):
            ns = 4
        for i in range(ns):
            if kc3[j, i] != kc2[j, i]:
                print('unmatched connection data:', i, j)
        for i in range(lb3[j, 0]):
            if lb3[j, i] != lb2[j, i]:
                print('unmatched label data:', i, j)

# -----------------

Test_fem2p()
Test_gmsh22()
