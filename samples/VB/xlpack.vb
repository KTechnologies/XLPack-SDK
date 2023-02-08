'*****************************************
'*                                       *
'*  Experimental VB interface to XLPack  *
'*  Version 6.1 (December 1, 2022)       *
'*  (C) 2014-2022  K Technologies        *
'*                                       *
'*****************************************

Imports System.Math
Imports System.Runtime.InteropServices

Class XLPack

#If macOS
    Const DLL = "/Library/Application Support/Microsoft/Office365/User Content.localized/Add-Ins.localized/XLPack.dylib"
#ElseIf Win64 Then
    Const DLL = "XLPack.dll"
#Else
    Const DLL = "XLPack_32.dll"
#End If

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_d1num")>
Public Shared Function D1num(I As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_factorial")>
Public Shared Function Factorial(IX As UInteger, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_li")>
Public Shared Function Li(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_ei")>
Public Shared Function Ei(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_e1")>
Public Shared Function E1(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_digamma")>
Public Shared Function Digamma(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besj0")>
Public Shared Function Besj0(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besj1")>
Public Shared Function Besj1(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besjnu")>
Public Shared Function Besjnu(Nu As Double, X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besy0")>
Public Shared Function Besy0(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besy1")>
Public Shared Function Besy1(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besynu")>
Public Shared Function Besynu(Nu As Double, X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besi0")>
Private Shared Function _Besi0(X As Double, ByRef Info As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besi0e")>
Private Shared Function _Besi0e(X As Double, ByRef Info As Integer) As Double
End Function

Public Shared Function Besi0(X As Double, ByRef Info As Integer, Optional Kode As Integer = 1) As Double
    If Kode = 2 Then
        Besi0 = _Besi0e(X, Info)
    Else
        Besi0 = _Besi0(X, Info)
    End If
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besi1")>
Private Shared Function _Besi1(X As Double, ByRef Info As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besi1e")>
Private Shared Function _Besi1e(X As Double, ByRef Info As Integer) As Double
End Function

Public Shared Function Besi1(X As Double, ByRef Info As Integer, Optional Kode As Integer = 1) As Double
    If Kode = 2 Then
        Besi1 = _Besi1e(X, Info)
    Else
        Besi1 = _Besi1(X, Info)
    End If
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besinu")>
Public Shared Function Besinu(Nu As Double, X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besk0")>
Private Shared Function _Besk0(X As Double, ByRef Info As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besk0")>
Private Shared Function _Besk0e(X As Double, ByRef Info As Integer) As Double
End Function

Public Shared Function Besk0(X As Double, ByRef Info As Integer, Optional Kode As Integer = 1) As Double
    If Kode = 2 Then
        Besk0 = _Besk0e(X, Info)
    Else
        Besk0 = _Besk0(X, Info)
    End If
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besk1")>
Private Shared Function _Besk1(X As Double, ByRef Info As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besk1")>
Private Shared Function _Besk1e(X As Double, ByRef Info As Integer) As Double
End Function

Public Shared Function Besk1(X As Double, ByRef Info As Integer, Optional Kode As Integer = 1) As Double
    If Kode = 2 Then
        Besk1 = _Besk1e(X, Info)
    Else
        Besk1 = _Besk1(X, Info)
    End If
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_besknu")>
Public Shared Function Besknu(Nu As Double, X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_celli1")>
Public Shared Function Celli1(K As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_celli2")>
Public Shared Function Celli2(K As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="x_celli3")>
Public Shared Function Celli3(N As Double, K As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dconst")>
Public Shared Function Dconst(I As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dlange")>
Private Shared Function _Dlange(Norm As Char, M As Integer, N As Integer, Lda As Integer, A(,) As Double, Work() As Double) As Double
End Function

Public Shared Function Dlange(Norm As Char, M As Integer, N As Integer, A(,) As Double) As Double
    Dim Lda As Integer
    Dim Work() As Double
    Lda = A.GetLength(1)
    Redim Work(Max(0, M - 1))
    Dlange = _Dlange(Norm, M, N, Lda, A, Work)
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dlansy")>
Private Shared Function _Dlansy(Norm As Char, Uplo As Char, N As Integer, Lda As Integer, A(,) As Double, Work() As Double) As Double
End Function

Public Shared Function Dlansy(Norm As Char, Uplo As Char, N As Integer, A(,) As Double) As Double
    Dim Lda As Integer
    Dim Work() As Double
    Lda = A.GetLength(1)
    Redim Work(Max(0, N - 1))
    Dlansy = _Dlansy(Norm, Uplo, N, Lda, A, Work)
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dgesv")>
Private Shared Sub _Dgesv(N As Integer, Nrhs As Integer, Lda As Integer, A(,) As Double, IPiv() As Integer, Ldb As Integer, B() As Double, ByRef Info As Integer)
End Sub

Public Shared Sub Dgesv(N As Integer, A(,) As Double, IPiv() As Integer, B() As Double, ByRef Info As Integer, Optional Nrhs As Integer = 1)
    Dim Lda As Integer, Ldb As Integer
    Lda = A.GetLength(1)
    Ldb = N
    Call _Dgesv(N, Nrhs, Lda, A, IPiv, Ldb, B, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dgecon")>
Private Shared Sub _Dgecon(Norm As Char, N As Integer, Lda As Integer, A(,) As Double, ANorm As Double, ByRef RCond As Double, Work() As Double, IWork() As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Dgecon(Norm As Char, N As Integer, A(,) As Double, ANorm As Double, ByRef RCond As Double, ByRef Info As Integer)
    Dim Lda As Integer
    Dim Work() As Double, IWork() As Integer
    Lda = A.GetLength(1)
    Redim Work(Max(0, 4 * N - 1)), IWork(Max(0, N - 1))
    Call _Dgecon(Norm, N, Lda, A, ANorm, RCond, Work, IWork, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dposv")>
Private Shared Sub _Dposv(Uplo As Char, N As Integer, Nrhs As Integer, Lda As Integer, A(,) As Double, Ldb As Integer, B() As Double, ByRef Info As Integer)
End Sub

Public Shared Sub Dposv(Uplo As Char, N As Integer, A(,) As Double, B() As Double, ByRef Info As Integer, Optional Nrhs As Integer = 1)
    Dim Lda As Integer, Ldb As Integer
    Lda = A.GetLength(1)
    Ldb = N
    Call _Dposv(Uplo, N, Nrhs, Lda, A, Ldb, B, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dpocon")>
Private Shared Sub _Dpocon(Uplo As Char, N As Integer, Lda As Integer, A(,) As Double, ANorm As Double, ByRef RCond As Double, Work() As Double, IWork() As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Dpocon(Uplo As Char, N As Integer, A(,) As Double, ANorm As Double, ByRef RCond As Double, ByRef Info As Integer)
    Dim Lda As Integer
    Dim Work() As Double, IWork() As Integer
    Lda = A.GetLength(1)
    Redim Work(Max(0, 3 * N - 1)), IWork(Max(0, N - 1))
    Call _Dpocon(Uplo, N, Lda, A, ANorm, RCond, Work, IWork, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dsyev")>
Private Shared Sub _Dsyev(Jobz As Char, Uplo As Char, N As Integer, Lda As Integer, A(,) As Double, W() As Double, Work() As Double, LWork As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Dsyev(Jobz As Char, Uplo As Char, N As Integer, A(,) As Double, W() As Double, ByRef Info As Integer)
    Dim Lda As Integer
    Dim Work() As Double, LWork As Integer
    Lda = A.GetLength(1)
    Redim Work(0)
    LWork = -1
    Call _Dsyev(Jobz, Uplo, N, Lda, A, W, Work, LWork, Info)
    if Info = 0 Then
        LWork = Work(0)
        Redim Work(LWork - 1)
        Call _Dsyev(Jobz, Uplo, N, Lda, A, W, Work, Lwork, Info)
    End If
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dgels")>
Private Shared Sub _Dgels(Trans As Char, M As Integer, N As Integer, Nrhs As Integer, Lda As Integer, A(,) As Double, Ldb As Integer, B() As Double, Work() As Double, LWork As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Dgels(Trans As Char, M As Integer, N As Integer, A(,) As Double, B() As Double, ByRef Info As Integer, Optional Nrhs As Integer = 1)
    Dim Lda As Integer, Ldb As Integer
    Dim Work() As Double, LWork As Integer
    Lda = A.GetLength(1)
    Ldb = M
    Redim Work(0)
    LWork = -1
    Call _Dgels(Trans, M, N, Nrhs, Lda, A, Ldb, B, Work, LWork, Info)
    if Info = 0 Then
        LWork = Work(0)
        Redim Work(LWork - 1)
        Call _Dgels(Trans, M, N, Nrhs, Lda, A, Ldb, B, Work, LWork, Info)
    End If
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dgecov")>
Private Shared Sub _Dgecov(Job As Integer, N As Integer, Lda As Integer, A(,) As Double, Ci() As Double, ByRef Info As Integer)
End Sub

Public Shared Sub Dgecov(Job As Integer, N As Integer, A(,) As Double, Ci() As Double, ByRef Info As Integer)
    Dim Lda As Integer
    Lda = A.GetLength(1)
    Call _Dgecov(Job, N, Lda, A, Ci, Info)
End Sub


<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_pchse")>
Private Shared Sub _Pchse(N As Integer, X() As Double, F() As Double, D() As Double, Incfd As Integer, Work() As Double, LWork As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Pchse(N As Integer, X() As Double, F() As Double, D() As Double, ByRef Info As Integer, Optional Incfd As Integer = 1)
    Dim Work() As Double, LWork As Integer
    LWork = Max(1, 2 * N)
    Redim Work(LWork - 1)
    Call _Pchse(N, X, F, D, Incfd, Work, LWork, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_pchfe")>
Private Shared Sub _Pchfe(N As Integer, X() As Double, F() As Double, D() As Double, Incfd As Integer, Skip As Integer, Ne As Integer, Xe() As Double, Fe() As Double, ByRef Info As Integer)
End Sub

Public Shared Sub Pchfe(N As Integer, X() As Double, F() As Double, D() As Double, Ne As Integer, Xe() As Double, Fe() As Double, ByRef Info As Integer, Optional Incfd As Integer = 1, Optional Skip As Boolean = False)
    Dim Skip1 As Integer
    Skip1 = 0
    If Skip Then Skip1 = 1
    Call _Pchfe(N, X, F, D, Incfd, Skip1, Ne, Xe, Fe, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_pchia")>
Private Shared Function _Pchia(N As Integer, X() As Double, F() As Double, D() As Double, Incfd As Integer, Skip As Integer, A As Double, B As Double, ByRef Info As Integer) As Double
End Function

Public Shared Function Pchia(N As Integer, X() As Double, F() As Double, D() As Double, A As Double, B As Double, ByRef Info As Integer, Optional Incfd As Integer = 1, Optional Skip As Boolean = False) As Double
    Dim Skip1 As Integer
    Skip1 = 0
    If Skip Then Skip1 = 1
    Pchia = _Pchia(N, X, F, D, Incfd, Skip1, A, B, Info)
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_rpzero2")>
Private Shared Sub _Rpzero2(N As Integer, A() As Double, Zr() As Double, Zi() As Double, Iflag As Integer, MaxIter As Integer, ByRef Iter As Integer, S() As Double, Work() As Double, ByRef Info As Integer)
End Sub

Public Shared Sub Rpzero2(N As Integer, A() As Double, Zr() As Double, Zi() As Double, S() As Double, ByRef Info As Integer, Optional ByRef Iter As Long = 0, Optional IFlag As Integer = 0, Optional MaxIter As Integer = 0)
    Dim MaxIter1 As Integer, Work() As Double, LWork As Integer
    MaxIter1 = 0
    If MaxIter <= 0 Then MaxIter1 = 25 * N
    LWork = Max(1, 8 * N + 6)
    Redim Work(LWork - 1)
    Call _Rpzero2(N, A, Zr, Zi, IFlag, MaxIter1, Iter, S, Work, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dfzero_r")>
Private Shared Sub _Dfzero_r(ByRef B As Double, ByRef C As Double, R As Double, Re As Double, Ae As Double, ByRef Info As Integer, ByRef XX As Double, YY As Double, ByRef IRev As Integer)
End Sub

Public Delegate Function DfzeroFunc(X As Double) As Double

Public Shared Sub Dfzero(F As DfzeroFunc, ByRef B As Double, ByRef C As Double, R As Double, ByRef Info As Integer, Optional Re As Double = 1.0e-10, Optional Ae As Double = 1.0e-10)
    Dim IRev As Integer, XX As Double, YY As Double
    IRev = 0
    Do
        Call _Dfzero_r(B, C, R, Re, Ae, Info, XX, YY, IRev)
        If IRev <> 0 Then YY = F(XX)
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_hybrd1_r")>
Private Shared Sub _Hybrd1_r(N As Integer, X() As Double, Fvec() As Double, Tol As Double, Work() As Double, LWork As Integer, ByRef Info As Integer, XX() As Double, YY() As Double, ByRef IRev As Integer)
End Sub

Public Delegate Sub Hybrd1Proc(N As Integer, X() As Double, Fvec() As Double, ByRef IFlag As Integer)

Public Shared Sub Hybrd1(F As Hybrd1Proc, N As Integer, X() As Double, Fvec() As Double, ByRef Info As Integer, Optional Tol As Double = 1.0e-10)
    Dim IFlag As Integer, IRev As Integer, XX() As Double, YY() As Double
    Dim Work() As Double, LWork As Integer
    If N < 1 Then
        Info = -2
        Exit Sub
    ElseIf Tol < 0 Then
        Info = -5
        Exit Sub
    End If
    LWork = N * (3 * N + 13) / 2
    Redim Work(LWork - 1), XX(N - 1), YY(N - 1)
    IRev = 0
    Do
        Call _Hybrd1_r(N, X, Fvec, Tol, Work, LWork, Info, XX, YY, IRev)
        If IRev = 1 Or IRev = 2 Then
            IFlag = 1
            Call F(N, XX, YY, IFlag)
        ElseIf IRev = 3 Or IRev = 4 Then
            IFlag = 2
            Call F(N, XX, YY, IFlag)
        End If
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dfmin_r")>
Private Shared Sub _Dfmin_r(A As Double, B As Double, Tol As Double, ByRef XX As Double, YY As Double, ByRef IRev As Integer)
End Sub

Public Delegate Function DfminFunc(X As Double) As Double

Public Shared Function Dfmin(A As Double, B As Double, F As DfminFunc, Tol As Double) As Double
    Dim IRev As Integer, XX As Double, YY As Double
    IRev = 0
    Do
        Call _Dfmin_r(A, B, Tol, XX, YY, IRev)
        If IRev <> 0 Then YY = F(XX)
    Loop While IRev <> 0
    Dfmin = XX
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_optif0_r")>
Private Shared Sub _Optif0_r(N As Integer, X() As Double, Xpls() As Double, ByRef Fpls As Double, Work() As Double, LWork As Integer, ByRef Info As Integer, XX() As Double, YY As Double, ByRef IRev As Integer)
End Sub

Public Delegate Sub Optif0Proc(N As Integer, X() As Double, ByRef Fval As Double)

Public Shared Sub Optif0(N As Integer, X() As Double, F As Optif0Proc, Xpls() As Double, ByRef Fpls As Double, ByRef Info As Integer)
    Dim IRev As Integer, XX() As Double, YY As Double
    Dim Work() As Double, LWork As Integer
    If N < 1 Then
        Info = -1
        Fpls = 0
        Exit Sub
    End If
    LWork = N * (N + 10)
    Redim Work(LWork - 1), XX(N - 1)
    IRev = 0
    Do
        Call _Optif0_r(N, X, Xpls, Fpls, Work, LWork, Info, XX, YY, IRev)
        If IRev >= 1 And IRev <= 20 Then Call F(N, XX, YY)
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_qk15_r")>
Private Shared Sub _Qk15_r(A As Double, B As Double, ByRef Result As Double, ByRef AbsErr As Double, ByRef ResAbs As Double, ByRef ResAsc As Double, ByRef XX As Double, YY As Double, ByRef IRev As Integer)
End Sub

Public Delegate Function Qk15Func(X As Double) As Double

Public Shared Sub Qk15(F As Qk15Func, A As Double, B As Double, ByRef Result As Double, ByRef AbsErr As Double)
    Dim ResAbs As Double, ResAsc As Double
    Dim IRev As Integer, XX As Double, YY As Double
    IRev = 0
    Do
        Call _Qk15_r(A, B, Result, AbsErr, ResAbs, ResAsc, XX, YY, IRev)
        If IRev >= 1 And IRev <= 5 Then YY = F(XX)
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_qag_r")>
Private Shared Sub _Qag_r(A As Double, B As Double, EpsAbs As Double, EpsRel As Double, Key As Integer, Limit As Integer, ByRef Result As Double, ByRef AbsErr As Double, ByRef Neval As Integer, ByRef Last As Integer, Work() As Double, LWork As Integer, IWork() As Integer, ByRef Info As Integer, ByRef XX As Double, YY As Double, ByRef IRev As Integer)
End Sub

Public Delegate Function QagFunc(X As Double) As Double

Public Shared Sub Qag(F As QagFunc, A As Double, B As Double, ByRef Result As Double, ByRef AbsErr As Double, ByRef Info As Integer, Optional ByRef Neval As Integer = 0, Optional ByRef Last As Integer = 0, Optional EpsAbs As Double = 1.0e-10, Optional EpsRel As Double = 1.0e-10, Optional Key As Integer = 1, Optional Limit As Integer = 100)
    Dim IRev As Integer, XX As Double, YY As Double
    Dim Work() As Double, LWork As Integer, IWork() As Integer, LIWork As Integer
    If Limit < 1 Then
        Info = -7
        Result = 0
        AbsErr = 0
        Neval = 0
        Last = 0
        Exit Sub
    End If
    LWork = 4 * Limit
    LIWork = Limit
    Redim Work(LWork - 1), IWork(LIWork - 1)
    IRev = 0
    Do
        Call _Qag_r(A, B, EpsAbs, EpsRel, Key, Limit, Result, AbsErr, Neval, Last, Work, LWork, IWork, Info, XX, YY, IRev)
        If IRev >= 1 And IRev <= 15 Then YY = F(XX)
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_qagi_r")>
Private Shared Sub _Qagi_r(Bound As Double, inf As Integer, EpsAbs As Double, EpsRel As Double, Limit As Integer, ByRef Result As Double, ByRef AbsErr As Double, ByRef Neval As Integer, ByRef Last As Integer, Work() As Double, LWork As Integer, IWork() As Integer, ByRef Info As Integer, ByRef XX As Double, YY As Double, ByRef IRev As Integer)
End Sub

Public Delegate Function QagiFunc(X As Double) As Double

Public Shared Sub Qagi(F As QagiFunc, Bound As Double, Inf As Integer, ByRef Result As Double, ByRef AbsErr As Double, ByRef Info As Integer, Optional ByRef Neval As Integer = 0, Optional ByRef Last As Integer = 0, Optional EpsAbs As Double = 1.0e-10, Optional EpsRel As Double = 1.0e-10, Optional Limit As Integer = 100)
    Dim IRev As Integer, XX As Double, YY As Double
    Dim Work() As Double, LWork As Integer, IWork() As Integer, LIWork As Integer
    If Limit < 1 Then
        Info = -6
        Result = 0
        AbsErr = 0
        Neval = 0
        Last = 0
        Exit Sub
    End If
    LWork = 4 * Limit
    LIWork = Limit
    Redim Work(LWork - 1), IWork(LIWork - 1)
    IRev = 0
    Do
        Call _Qagi_r(Bound, Inf, EpsAbs, EpsRel, Limit, Result, AbsErr, Neval, Last, Work, LWork, IWork, Info, XX, YY, IRev)
        If IRev >= 1 And IRev <= 18 Then YY = F(XX)
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_derkf_r")>
Private Shared Sub _Derkf_r(N As Integer, ByRef T As Double, Y() As Double, Tout As Double, ByRef RTol As Double, ByRef ATol As Double, ITol As Integer, Mode As Integer, Work() As Double, LWork As Integer, IWork() As Integer, LIWork As Integer, ByRef Info As Integer, ByRef TT As Double, YY() As Double, YYp() As Double, ByRef IRev As Integer)
End Sub

Public Delegate Sub DerkfProc(N As Integer, T As Double, Y() As Double, Yp() As Double)

Public Shared Sub Derkf(N As Integer, F As DerkfProc, ByRef T As Double, Y() As Double, Tout As Double, Work() As Double, IWork() As Integer, ByRef Info As Integer, Optional RTol As Double = 1.0e-10, Optional ATol As Double = 0.0, Optional Mode As Integer = 0)
    Dim ITol As Integer, LWork As Integer, LIWork As Integer
    Dim IRev As Integer, TT As Double, YY() As Double, YYp() As Double
    If N < 1 Then
        Info = -1
        Exit Sub
    End If
    ITol = 0
    LWork = Work.Length
    LIWork = IWork.Length
    Redim YY(N - 1), YYp(N - 1)
    IRev = 0
    Do
        Call _Derkf_r(N, T, Y, Tout, RTol, ATol, ITol, Mode, Work, LWork, IWork, LIWork, Info, TT, YY, YYp, IRev)
        If IRev >= 1 And IRev <= 11 Then Call F(N, TT, YY, YYp)
    Loop While IRev <> 0
    If (Info = -2 Or Info = -4) Then
        Info = Info - 1
    ElseIf (Info = -5 Or Info = -6) Then
        Info = Info - 4
    ElseIf (Info = -10) Then
        Info = -6
    ElseIf (Info < 0) Then
        Info = Info + 5
    End If
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_derkf_int")>
Private Shared Sub _Derkf_int(N As Integer, T As Double, Y() As Double, Work() As Double)
End Sub

Public Shared Sub DerkfInt(N As Integer, T As Double, Y() As Double, Work() As Double)
    If N < 1 Or Y.Length < N Or Work.Length < 9*N + 20 Then Exit Sub
    Call _Derkf_int(N, T, Y, Work)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_rfft1f")>
Private Shared Sub _Rfft1f(N As Integer, Inc As Integer, R() As Double, LR As Integer, Wsave() As Double, LWsave As Integer, Work() As Double, LWork As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Rfft1f(N As Integer, R() As Double, Wsave() As Double, ByRef Info As Integer, Optional Inc As Integer = 1)
    Dim LR As Integer, LWsave As Integer
    Dim LWork As Integer, Work() As Double
    If N < 1 Then
        Info = -1
        Exit Sub
    End If
    LR = R.Length
    LWsave = Wsave.Length
    LWork = N
    Redim Work(LWork - 1)
    Call _Rfft1f(N, Inc, R, LR, Wsave, LWsave, Work, LWork, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_rfft1b")>
Private Shared Sub _Rfft1b(N As Integer, Inc As Integer, R() As Double, LR As Integer, Wsave() As Double, LWsave As Integer, Work() As Double, LWork As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Rfft1b(N As Integer, R() As Double, Wsave() As Double, ByRef Info As Integer, Optional Inc As Integer = 1)
    Dim LR As Integer, LWsave As Integer
    Dim LWork As Integer, Work() As Double
    If N < 1 Then
        Info = -1
        Exit Sub
    End If
    LR = R.Length
    LWsave = Wsave.Length
    LWork = N
    Redim Work(LWork - 1)
    Call _Rfft1b(N, Inc, R, LR, Wsave, LWsave, Work, LWork, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_rfft1i")>
Private Shared Sub _Rfft1i(N As Integer, Wsave() As Double, LWsave As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Rfft1i(N As Integer, Wsave() As Double, ByRef Info As Integer)
    Dim LWsave As Integer
    LWsave = Wsave.Length
    Call _Rfft1i(N, Wsave, LWsave, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_lmdif1_r")>
Private Shared Sub _Lmdif1_r(M As Integer, N As Integer, X() As Double, Fvec() As Double, Tol As Double, Work() As Double, LWork As Integer, IWork() As Integer, ByRef Info As Integer, XX() As Double, YY() As Double, ByRef IRev As Integer)
End Sub

Public Delegate Sub Lmdif1Proc(M As Integer, N As Integer, X() As Double, Fvec() As Double, ByRef IFlag As Integer)

Public Shared Sub Lmdif1(F As Lmdif1Proc, M As Integer, N As Integer, X() As Double, Fvec() As Double, ByRef Info As Integer, Optional Tol As Double = 1.0e-10)
    Dim IFlag As Integer, IRev As Integer, XX() As Double, YY() As Double
    Dim Work() As Double, LWork As Integer, IWork() As Integer, LIWork As Integer
    If M < N Then
        Info = -2
        Exit Sub
    ElseIf N < 1 Then
        Info = -3
        Exit Sub
    ElseIf Tol < 0 Then
        Info = -6
        Exit Sub
    End If
    LWork = N * (M + 6) + 2 * M
    LIWork = N
    Redim XX(N - 1), YY(M - 1), Work(LWork - 1), IWork(LIWork - 1)
    IRev = 0
    Do
        Call _Lmdif1_r(M, N, X, Fvec, Tol, Work, LWork, IWork, Info, XX, YY, IRev)
        If IRev = 1 Or IRev = 2 Then
            IFlag = 1
            Call F(M, N, XX, YY, IFlag)
        ElseIf IRev = 3 Then
            IFlag = 2
            Call F(M, N, XX, YY, IFlag)
        End If
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_init_genrand")>
Public Shared Sub InitGenrand(S As UInteger)
End Sub
<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_genrand_int32")>
Public Shared Function GenrandInt32() As UInteger
End Function
<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_genrand_int31")>
Public Shared Function GenrandInt31() As Integer
End Function
<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_genrand_res53")>
Public Shared Function GenrandRes53() As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_dlamch")>
Public Shared Function Dlamch(Cmach As Char) As Double
End Function

End Class
