'*****************************************
'*                                       *
'*  Experimental VB interface to XLPack  *
'*  Version 6.0 (June 14, 2022)          *
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

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="d1num")>
Public Shared Function D1num(I As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_factorial")>
Public Shared Function Factorial(IX As UInteger, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_li")>
Public Shared Function Li(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_ei")>
Public Shared Function Ei(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_e1")>
Public Shared Function E1(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_digamma")>
Public Shared Function Digamma(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_besj0")>
Public Shared Function Besj0(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_besj1")>
Public Shared Function Besj1(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_besjnu")>
Public Shared Function Besjnu(Nu As Double, X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_besy0")>
Public Shared Function Besy0(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_besy1")>
Public Shared Function Besy1(X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_besynu")>
Public Shared Function Besynu(Nu As Double, X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl)>
Private Shared Function _besi0(X As Double, ByRef Info As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl)>
Private Shared Function _besi0e(X As Double, ByRef Info As Integer) As Double
End Function

Public Shared Function Besi0(X As Double, ByRef Info As Integer, Optional Kode As Integer = 1) As Double
    If Kode = 2 Then
        Besi0 = _besi0e(X, Info)
    Else
        Besi0 = _besi0(X, Info)
    End If
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl)>
Private Shared Function _besi1(X As Double, ByRef Info As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl)>
Private Shared Function _besi1e(X As Double, ByRef Info As Integer) As Double
End Function

Public Shared Function Besi1(X As Double, ByRef Info As Integer, Optional Kode As Integer = 1) As Double
    If Kode = 2 Then
        Besi1 = _besi1e(X, Info)
    Else
        Besi1 = _besi1(X, Info)
    End If
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_besinu")>
Public Shared Function Besinu(Nu As Double, X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl)>
Private Shared Function _besk0(X As Double, ByRef Info As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl)>
Private Shared Function _besk0e(X As Double, ByRef Info As Integer) As Double
End Function

Public Shared Function Besk0(X As Double, ByRef Info As Integer, Optional Kode As Integer = 1) As Double
    If Kode = 2 Then
        Besk0 = _besk0e(X, Info)
    Else
        Besk0 = _besk0(X, Info)
    End If
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl)>
Private Shared Function _besk1(X As Double, ByRef Info As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl)>
Private Shared Function _besk1e(X As Double, ByRef Info As Integer) As Double
End Function

Public Shared Function Besk1(X As Double, ByRef Info As Integer, Optional Kode As Integer = 1) As Double
    If Kode = 2 Then
        Besk1 = _besk1e(X, Info)
    Else
        Besk1 = _besk1(X, Info)
    End If
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_besknu")>
Public Shared Function Besknu(Nu As Double, X As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_celli1")>
Public Shared Function Celli1(K As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_celli2")>
Public Shared Function Celli2(K As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="_celli3")>
Public Shared Function Celli3(N As Double, K As Double, ByRef Errno As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dconst")>
Public Shared Function Dconst(I As Integer) As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dlange")>
Private Shared Function _dlange(Norm As Char, M As Integer, N As Integer, Lda As Integer, A(,) As Double, Work() As Double) As Double
End Function

Public Shared Function Dlange(Norm As Char, M As Integer, N As Integer, A(,) As Double) As Double
    Dim Lda As Integer
    Dim Work() As Double
    Lda = A.GetLength(1)
    Redim Work(Max(0, M - 1))
    Dlange = _dlange(Norm, M, N, Lda, A, Work)
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dlansy")>
Private Shared Function _dlansy(Norm As Char, Uplo As Char, N As Integer, Lda As Integer, A(,) As Double, Work() As Double) As Double
End Function

Public Shared Function Dlansy(Norm As Char, Uplo As Char, N As Integer, A(,) As Double) As Double
    Dim Lda As Integer
    Dim Work() As Double
    Lda = A.GetLength(1)
    Redim Work(Max(0, N - 1))
    Dlansy = _dlansy(Norm, Uplo, N, Lda, A, Work)
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dgesv")>
Private Shared Sub _dgesv(N As Integer, Nrhs As Integer, Lda As Integer, A(,) As Double, IPiv() As Integer, Ldb As Integer, B() As Double, ByRef Info As Integer)
End Sub

Public Shared Sub Dgesv(N As Integer, A(,) As Double, IPiv() As Integer, B() As Double, ByRef Info As Integer, Optional Nrhs As Integer = 1)
    Dim Lda As Integer, Ldb As Integer
    Lda = A.GetLength(1)
    Ldb = N
    Call _dgesv(N, Nrhs, Lda, A, IPiv, Ldb, B, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dgecon")>
Private Shared Sub _dgecon(Norm As Char, N As Integer, Lda As Integer, A(,) As Double, ANorm As Double, ByRef RCond As Double, Work() As Double, IWork() As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Dgecon(Norm As Char, N As Integer, A(,) As Double, ANorm As Double, ByRef RCond As Double, ByRef Info As Integer)
    Dim Lda As Integer
    Dim Work() As Double, IWork() As Integer
    Lda = A.GetLength(1)
    Redim Work(Max(0, 4 * N - 1)), IWork(Max(0, N - 1))
    Call _dgecon(Norm, N, Lda, A, ANorm, RCond, Work, IWork, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dposv")>
Private Shared Sub _dposv(Uplo As Char, N As Integer, Nrhs As Integer, Lda As Integer, A(,) As Double, Ldb As Integer, B() As Double, ByRef Info As Integer)
End Sub

Public Shared Sub Dposv(Uplo As Char, N As Integer, A(,) As Double, B() As Double, ByRef Info As Integer, Optional Nrhs As Integer = 1)
    Dim Lda As Integer, Ldb As Integer
    Lda = A.GetLength(1)
    Ldb = N
    Call _dposv(Uplo, N, Nrhs, Lda, A, Ldb, B, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dpocon")>
Private Shared Sub _dpocon(Uplo As Char, N As Integer, Lda As Integer, A(,) As Double, ANorm As Double, ByRef RCond As Double, Work() As Double, IWork() As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Dpocon(Uplo As Char, N As Integer, A(,) As Double, ANorm As Double, ByRef RCond As Double, ByRef Info As Integer)
    Dim Lda As Integer
    Dim Work() As Double, IWork() As Integer
    Lda = A.GetLength(1)
    Redim Work(Max(0, 3 * N - 1)), IWork(Max(0, N - 1))
    Call _dpocon(Uplo, N, Lda, A, ANorm, RCond, Work, IWork, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dsyev")>
Private Shared Sub _dsyev(Jobz As Char, Uplo As Char, N As Integer, Lda As Integer, A(,) As Double, W() As Double, Work() As Double, LWork As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Dsyev(Jobz As Char, Uplo As Char, N As Integer, A(,) As Double, W() As Double, ByRef Info As Integer)
    Dim Lda As Integer
    Dim Work() As Double, LWork As Integer
    Lda = A.GetLength(1)
    Redim Work(0)
    LWork = -1
    Call _dsyev(Jobz, Uplo, N, Lda, A, W, Work, LWork, Info)
    if Info = 0 Then
        LWork = Work(0)
        Redim Work(LWork - 1)
        Call _dsyev(Jobz, Uplo, N, Lda, A, W, Work, Lwork, Info)
    End If
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dgels")>
Private Shared Sub _dgels(Trans As Char, M As Integer, N As Integer, Nrhs As Integer, Lda As Integer, A(,) As Double, Ldb As Integer, B() As Double, Work() As Double, LWork As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Dgels(Trans As Char, M As Integer, N As Integer, A(,) As Double, B() As Double, ByRef Info As Integer, Optional Nrhs As Integer = 1)
    Dim Lda As Integer, Ldb As Integer
    Dim Work() As Double, LWork As Integer
    Lda = A.GetLength(1)
    Ldb = M
    Redim Work(0)
    LWork = -1
    Call _dgels(Trans, M, N, Nrhs, Lda, A, Ldb, B, Work, LWork, Info)
    if Info = 0 Then
        LWork = Work(0)
        Redim Work(LWork - 1)
        Call _dgels(Trans, M, N, Nrhs, Lda, A, Ldb, B, Work, LWork, Info)
    End If
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dgecov")>
Private Shared Sub _dgecov(Job As Integer, N As Integer, Lda As Integer, A(,) As Double, Ci() As Double, ByRef Info As Integer)
End Sub

Public Shared Sub Dgecov(Job As Integer, N As Integer, A(,) As Double, Ci() As Double, ByRef Info As Integer)
    Dim Lda As Integer
    Lda = A.GetLength(1)
    Call _dgecov(Job, N, Lda, A, Ci, Info)
End Sub


<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="pchse")>
Private Shared Sub _pchse(N As Integer, X() As Double, F() As Double, D() As Double, Incfd As Integer, Work() As Double, LWork As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Pchse(N As Integer, X() As Double, F() As Double, D() As Double, ByRef Info As Integer, Optional Incfd As Integer = 1)
    Dim Work() As Double, LWork As Integer
    LWork = Max(1, 2 * N)
    Redim Work(LWork - 1)
    Call _pchse(N, X, F, D, Incfd, Work, LWork, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="pchfe")>
Private Shared Sub _pchfe(N As Integer, X() As Double, F() As Double, D() As Double, Incfd As Integer, Skip As Integer, Ne As Integer, Xe() As Double, Fe() As Double, ByRef Info As Integer)
End Sub

Public Shared Sub Pchfe(N As Integer, X() As Double, F() As Double, D() As Double, Ne As Integer, Xe() As Double, Fe() As Double, ByRef Info As Integer, Optional Incfd As Integer = 1, Optional Skip As Boolean = False)
    Dim Skip1 As Integer
    Skip1 = 0
    If Skip Then Skip1 = 1
    Call _pchfe(N, X, F, D, Incfd, Skip1, Ne, Xe, Fe, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="pchia")>
Private Shared Function _pchia(N As Integer, X() As Double, F() As Double, D() As Double, Incfd As Integer, Skip As Integer, A As Double, B As Double, ByRef Info As Integer) As Double
End Function

Public Shared Function Pchia(N As Integer, X() As Double, F() As Double, D() As Double, A As Double, B As Double, ByRef Info As Integer, Optional Incfd As Integer = 1, Optional Skip As Boolean = False) As Double
    Dim Skip1 As Integer
    Skip1 = 0
    If Skip Then Skip1 = 1
    Pchia = _pchia(N, X, F, D, Incfd, Skip1, A, B, Info)
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="rpzero2")>
Private Shared Sub _rpzero2(N As Integer, A() As Double, Zr() As Double, Zi() As Double, Iflag As Integer, MaxIter As Integer, ByRef Iter As Integer, S() As Double, Work() As Double, ByRef Info As Integer)
End Sub

Public Shared Sub Rpzero2(N As Integer, A() As Double, Zr() As Double, Zi() As Double, S() As Double, ByRef Info As Integer, Optional ByRef Iter As Long = 0, Optional IFlag As Integer = 0, Optional MaxIter As Integer = 0)
    Dim MaxIter1 As Integer, Work() As Double, LWork As Integer
    MaxIter1 = 0
    If MaxIter <= 0 Then MaxIter1 = 25 * N
    LWork = Max(1, 8 * N + 6)
    Redim Work(LWork - 1)
    Call _rpzero2(N, A, Zr, Zi, IFlag, MaxIter1, Iter, S, Work, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dfzero_r")>
Private Shared Sub _dfzero_r(ByRef B As Double, ByRef C As Double, R As Double, Re As Double, Ae As Double, ByRef Info As Integer, ByRef XX As Double, YY As Double, ByRef IRev As Integer)
End Sub

Public Delegate Function DfzeroFunc(X As Double) As Double

Public Shared Sub Dfzero(F As DfzeroFunc, ByRef B As Double, ByRef C As Double, R As Double, ByRef Info As Integer, Optional Re As Double = 1.0e-10, Optional Ae As Double = 1.0e-10)
    Dim IRev As Integer, XX As Double, YY As Double
    IRev = 0
    Do
        Call _dfzero_r(B, C, R, Re, Ae, Info, XX, YY, IRev)
        If IRev <> 0 Then YY = F(XX)
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="hybrd1_r")>
Private Shared Sub _hybrd1_r(N As Integer, X() As Double, Fvec() As Double, Tol As Double, Work() As Double, LWork As Integer, ByRef Info As Integer, XX() As Double, YY() As Double, ByRef IRev As Integer)
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
        Call _hybrd1_r(N, X, Fvec, Tol, Work, LWork, Info, XX, YY, IRev)
        If IRev = 1 Or IRev = 2 Then
            IFlag = 1
            Call F(N, XX, YY, IFlag)
        ElseIf IRev = 3 Or IRev = 4 Then
            IFlag = 2
            Call F(N, XX, YY, IFlag)
        End If
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dfmin_r")>
Private Shared Sub _dfmin_r(A As Double, B As Double, Tol As Double, ByRef XX As Double, YY As Double, ByRef IRev As Integer)
End Sub

Public Delegate Function DfminFunc(X As Double) As Double

Public Shared Function Dfmin(A As Double, B As Double, F As DfminFunc, Tol As Double) As Double
    Dim IRev As Integer, XX As Double, YY As Double
    IRev = 0
    Do
        Call _dfmin_r(A, B, Tol, XX, YY, IRev)
        If IRev <> 0 Then YY = F(XX)
    Loop While IRev <> 0
    Dfmin = XX
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="optif0_r")>
Private Shared Sub _optif0_r(N As Integer, X() As Double, Xpls() As Double, ByRef Fpls As Double, Work() As Double, LWork As Integer, ByRef Info As Integer, XX() As Double, YY As Double, ByRef IRev As Integer)
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
        Call _optif0_r(N, X, Xpls, Fpls, Work, LWork, Info, XX, YY, IRev)
        If IRev >= 1 And IRev <= 20 Then Call F(N, XX, YY)
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="qk15_r")>
Private Shared Sub _qk15_r(A As Double, B As Double, ByRef Result As Double, ByRef AbsErr As Double, ByRef ResAbs As Double, ByRef ResAsc As Double, ByRef XX As Double, YY As Double, ByRef IRev As Integer)
End Sub

Public Delegate Function Qk15Func(X As Double) As Double

Public Shared Sub Qk15(F As Qk15Func, A As Double, B As Double, ByRef Result As Double, ByRef AbsErr As Double)
    Dim ResAbs As Double, ResAsc As Double
    Dim IRev As Integer, XX As Double, YY As Double
    IRev = 0
    Do
        Call _qk15_r(A, B, Result, AbsErr, ResAbs, ResAsc, XX, YY, IRev)
        If IRev >= 1 And IRev <= 5 Then YY = F(XX)
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="qag_r")>
Private Shared Sub _qag_r(A As Double, B As Double, EpsAbs As Double, EpsRel As Double, Key As Integer, Limit As Integer, ByRef Result As Double, ByRef AbsErr As Double, ByRef Neval As Integer, ByRef Last As Integer, Work() As Double, LWork As Integer, IWork() As Integer, ByRef Info As Integer, ByRef XX As Double, YY As Double, ByRef IRev As Integer)
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
        Call _qag_r(A, B, EpsAbs, EpsRel, Key, Limit, Result, AbsErr, Neval, Last, Work, LWork, IWork, Info, XX, YY, IRev)
        If IRev >= 1 And IRev <= 15 Then YY = F(XX)
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="qagi_r")>
Private Shared Sub _qagi_r(Bound As Double, inf As Integer, EpsAbs As Double, EpsRel As Double, Limit As Integer, ByRef Result As Double, ByRef AbsErr As Double, ByRef Neval As Integer, ByRef Last As Integer, Work() As Double, LWork As Integer, IWork() As Integer, ByRef Info As Integer, ByRef XX As Double, YY As Double, ByRef IRev As Integer)
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
        Call _qagi_r(Bound, Inf, EpsAbs, EpsRel, Limit, Result, AbsErr, Neval, Last, Work, LWork, IWork, Info, XX, YY, IRev)
        If IRev >= 1 And IRev <= 18 Then YY = F(XX)
    Loop While IRev <> 0
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="derkf_r")>
Private Shared Sub _derkf_r(N As Integer, ByRef T As Double, Y() As Double, Tout As Double, ByRef RTol As Double, ByRef ATol As Double, ITol As Integer, Mode As Integer, Work() As Double, LWork As Integer, IWork() As Integer, LIWork As Integer, ByRef Info As Integer, ByRef TT As Double, YY() As Double, YYp() As Double, ByRef IRev As Integer)
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
        Call _derkf_r(N, T, Y, Tout, RTol, ATol, ITol, Mode, Work, LWork, IWork, LIWork, Info, TT, YY, YYp, IRev)
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

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="derkf_int")>
Private Shared Sub _derkf_int(N As Integer, T As Double, Y() As Double, Work() As Double)
End Sub

Public Shared Sub DerkfInt(N As Integer, T As Double, Y() As Double, Work() As Double)
    If N < 1 Or Y.Length < N Or Work.Length < 9*N + 20 Then Exit Sub
    Call _derkf_int(N, T, Y, Work)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="rfft1f")>
Private Shared Sub _rfft1f(N As Integer, Inc As Integer, R() As Double, LR As Integer, Wsave() As Double, LWsave As Integer, Work() As Double, LWork As Integer, ByRef Info As Integer)
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
    Call _rfft1f(N, Inc, R, LR, Wsave, LWsave, Work, LWork, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="rfft1b")>
Private Shared Sub _rfft1b(N As Integer, Inc As Integer, R() As Double, LR As Integer, Wsave() As Double, LWsave As Integer, Work() As Double, LWork As Integer, ByRef Info As Integer)
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
    Call _rfft1b(N, Inc, R, LR, Wsave, LWsave, Work, LWork, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="rfft1i")>
Private Shared Sub _rfft1i(N As Integer, Wsave() As Double, LWsave As Integer, ByRef Info As Integer)
End Sub

Public Shared Sub Rfft1i(N As Integer, Wsave() As Double, ByRef Info As Integer)
    Dim LWsave As Integer
    LWsave = Wsave.Length
    Call _rfft1i(N, Wsave, LWsave, Info)
End Sub

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="lmdif1_r")>
Private Shared Sub _lmdif1_r(M As Integer, N As Integer, X() As Double, Fvec() As Double, Tol As Double, Work() As Double, LWork As Integer, IWork() As Integer, ByRef Info As Integer, XX() As Double, YY() As Double, ByRef IRev As Integer)
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

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="init_genrand")>
Public Shared Sub InitGenrand(S As UInteger)
End Sub
<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="genrand_int32")>
Public Shared Function GenrandInt32() As UInteger
End Function
<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="genrand_int31")>
Public Shared Function GenrandInt31() As Integer
End Function
<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="genrand_res53")>
Public Shared Function GenrandRes53() As Double
End Function

<DllImport(DLL, CallingConvention:=CallingConvention.Cdecl, EntryPoint:="dlamch")>
Public Shared Function Dlamch(Cmach As Char) As Double
End Function

End Class
