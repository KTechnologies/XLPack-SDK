'*****************************************
'*                                       *
'*  Experimental VB interface to XLPack  *
'*  Test program                         *
'*  Version 5.4 (November 20, 2020)      *
'*  (C) 2014-2020  K Technologies        *
'*                                       *
'*****************************************

Imports System.Math
Imports XLPack

Class Test

Shared Sub TestSf
	Dim X As Double, X2 As Double, Y As Double, Nu As Double, IX As UInteger, Errno As Integer

	Console.WriteLine("** D1num")
	Console.WriteLine("i = 1: {0}", D1num(1))
	Console.WriteLine("i = 2: {0}", D1num(2))
	Console.WriteLine("i = 3: {0}", D1num(3))
	Console.WriteLine("i = 4: {0}", D1num(4))

	Console.WriteLine("** Factorial")
	IX = 10
	Y = Factorial(IX, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Li")
	X = 2
	Y = Li(2, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Ei")
	X = 1
	Y = Ei(X, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** E1")
	X = 1
	Y = E1(X, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Digamma")
	X = 3.5
	Y = Digamma(X, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Besj0")
	X = 1
	Y = Besj0(X, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Besj1")
	X = 1
	Y = Besj1(X, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Besjnu")
	Nu = 1
	X = 1
	Y = Besjnu(Nu, X, Errno)
	Console.WriteLine("Nu = {0}, X = {1}, Y = {2}, Errno = {3}", Nu, X, Y, Errno)

	Console.WriteLine("** Besy0")
	X = 1
	Y = Besy0(X, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Besy1")
	X = 1
	Y = Besy1(X, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Besynu")
	Nu = 1
	X = 1
	Y = Besynu(Nu, X, Errno)
	Console.WriteLine("Nu = {0}, X = {1}, Y = {2}, Errno = {3}", Nu, X, Y, Errno)

	Console.WriteLine("** Besi0")
	X = 1
	Y = Besi0(X, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)
	Y = Besi0(X, Errno, 2)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Besi1")
	X = 1
	Y = Besi1(X, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)
	Y = Besi1(X, Errno, 2)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Besinu")
	Nu = 1
	X = 1
	Y = Besinu(Nu, X, Errno)
	Console.WriteLine("Nu = {0}, X = {1}, Y = {2}, Errno = {3}", Nu, X, Y, Errno)

	Console.WriteLine("** Besk0")
	X = 1
	Y = Besk0(X, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)
	Y = Besk0(X, Errno, 2)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Besk1")
	X = 1
	Y = Besk1(X, Errno)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)
	Y = Besk1(X, Errno, 2)
	Console.WriteLine("X = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Besknu")
	Nu = 1
	X = 1
	Y = Besknu(Nu, X, Errno)
	Console.WriteLine("Nu = {0}, X = {1}, Y = {2}, Errno = {3}", Nu, X, Y, Errno)

	Console.WriteLine("** Celli1")
	X = 0.5
	Y = Celli1(X, Errno)
	Console.WriteLine("K = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Celli2")
	X = 0.5
	Y = Celli2(X, Errno)
	Console.WriteLine("K = {0}, Y = {1}, Errno = {2}", X, Y, Errno)

	Console.WriteLine("** Celli3")
	X = 0.7
	X2 = 0.5
	Y = Celli3(X, X2, Errno)
	Console.WriteLine("N = {0}, K = {1}, Y = {2}, Errno = {3}", X, X2, Y, Errno)
End Sub

Shared Sub TestDconst
	Dim X As Double, I As Integer

	Console.WriteLine("** Dconst")
	For I = 0 To 35
		X = Dconst(I)
		Console.WriteLine("{0}: {1}", I, X)
	Next
End Sub

Shared Sub TestDgesv
	Const N = 3
	Dim A(,) As Double = {
		{ 0.2, -0.32, -0.8 },
		{ -0.11, 0.81, -0.92 },
		{ -0.93, 0.37, -0.29 } }
	Dim B() As Double = { -0.3727, 0.4319, -1.4247 }
	Dim IPiv(N - 1) As Integer
	Dim ANorm As Double, RCond As Double
	Dim Info As Integer

	ANorm = Dlange("1", N, N, A)
	Call Dgesv(N, A, IPiv, B, Info)
	If Info = 0 Then Call Dgecon("1", N, A, ANorm, RCond, Info)
	Console.WriteLine("** Dgesv")
	Console.WriteLine("X = {0}  {1}  {2}", B(0), B(1), B(2))
	Console.WriteLine("RCond = {0}, Info = {1}", RCond, Info)
End Sub

Shared Sub TestDposv
	Const N = 3
	Dim A(,) As Double = {
		{ 2.2, 0.0, 0.0 },
		{ -0.11, 2.93, 0.0 },
		{ -0.32, 0.81, 2.37 } }
	Dim B() As Double = { -1.566, -2.8425, -1.1765 }
	Dim ANorm As Double, RCond As Double
	Dim Info As Integer

	ANorm = Dlansy("1", "U", N, A)
	Call Dposv("U", N, A, B, Info)
	If Info = 0 Then Call Dpocon("U", N, A, ANorm, RCond, Info)
	Console.WriteLine("** Dposv")
	Console.WriteLine("X = {0}  {1}  {2}", B(0), B(1), B(2))
	Console.WriteLine("RCond = {0}, Info = {1}", RCond, Info)
End Sub

Shared Sub TestDsyev
	Const N = 3
	Dim A(,) As Double = {
		{ 2.2, 0.0, 0.0 },
		{ -0.11, 2.93, 0.0 },
		{ -0.32, 0.81, 2.37 } }
	Dim W(N - 1) As Double
	Dim Info As Integer

	Call Dsyev("V", "U", N, A, W, Info)
	Console.WriteLine("** Dsyev")
	Console.WriteLine("Eigenvalues =")
	Console.WriteLine("{0}  {1}  {2}", W(0), W(1), W(2))
	Console.WriteLine("Eigenvectors =")
	Console.WriteLine("{0}  {1}  {2}", A(0, 0), A(1, 0), A(2, 0))
	Console.WriteLine("{0}  {1}  {2}", A(0, 1), A(1, 1), A(2, 1))
	Console.WriteLine("{0}  {1}  {2}", A(0, 2), A(1, 2), A(2, 2))
	Console.WriteLine("Info = {0}", Info)
End Sub

Shared Sub TestDgels
	Const M = 4, N = 2
	Dim X() As Double = { 0.2, 118.2, 337.4, 884.6 }
	Dim Y() As Double = { 0.1, 118.1, 338.8, 888.0 }
	Dim A(N - 1, M - 1) As Double, Ci(N - 1) As Double, S As Double
	Dim I As Integer, Info As Integer

	For I = 0 To M - 1
		A(0, I) = 1
		A(1, I) = X(I)
	Next
	Call Dgels("N", M, N, A, Y, Info)
	Console.WriteLine("** Dgels")
	Console.WriteLine("a0 = {0}, a1 = {1}", Y(0), Y(1))
	Console.WriteLine("Info = {0}", Info)
	If Info = 0 Then
		Call Dgecov(0, N, A, Ci, Info)
		S = 0
		For I = N To M - 1
			S = S + Y(I) ^ 2
		Next
		S = S / (M - N)
		Console.WriteLine("Std. dev. = {0}, {1}", Sqrt(S*Ci(0)), Sqrt(S*Ci(1)))
		Console.WriteLine("Info = {0}", Info)
	End If
End Sub

Shared Sub TestPchse
	Const N = 4, Ne = 2
	Dim X() As Double = { 0.1, 0.11, 0.12, 0.13 }
	Dim Y() As Double = { 2.3026, 2.2073, 2.1203, 2.0402 }
	Dim D(N - 1) As Double, Xe(Ne - 1) As Double, Ye(Ne - 1) As Double
	Dim Info As Integer

	Call Pchse(N, X, Y, D, Info)
	Console.WriteLine("** Pchse")
	Console.WriteLine("Info = {0}", Info)
	If Info = 0 Then
		Xe(0) = 0.115: Xe(1) = 0.125
		Call Pchfe(N, X, Y, D, Ne, Xe, Ye, Info)
		Console.WriteLine("ln({0}) = {1}", Xe(0), Ye(0))
		Console.WriteLine("ln({0}) = {1}", Xe(1), Ye(1))
		Console.WriteLine("Info = {0}", Info)
	End If
End Sub

Shared Sub TestPchia
	Const N = 7
	Dim X() As Double = { -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 }
	Dim Y() As Double = { 0.5, 1.0, 0.5, 0.2, 0.1, 0.05882, 0.03846 }
	Dim D(N - 1) As Double
	Dim A As Double, B As Double, S As Double
	Dim Info As Integer

	Call Pchse(N, X, Y, D, Info)
	Console.WriteLine("** Pchse")
	Console.WriteLine("Info = {0}", Info)
	If Info = 0 Then
		A = 0: B = 4
		S = Pchia(N, X, Y, D, A, B, Info)
		Console.WriteLine("** Pchia")
		Console.WriteLine("S = {0}, S(true) = {1}", S, Atan(4.0))
		Console.WriteLine("Info = {0}", Info)
	End If
End Sub

Shared Sub TestRpzero2
	const N = 5
	Dim A() As Double = { 1.0, 0.0, 2.0, 2.0, -15.0, 10.0 }
	DIm Zr(N - 1) As Double, Zi(N - 1) As Double, S(N - 1) As Double
	Dim Iter As Integer, I As Integer, Info As Integer

	Call Rpzero2(N, A, Zr, Zi, S, Info, Iter)
	Console.WriteLine("** Rpzero2")
	For I = 0 To N - 1
		Console.WriteLine("{0,22}{1,22}{2,22}", Zr(I), Zi(I), S(I))
	Next
	Console.WriteLine("Iter = {0}, Info = {1}", Iter, Info)
End Sub

Shared Function FDfzero(X As Double) As Double
	FDfzero = (X ^ 2 - 2) * X - 5
End Function

Shared Sub TestDfzero
	Dim B As Double, C As Double, R As Double
	Dim Info As Integer

	B = 1: C = 3: R = B
	Call Dfzero(AddressOf FDfzero, B, C, R, Info)
	Console.WriteLine("** Dfzero")
	Console.WriteLine("X = {0}, Info = {1}", B, Info)
End Sub

Shared Sub FHybrd1(N As Integer, X() As Double, Fvec() As Double, ByRef IFlag As Integer)
	Fvec(0) = 4 * X(0) ^ 2 + X(1) ^ 2 - 16
	Fvec(1) = X(0) ^ 2 + X(1) ^ 2 - 9
End Sub

Shared Sub TestHybrd1
	Const N = 2
	Dim X(N - 1) As Double, Fvec(N - 1) As Double
	Dim Info As Integer

	X(0) = 1: X(1) = 2
	Call Hybrd1(AddressOf FHybrd1, N, X, Fvec, Info)
	Console.WriteLine("** Hybrd1")
	Console.WriteLine("X(0) = {0}, X(1) = {1}", X(0), X(1))
	Console.WriteLine("Info = {0}", Info)
End Sub

Shared Function FDfmin(X As Double) As Double
	FDfmin = (X ^ 2 - 2) * X - 5
End Function

Shared Sub TestDfmin
	Dim A As Double, B As Double, Tol As Double, X As Double

	A = 0: B = 1
	Tol = 1.0e-8
	X = Dfmin(A, B, AddressOf FDfzero, Tol)
	Console.WriteLine("** Dfmin")
	Console.WriteLine("X = {0}", X)
End Sub

Shared Sub FOptif0(N As Integer, X() As Double, ByRef Fval As Double)
	Fval = 100 * (X(1) - X(0) ^ 2)*(X(1) - X(0) ^ 2) + (1 - X(0)) ^ 2
End Sub

Shared Sub TestOptif0
	Const N = 2
	Dim X(N - 1) As Double, Xpls(N - 1) As Double, Fpls As Double
	Dim Info As Integer

	X(0) = -1.2: X(1) = 1
	Call Optif0(N, X, AddressOf FOptif0, Xpls, Fpls, Info)
	Console.WriteLine("** Optif0")
	Console.WriteLine("Xpls(0) = {0}, Xpls(1) = {1}", Xpls(0), Xpls(1))
	Console.WriteLine("Fpls = {0}", Fpls)
	Console.WriteLine("Info = {0}", Info)
End Sub

Shared Function FQk15(X As Double) As Double
	FQk15 = 1/(1 + X ^ 2)
End Function

Shared Sub TestQk15
	Dim A As Double, B As Double, Result As Double, AbsErr As Double

	' int(1/(1+x**2)) [0, 4] = atan(4)
	A = 0: B = 4
	Call Qk15(AddressOf FQk15, A, B, Result, AbsErr)
	Console.WriteLine("** Qk15")
	Console.WriteLine("Result = {0}, AbsErr = {1}", Result, AbsErr)
End Sub

Shared Function FQag(X As Double) As Double
	FQag = 1/(1 + X ^ 2)
End Function

Shared Sub TestQag
	Dim A As Double, B As Double
	Dim Result As Double, AbsErr As Double
	Dim Key As Integer, Neval As Integer, Last As Integer, Info As Integer

	' int(1/(1+x**2)) [0, 4] = atan(4)
	A = 0: B = 4
	For I = 1 To 6
		Key = I
		Call Qag(AddressOf FQag, A, B, Result, AbsErr, Info, Neval, Last, , , Key)
		Console.WriteLine("** Qag (Key = {0})", Key)
		Console.WriteLine("Result = {0}, AbsErr = {1}, Neval = {2}, Last = {3}", Result, AbsErr, Neval, Last)
	Next
End Sub

Shared Function FQagi(X As Double) As Double
	FQagi = 1/(1 + X ^ 2)
End Function

Shared Sub TestQagi
	Dim Bound As Double, Inf As Integer
	Dim Result As Double, AbsErr As Double
	Dim Neval As Integer, Last As Integer, Info As Integer

	' int(1/(1+x**2)) [-inf, +inf] = pi
	' int(1/(1+x**2)) [0, +inf] = pi / 2
	' int(1/(1+x**2)) [-inf, 0] = pi / 2
	Bound = 0
	Inf = 2
	Call Qagi(AddressOf FQagi, Bound, Inf, Result, AbsErr, Info, Neval, Last)
	Console.WriteLine("** Qagi [-inf, +inf]")
	Console.WriteLine("Result = {0}, AbsErr = {1}, Neval = {2}, Last = {3}", Result, AbsErr, Neval, Last)
	Inf = 1
	Call Qagi(AddressOf FQagi, Bound, Inf, Result, AbsErr, Info, Neval, Last)
	Console.WriteLine("** Qagi [0, +inf]")
	Console.WriteLine("Result = {0}, AbsErr = {1}, Neval = {2}, Last = {3}", Result, AbsErr, Neval, Last)
	Inf = -1
	Call Qagi(AddressOf FQagi, Bound, Inf, Result, AbsErr, Info, Neval, Last)
	Console.WriteLine("** Qagi [-inf, 0]")
	Console.WriteLine("Result = {0}, AbsErr = {1}, Neval = {2}, Last = {3}", Result, AbsErr, Neval, Last)
End Sub

Shared AlfaSq As Double

Shared Sub FDerkf(N As Integer, T As Double, Y() As Double, Yp() As Double)
	Dim R As Double
	R = Y(0) ^ 2 + Y(1) ^ 2
	R = R * Sqrt(R) / AlfaSq
	Yp(0) = Y(2)
	Yp(1) = Y(3)
	Yp(2) = -Y(0) / R
	Yp(3) = -Y(1) / R
End Sub

Shared Sub TestDerkf
	Const N = 4, LWork = 7 * N + 20, LIWork = 20
	Dim Ecc As Double, Alfa As Double
	Dim T As Double, Tout As Double, Y(N - 1) As Double
	Dim Info As Integer
	Dim Work(LWork - 1) As Double, IWork(LIWork - 1) As Integer
	Dim TFinal As Double, TPrint As Double

	Ecc = 0.25
	Alfa = Pi / 4
	AlfaSq = Alfa ^ 2
	T = 0
	Y(0) = 1 - Ecc
	Y(1) = 0
	Y(2) = 0
	Y(3) = Alfa * Sqrt((1 + Ecc) / (1 - Ecc))
	TFinal = 12
	TPrint = 1
	Console.WriteLine("** Derkf")
	Info = 0
	Do
		Tout = T + TPrint
		Call Derkf(N, AddressOf FDerkf, T, Y, Tout, Work, IWork, Info)
		If Info <> 1 Then Exit Do
		Console.WriteLine("T = {0}, Y(0) = {1}, Y(1) = {2}, Y(2) = {3}, Y(3) = {4}", T, Y(0), Y(1), Y(2), Y(3))
	Loop While T < TFinal
	Console.WriteLine("Info = {0}", Info)
End Sub

Shared Sub TestDerkf_2
	Const N = 4, LWork = 9 * N + 20, LIWork = 20
	Dim Ecc As Double, Alfa As Double
	Dim T As Double, Tout As Double, Y(N - 1) As Double
	Dim Info As Integer
	Dim Work(LWork - 1) As Double, IWork(LIWork - 1) As Integer
	Dim TFinal As Double, TPrint As Double

	Ecc = 0.25
	Alfa = Pi / 4
	AlfaSq = Alfa ^ 2
	T = 0
	Y(0) = 1 - Ecc
	Y(1) = 0
	Y(2) = 0
	Y(3) = Alfa * Sqrt((1 + Ecc) / (1 - Ecc))
	TFinal = 12
	TPrint = 1
	Console.WriteLine("** Derkf (2) (dense output)")
	Tout = T + TPrint
	Info = 0
	Do
		Call Derkf(N, AddressOf FDerkf, T, Y, TFinal, Work, IWork, Info, Mode:=2)
		If Info = 1 Or Info = 2 Then
			While T >= Tout
				Dim Y1(N - 1) As Double
				Call DerkfInt(N, Tout, Y1, Work)
				Console.WriteLine("T = {0}, Y(0) = {1}, Y(1) = {2}, Y(2) = {3}, Y(3) = {4}", Tout, Y1(0), Y1(1), Y1(2), Y1(3))
				Tout = Tout + Tprint
			End While
		Else
			Exit Do
		End If
	Loop While T < TFinal
	Console.WriteLine("Info = {0}", Info)
End Sub

Shared Sub TestRfft1
	Dim Wsave() As Double, R() As Double, Rcopy() As Double, Diff As Double
	Dim N As Integer, LWsave As Integer
	Dim Info As Integer, I As Integer
	Dim Seed As UInteger

	' Initialization
	Console.WriteLine("** Rfft1")
	Seed = 13
	Call InitGenrand(Seed)
	N = 10
	LWsave = N + Log(N)/Log(2) + 4
	Redim Wsave(LWsave - 1)
	Call Rfft1i(N, Wsave, Info)
	If Info <> 0 Then
		Console.WriteLine("Error during initialization")
		Exit Sub
	End If
	' Generate test data
	Redim R(N - 1), Rcopy(N - 1)
	For I = 0 To N - 1
		R(I) = GenrandRes53()
		Rcopy(I) = R(I)
	Next
	' Forward transform
	Rfft1f(N, R, Wsave, Info)
	If Info <> 0 Then
		Console.WriteLine("Error in Rfft1f")
		Exit Sub
	End If
	' Backward transform
	Rfft1b(N, R, Wsave, Info)
	If Info <> 0 Then
		Console.WriteLine("Error in Rfft1b")
		Exit Sub
	End If
	' Check results
	Diff = 0
	For I = 0 To N - 1
		If Abs(R(I) - Rcopy(I)) > Diff Then
			Diff = Abs(R(I) - Rcopy(I))
		End If
		Console.WriteLine("{0,22}{1,22}{2,22}", Rcopy(I), R(I), Abs(R(I) - Rcopy(I)))
	Next
	Console.WriteLine("diff(max) = {0}", Diff)
End Sub

Shared Sub FLmdif1(M As Integer, N As Integer, X() As Double, Fvec() As Double, ByRef IFlag As Integer)
	Dim XData() As Double = { 77.6, 114.9, 141.1, 190.8, 239.9, 289.0, 332.8, 378.4, 434.8, 477.3, 536.8, 593.1, 689.1, 760.0 }
	Dim YData() As Double = { 10.07, 14.73, 17.94, 23.93, 29.61, 35.18, 40.02, 44.82, 50.76, 55.05, 61.01,  66.4, 75.47, 81.78 }
	Dim I As Integer
	For I = 0 To M - 1
		Fvec(I) = YData(I) - X(0) * (1 - Exp(-XData(I) * X(1)))
	Next
End Sub

Shared Sub TestLmdif1
	Const M = 14, N = 2: 
	Dim X(N - 1) As Double, Fvec(M - 1) As Double
	Dim Info As Integer

	X(0) = 500: X(1) = 0.0001
	Call Lmdif1(AddressOf FLmdif1, M, N, X, Fvec, Info)
	Console.WriteLine("** Lmdif1")
	Console.WriteLine("X(0) = {0}, X(1) = {1}", X(0), X(1))
	Console.WriteLine("Info = {0}", Info)
End Sub

Shared Sub TestRand
	Dim Seed As UInteger, R32 As UInteger, R31 As Integer, R53 As Double, I As Integer

	Seed = 11
	Console.WriteLine("** Random numbers: Seed = {0}", Seed)
	Call InitGenrand(Seed)
	For I = 0 To 9
		R32 = GenrandInt32()
		R31 = GenrandInt31()
		R53 = GenrandRes53()
		Console.WriteLine("{0,12}{1,12}{2,22}", R32, R31, R53)
	Next
End Sub

Shared Sub TestDlamch
	Console.WriteLine("** Dlamch")
	Console.WriteLine("e: {0}", Dlamch("e"))
	Console.WriteLine("s: {0}", Dlamch("s"))
	Console.WriteLine("b: {0}", Dlamch("b"))
	Console.WriteLine("p: {0}", Dlamch("p"))
	Console.WriteLine("n: {0}", Dlamch("n"))
	Console.WriteLine("r: {0}", Dlamch("r"))
	Console.WriteLine("m: {0}", Dlamch("m"))
	Console.WriteLine("u: {0}", Dlamch("u"))
	Console.WriteLine("l: {0}", Dlamch("l"))
	Console.WriteLine("o: {0}", Dlamch("o"))
End Sub

Public Shared Sub Main(ByVal args() As String)
	Call TestSf
	Call TestDconst
	Call TestDgesv
	Call TestDposv
	Call TestDsyev
	Call TestDgels
	Call TestPchse
	Call TestPchia
	Call TestRpzero2
	Call TestDfzero
	Call TestHybrd1
	Call TestDfmin
	Call TestOptif0
	Call TestQk15
	Call TestQag
	Call TestQagi
	Call TestDerkf
	Call TestDerkf_2
	Call TestRfft1
	Call TestLmdif1
	Call TestRand: Call TestRand
	Call TestDlamch
End Sub

End Class
