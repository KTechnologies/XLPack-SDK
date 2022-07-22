{*********************************************
 *                                           *
 *  Experimental Delphi interface to XLPack  *
 *  Test program                             *
 *  Version 6.0 (February 12, 2022)          *
 *  (C) 2014-2022  K Technologies            *
 *                                           *
 *********************************************}
{$Ifdef FPC}
{$Mode Delphi}
{$EndIf}
{$AppType CONSOLE}
{$J+}

program test_xlpack;

uses XLPack;

{ besj0 by Chebyshev series approximation (0 <= x <= 4) }
function Bj0(x: Double; var errno: Integer): Double;
const
	n = 13;
	bj0cs: array[1..n] of Double = (
		0.100254161968939137,
		-0.665223007764405132,
		0.248983703498281314,
		-0.0332527231700357697,
		0.0023114179304694015,
		-0.0000991127741995080,
		0.0000028916708643998,
		-0.0000000612108586630,
		0.0000000009838650793,
		-0.0000000000124235515,
		0.0000000000001265433,
		-0.0000000000000010619,
		0.0000000000000000074);
begin
	Bj0 := Chebs(bj0cs, n, 0.125*x*x - 1, errno);
end;

procedure test_sf;
var
	x, x2, y, nu: Double; ix: Cardinal; errno : Integer;
begin
	Writeln('** D1num');
	Writeln('i = 1:', D1num(1));
	Writeln('i = 2:', D1num(2));
	Writeln('i = 3:', D1num(3));
	Writeln('i = 4:', D1num(4));

	Writeln('** Factorial');
	ix := 10;
	y := Factorial(ix, errno);
	Writeln('x = ', ix, ', y = ', y, ', errno = ', errno);

	Writeln('** Chebs');
	x := 1;
	y := Bj0(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Sqrt1pm1');
	x := 2.0e-13;
	y := Sqrt1pm1(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Powm1');
	x := 2.0e-13;
	x2 := 2.0;
	y := Powm1(x2, x, errno);
	Writeln('x = ', x2, x, ', y = ', y, ', errno = ', errno);

	Writeln('** Cospi');
	x := 9000.5;
	y := Cospi(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Sinpi');
	x := 9000.0;
	y := Sinpi(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Li');
	x := 2;
	y := Li(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Ei');
	x := 1;
	y := Ei(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** E1');
	x := 1;
	y := E1(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Digamma');
	x := 3.5;
	y := Digamma(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besj0');
	x := 1;
	y := Besj0(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besj1');
	x := 1;
	y := Besj1(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besjnu');
	nu := 1; x := 1;
	y := Besjnu(nu, x, errno);
	Writeln('nu = ', nu, 'x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besy0');
	x := 1;
	y := Besy0(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besy1');
	x := 1;
	y := Besy1(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besynu');
	nu := 1; x := 1;
	y := Besynu(nu, x, errno);
	Writeln('nu = ', nu, 'x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besi0');
	x := 1;
	y := Besi0(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);
	y := Besi0(x, errno, 2);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besi1');
	x := 1;
	y := Besi1(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);
	y := Besi1(x, errno, 2);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besinu');
	nu := 1; x := 1;
	y := Besinu(nu, x, errno);
	Writeln('nu = ', nu, 'x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besk0');
	x := 1;
	y := Besk0(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);
	y := Besk0(x, errno, 2);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besk1');
	x := 1;
	y := Besk1(x, errno);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);
	y := Besk1(x, errno, 2);
	Writeln('x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Besknu');
	nu := 1; x := 1;
	y := Besknu(nu, x, errno);
	Writeln('nu = ', nu, 'x = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Celli1');
	x := 0.5;
	y := Celli1(x, errno);
	Writeln('k = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Celli2');
	x := 0.5;
	y := Celli2(x, errno);
	Writeln('k = ', x, ', y = ', y, ', errno = ', errno);

	Writeln('** Celli3');
	x := 0.7; x2 := 0.5;
	y := Celli3(x, x2, errno);
	Writeln('n = ', x, 'k = ', x2, ', y = ', y, ', errno = ', errno);
end;

procedure test_dconst;
var
	x: Double; i: Integer;
begin
	Writeln('** Dconst');
	for i := 0 to 35 do
	begin
		x := Dconst(i);
		Writeln(i, '  ', x);
	end;
end;

procedure test_dgesv;
const
	n = 3; lda = n; ldb = n;
	a: array[1..n*n] of Double =
		(0.2,  -0.32, -0.8,
		-0.11, 0.81, -0.92,
		-0.93, 0.37, -0.29);
	b: array[1..n] of Double =
		(-0.3727, 0.4319, -1.4247);
var
	ipiv: array[1..n] of Integer;
	anorm, rcond: Double; info: Integer;
begin
	anorm := Dlange('1', n, n, lda, a);
	Dgesv(n, lda, a, ipiv, ldb, b, info);
	if info = 0 then
		Dgecon('1', n, lda, a, anorm, rcond, info);
	Writeln('** Dgesv');
	Writeln('x =');
	Writeln(b[1], '  ', b[2], '  ', b[3]);
	Writeln('rcond = ', rcond, ', info = ', info);
end;

procedure test_dposv;
const
	n = 3; lda = n; ldb = n;
	a: array[1..n*n] of Double =
		(2.2, 0.0, 0.0,
		-0.11, 2.93, 0.0,
		-0.32, 0.81, 2.37);
	b: array[1..n] of Double =
		(-1.566, -2.8425, -1.1765);
var
	anorm, rcond: Double; info: Integer;
begin
	anorm := Dlansy('1', 'U', n, lda, a);
	Dposv('U', n, lda, a, ldb, b, info);
	if info = 0 then
		Dpocon('U', n, lda, a, anorm, rcond, info);
	Writeln('** Dposv');
	Writeln('x =');
	Writeln(b[1], '  ', b[2], '  ', b[3]);
	Writeln('rcond = ', rcond, ', info = ', info);
end;

procedure test_dsyev;
const
	n = 3; lda = n;
	a: array[1..n*n] of Double =
		(2.2, 0.0, 0.0,
		-0.11, 2.93, 0.0,
		-0.32, 0.81, 2.37);
var
	w: array[1..n] of Double;
	info: Integer;
begin
	Dsyev('V', 'U', n, lda, a, w, info);
	Writeln('** Dsyev');
	Writeln('Eigenvalues =');
	Writeln(w[1], '  ', w[2], '  ', w[3]);
	Writeln('Eigenvectors =');
	Writeln(a[1], '  ', a[4], '  ', a[7]);
	Writeln(a[2], '  ', a[5], '  ', a[8]);
	Writeln(a[3], '  ', a[6], '  ', a[9]);
    Writeln('info = ', info);
end;

procedure test_dgels;
const
	m= 4; n = 2; lda = m; ldb = m;
	x: array[1..m] of Double =
		(0.2, 118.2, 337.4, 884.6);
	y: array[1..m] of Double =
		(0.1, 118.1, 338.8, 888.0);
var
	a: array[1..m*n] of Double; ci: array[1..n] of Double;
	s: Double; i, info: Integer;
begin
	for i := 1 to m do
	begin
		a[i] := 1.0; a[i+lda] := x[i];
	end;
	Dgels('N', m, n, lda, a, ldb, y, info);
	Writeln('** Dgels');
	Writeln('a0 = ', y[1], ', a1 = ', y[2]);
	Writeln('info = ', info);
	if info = 0 then
	begin
		Dgecov(0, n, lda, a, ci, info);
		s := 0.0;
		for i := n+1 to m do
			s := s + Sqr(y[i]);
		s := s / (m - n);
		Writeln('Std. dev. = ', Sqrt(s*ci[1]), ', ', Sqrt(s*ci[2]));
		Writeln('info = ', info);
	end;
end;

procedure test_pchse;
const
	n = 4; ne = 2;
	x: array[1..n] of Double =
		(0.1, 0.11, 0.12, 0.13);
	y: array[1..n] of Double =
		(2.3026, 2.2073, 2.1203, 2.0402);
var
	d: array[1..n] of Double;
	xe: array[1..ne] of Double; ye: array[1..ne] of Double;
	info: Integer;
begin
	Pchse(n, x, y, d, info);
	Writeln('** Pchse');
	Writeln('info = ', info);
	if info = 0 then
	begin
		xe[1] := 0.115; xe[2] := 0.125;
		Pchfe(n, x, y, d, ne, xe, ye, info);
		Writeln('ln(', xe[1], ') = ', ye[1]);
		Writeln('ln(', xe[2], ') = ', ye[2]);
		Writeln('info = ', info);
	end;
end;

procedure test_pchia;
const
	n = 7;
	x: array[1..n] of Double =
		(-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0);
	y: array[1..n] of Double =
		(0.5, 1.0, 0.5, 0.2, 0.1, 0.05882, 0.03846);
var
	d: array[1..n] of Double;
	a, b, s: Double;
	info: Integer;
begin
	Pchse(n, x, y, d, info);
	Writeln('** Pchse');
	Writeln('info = ', info);
	if info = 0 then
	begin
		a := 0.0; b := 4.0;
		s := Pchia(n, x, y, d, a, b, info);
		Writeln('** Pchia');
		Writeln('S = ', s, ', S(true) = ', ArcTan(4.0));
		Writeln('info = ', info);
	end;
end;

procedure test_rpzero2;
const
	n = 5;
	a: array[0..n] of Double =
		(1.0, 0.0, 2.0, 2.0, -15.0, 10.0);
var
	zr: array[1..n] of Double; zi: array[1..n] of Double;
	s: array[1..n] of Double;
	iter, info, i: Integer;
begin
	Rpzero2(n, a, zr, zi, iter, s, info);
	Writeln('** Rpzero2');
	for i := 1 to n do
	begin
		Writeln(zr[i], '  ', zi[i], '  ', s[i]);
	end;
	Writeln('iter = ', iter, ', info = ', info);
end;

function f_dfzero(x: Double): Double;
begin
    f_dfzero := (x*x - 2)*x - 5;
end;

procedure test_dfzero;
var
	b, c, r: Double;
	info: Integer;
begin
	b := 1.0; c := 3.0; r := b;
	Dfzero(f_dfzero, b, c, r, info);
	Writeln('** Dfzero');
	Writeln('x = ', b, ', info = ', info);
end;

procedure f_hybrd1(n: Integer; var x, fvec: array of Double; var iflag: Integer);
begin
	if (iflag = 1) or (iflag = 2) then
	begin
		fvec[0] := 4*Sqr(x[0]) + Sqr(x[1]) - 16;
		fvec[1] := Sqr(x[0]) + Sqr(x[1]) - 9;
	end;
end;

procedure test_hybrd1;
const
	n = 2;
var
	x, fvec: array[1..n] of Double;
	info: Integer;
begin
	x[1] := 1.0; x[2] := 2.0;
	Hybrd1(f_hybrd1, n, x, fvec, info);
	Writeln('** Hybrd1');
	Writeln('x[1] = ', x[1], ', x[2] = ', x[2]);
	Writeln('info = ', info);
end;

function f_dfmin(x: Double): Double;
begin
    f_dfmin := (x*x - 2)*x - 5;
end;

procedure test_dfmin;
var
	a, b, x: Double;
begin
	a := 0.0; b := 1.0;
	x := Dfmin(a, b, f_dfmin);
	Writeln('** Dfmin');
	Writeln('x = ', x);
end;

procedure f_optif0(n: Integer; var x: array of Double; var fval: Double);
begin
	{ Rosenbrock function }
	fval := 100*(x[1] - x[0]*x[0])*(x[1] - x[0]*x[0]) + (1 - x[0])*(1 - x[0]);
end;

procedure test_optif0;
const
	n = 2;
var
	x, xpls: array[1..n] of Double; fpls: Double;
	info: Integer;
begin
	x[1] := -1.2; x[2] := 1.0;
	Optif0(n, x, f_optif0, xpls, fpls, info);
	Writeln('** Optif0');
	Writeln('xpls[1] = ', xpls[1], ', xpls[2] = ', xpls[2]);
	Writeln('fpls = ', fpls);
	Writeln('info = ', info);
end;

function f_qk15(x: Double): Double;
begin
	f_qk15 := 1/(1 + x*x);
end;

procedure test_qk15;
var
	a, b, result, abserr: Double;
begin
	{ int(1/(1+x**2)) [0, 4] = atan(4) }
	a := 0; b := 4;
	Qk15(f_qk15, a, b, result, abserr);
	Writeln('** Qk15');
	Writeln('result = ', result, ', abserr = ', abserr);
end;

function f_qag(x: Double): Double;
begin
	f_qag := 1/(1 + x*x);
end;

procedure test_qag;
var
	a, b, epsabs, epsrel, result, abserr: Double;
	key, neval, last, info, i: Integer;
begin
	{ int(1/(1+x**2)) [0, 4] = atan(4) }
	a := 0; b := 4;
	epsabs := 1.0e-10; epsrel := 1.0e-10;
	for i := 1 to 6 do
	begin
		key := i;
		Qag(f_qag, a, b, result, abserr, neval, last, info, epsabs, epsrel, key);
		Writeln('** Qag (key = ', key, ')');
		Writeln('result = ', result, ', abserr = ', abserr, ', neval = ', neval, ', last = ', last);
	end;
end;

function f_qagi(x: Double): Double;
begin
	f_qagi := 1/(1 + x*x);
end;

procedure test_qagi;
var
	bound, result, abserr: Double;
	inf, neval, last, info: Integer;
begin
	{ int(1/(1+x**2)) [-inf, +inf] = pi }
	{ int(1/(1+x**2)) [0, +inf] = pi / 2 }
	{ int(1/(1+x**2)) [-inf, 0] = pi / 2 }
	bound := 0;
	inf := 2;
	Qagi(f_qagi, bound, inf, result, abserr, neval, last, info);
	Writeln('** Qagi [-inf, +inf]');
	Writeln('result = ', result, ', abserr = ', abserr, ', neval = ', neval, ', last = ', last);
	inf := 1;
	Qagi(f_qagi, bound, inf, result, abserr, neval, last, info);
	Writeln('** Qagi [0, +inf]');
	Writeln('result = ', result, ', abserr = ', abserr, ', neval = ', neval, ', last = ', last);
	inf := -1;
	Qagi(f_qagi, bound, inf, result, abserr, neval, last, info);
	Writeln('** Qagi [-inf, 0]');
	Writeln('result = ', result, ', abserr = ', abserr, ', neval = ', neval, ', last = ', last);
end;

var alfasq: Double;

procedure f_derkf(n: Integer; t: Double; var y, yp: array of Double);
var
	r: Double;
begin
	r := y[0]*y[0] + y[1]*y[1];
	r := r*Sqrt(r)/alfasq;
	yp[0] := y[2];
	yp[1] := y[3];
	yp[2] := -y[0]/r;
	yp[3] := -y[1]/r;
end;

procedure test_derkf;
const
	ecc = 0.25; alfa = pi/4;
	n = 4;
	lwork = 7*n + 20; liwork = 20;
var
	t, tout: Double;
	y: array[1..n] of Double;
	info: Integer;
	work: array[1..lwork] of Double; iwork: array[1..liwork] of Integer;
	tfinal, tprint: Double;
begin
	alfasq := alfa*alfa;
	t := 0;
	y[1] := 1 - ecc;
	y[2] := 0;
	y[3] := 0;
	y[4] := alfa*Sqrt((1 + ecc)/(1 - ecc));
	tfinal := 12;
	tprint := 1;
	Writeln('** Derkf');
	info := 0;
	repeat
		tout := t + tprint;
		Derkf(n, f_derkf, t, y, tout, work, iwork, info);
		if info <> 1 then
			break;
		Writeln('t = ', t, ', y[1] = ', y[1], ', y[2] = ', y[2], ', y[3] = ', y[3], ', y[4] = ', y[4]);
	until t >= tfinal;
	Writeln('info = ', info);
end;

procedure test_derkf_2;
const
	ecc = 0.25; alfa = pi/4;
	n = 4;
	lwork = 9*n + 20; liwork = 20;
var
	t, tout: Double;
	y: array[1..n] of Double;
	info: Integer;
	work: array[1..lwork] of Double; iwork: array[1..liwork] of Integer;
	tfinal, tprint: Double;
	y1: array[1..n] of Double;
begin
	alfasq := alfa*alfa;
	t := 0;
	y[1] := 1 - ecc;
	y[2] := 0;
	y[3] := 0;
	y[4] := alfa*Sqrt((1 + ecc)/(1 - ecc));
	tfinal := 12;
	tprint := 1;
	Writeln('** Derkf (2) (dense output)');
	tout := t + tprint;
	info := 0;
	repeat
		Derkf(n, f_derkf, t, y, tfinal, work, iwork, info, 1.0e-10, 1.0e-10, 2);
		if (info = 1) or (info = 2) then
			while t >= tout do begin
				DerkfInt(n, tout, y1, work);
				Writeln('t = ', tout, ', y[1] = ', y1[1], ', y[2] = ', y1[2], ', y[3] = ', y1[3], ', y[4] = ', y1[4]);
				tout := tout + tprint;
			end
		else
			break;
	until t >= tfinal;
	Writeln('info = ', info);
end;

procedure test_rfft1;
var
	wsave, r, rcopy: array of Double;
	diff: Double;
	n, lwsave, info, i: Integer;
	seed: Cardinal;
begin
	{ Initialization }
	Writeln('** Rfft1');
	seed := 13;
	InitGenRand(seed);
	n := 10;
	lwsave := n + Trunc(Ln(n)/Ln(2)) + 4;
	SetLength(wsave, lwsave);
	Rfft1i(n, wsave, info);
	if info <> 0 then
	begin
		Writeln('Error during initialization');
		Exit;
	end;
	{ Generate test data }
	SetLength(r, n);
	SetLength(rcopy, n);
	for i := 0 to n-1 do
	begin
		r[i] := GenRandRes53;
		rcopy[i] := r[i];
	end;
	{ Forward transform }
	Rfft1f(n, r, wsave, info);
	if info <> 0 then
	begin
		Writeln('Error in Rfft1f');
		Exit;
	end;
	{ Backward transform }
	Rfft1b(n, r, wsave, info);
	if info <> 0 then
	begin
		Writeln('Error in Rfft1b');
		Exit;
	end;
	{ Check results }
	diff := 0;
	for i := 0 to n-1 do
	begin
		if abs(r[i] - rcopy[i]) > diff then
			diff := abs(r[i] - rcopy[i]);
		Writeln(rcopy[i], '  ', r[i], '  ', abs(r[i] - rcopy[i]));
	end;
	Writeln('diff(max) = ', diff);
end;

procedure f_lmdif1(m, n: Integer; var x, fvec: array of Double; var iflag: Integer);
const
	xdata: array[1..14] of Double = (77.6, 114.9, 141.1, 190.8, 239.9,
		289.0, 332.8, 378.4, 434.8, 477.3, 536.8, 593.1, 689.1, 760.0);
	ydata: array[1..14] of Double = (10.07, 14.73, 17.94, 23.93, 29.61,
		35.18, 40.02, 44.82, 50.76, 55.05, 61.01,  66.4, 75.47, 81.78);
var
	i: Integer;
begin
	for i := 1 to m do
		fvec[i - 1] := ydata[i] - x[0]*(1 - exp(-xdata[i]*x[1]));
end;

procedure test_lmdif1;
const
	m = 14; n = 2;
var
	x: array[1..n] of Double; fvec: array[1..m] of Double;
	info: Integer;
begin
	x[1] := 500; x[2] := 0.0001;
	Lmdif1(f_lmdif1, m, n, x, fvec, info);
	Writeln('** Lmdif1');
	Writeln('x[1] = ', x[1], ', x[2] = ', x[2]);
	Writeln('info = ', info);
end;

procedure test_rand;
var
	seed, r32: Cardinal; r31: Integer; r53: Double;
	i: Integer;
begin
	seed := 11;
	Writeln('** Random numbers: seed = ', seed);
	InitGenRand(seed);
	for i := 1 to 10 do
	begin
		r32 := GenRandInt32;
		r31 := GenRandInt31;
		r53 := GenRandRes53;
		Writeln(r32, '  ', r31, '  ', r53);
	end;
end;

procedure test_dlamch;
begin
	Writeln('** Dlamch');
	Writeln('e: ', Dlamch('e'));
	Writeln('s: ', Dlamch('s'));
	Writeln('b: ', Dlamch('b'));
	Writeln('p: ', Dlamch('p'));
	Writeln('n: ', Dlamch('n'));
	Writeln('r: ', Dlamch('r'));
	Writeln('m: ', Dlamch('m'));
	Writeln('u: ', Dlamch('u'));
	Writeln('l: ', Dlamch('l'));
	Writeln('o: ', Dlamch('o'));
end;

begin
	test_sf;
	test_dconst;
	test_dgesv;
	test_dposv;
	test_dsyev;
	test_dgels;
	test_pchse;
	test_pchia;
	test_rpzero2;
	test_dfzero;
	test_hybrd1;
	test_dfmin;
	test_optif0;
	test_qk15;
	test_qag;
	test_qagi;
	test_derkf;
	test_derkf_2;
	test_rfft1;
	test_lmdif1;
	test_rand; test_rand;
	test_dlamch;
end.
