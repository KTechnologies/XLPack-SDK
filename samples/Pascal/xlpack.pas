{*********************************************
 *                                           *
 *  Experimental Delphi interface to XLPack  *
 *  Version 6.1 (December 1, 2022)           *
 *  (C) 2014-2022  K Technologies            *
 *                                           *
 *********************************************}
{$Ifdef FPC}
{$Mode Delphi}
{$Warn 5036 off}
{$EndIf}

unit XLPack;

{****************}
   interface
{****************}

const
{$Ifdef Win64}
	dll = 'XLPack.dll';
{$Else}
	dll = 'XLPack_32.dll';
{$Endif}

function D1num(i: Integer): Double; cdecl; external dll name '_d1num';

function Factorial(x: Cardinal; var errno: Integer): Double; cdecl; external dll name 'x_factorial';

function Chebs(var c: array of Double; n: Integer; x: Double; var errno: Integer): Double;

function Sqrt1pm1(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_sqrt1pm1';
function Powm1(x, y: Double; var errno: Integer): Double; cdecl; external dll name 'x_powm1';
function Cospi(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_cos_pi';
function Sinpi(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_sin_pi';

function Li(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_li';
function Ei(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_ei';
function E1(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_e1';

function Digamma(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_digamma';

function Besj0(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besj0';
function Besj1(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besj1';
function Besjnu(nu: Double; x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besjnu';
function Besy0(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besy0';
function Besy1(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besy1';
function Besynu(nu: Double; x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besynu';
function Besi0(x: Double; var errno: Integer; kode: Integer = 1): Double;
function Besi1(x: Double; var errno: Integer; kode: Integer = 1): Double;
function Besinu(nu: Double; x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besinu';
function Besk0(x: Double; var errno: Integer; kode: Integer = 1): Double;
function Besk1(x: Double; var errno: Integer; kode: Integer = 1): Double;
function Besknu(nu: Double; x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besknu';

function Celli1(k: Double; var errno: Integer): Double; cdecl; external dll name 'x_celli1';
function Celli2(k: Double; var errno: Integer): Double; cdecl; external dll name 'x_celli2';
function Celli3(n: Double; k: Double; var errno: Integer): Double; cdecl; external dll name 'x_celli3';

function Dconst(i: Integer): Double; cdecl; external dll name '_dconst';

function Dlange(norm: Char; m, n, lda: Integer; var a: array of Double): Double;
function Dlansy(norm, uplo: Char; n, lda: Integer; var a: array of Double): Double;

procedure Dgesv(n, lda: Integer; var a: array of Double; var ipiv: array of Integer; ldb: Integer; var b: array of Double; var info: Integer; nrhs: Integer = 1);
procedure Dgecon(norm: Char; n, lda: Integer; var a: array of Double; anorm: Double; var rcond: Double; var info: Integer);

procedure Dposv(uplo: Char; n, lda: Integer; var a: array of Double; ldb: Integer; var b: array of Double; var info: Integer; nrhs: Integer = 1);
procedure Dpocon(uplo: Char; n, lda: Integer; var a: array of Double; anorm: Double; var rcond: Double; var info: Integer);

procedure Dsyev(jobz, uplo: Char; n, lda: Integer; var a: array of Double; var w: array of Double; var info: Integer);

procedure Dgels(trans: Char; m, n, lda: Integer; var a: array of Double; ldb: Integer; var b: array of Double; var info: Integer; nrhs: Integer = 1);
procedure Dgecov(job, n, lda: Integer; var a: array of Double; var ci: array of Double; var info: Integer);

procedure Pchse(n: Integer; var x, f, d: array of Double; var info: Integer; incfd: Integer = 1);
procedure Pchfe(n: Integer; var x, f, d: array of Double; ne: Integer; var xe, fe: array of Double; var info: Integer; incfd: Integer = 1; skip: Integer = 0);
function Pchia(n: Integer; var x, f, d: array of Double; a, b: Double; var info: Integer; incfd: Integer = 1; skip: Integer = 0): Double;

procedure Rpzero2(n: Integer; var a, zr, zi: array of Double; var iter: Integer; var s: array of Double; var info: Integer; iflag: Integer = 0; maxiter: Integer = 0);

type Dfzero_Func = function (x: Double): Double;
procedure Dfzero(f: Dfzero_Func; var b, c: Double; r: Double; var info: Integer; re: Double = 1.0e-10; ae: Double = 1.0e-10);
type Hybrd1_Proc = procedure (n: Integer; var x, fvec: array of Double; var iflag: Integer);
procedure Hybrd1(f: Hybrd1_Proc; n: Integer; var x, fvec: array of Double; var info: Integer; xtol: Double = 1.0e-10);

type Dfmin_Func = function (x: Double): Double;
function Dfmin(a, b: Double; f: Dfmin_Func; tol: Double = 1.0e-10): Double;
type Optif0_Proc = procedure (n: Integer; var x: array of Double; var fval: Double);
procedure Optif0(n: Integer; var x: array of Double; f: Optif0_Proc; var xpls: array of Double; var fpls: Double; var info: Integer);

type Qk15_Func = function (x: Double): Double;
procedure Qk15(f: Qk15_Func; a, b: Double; var result, abserr: Double);
type Qag_Func = function (x: Double): Double;
procedure Qag(f: Qag_Func; a, b: Double; var result, abserr: Double; var neval, last, info: Integer; epsabs: Double = 1.0e-10; epsrel: Double = 1.0e-10; key: Integer = 1; limit: Integer = 100);
type Qagi_Func = function (x: Double): Double;
procedure Qagi(f: Qagi_Func; bound: Double; inf: Integer; var result, abserr: Double; var neval, last, info: Integer; epsabs: Double = 1.0e-10; epsrel: Double = 1.0e-10; limit: Integer = 100);

type Derkf_Proc = procedure (n: Integer; t: Double; var y, yp: array of Double);
procedure Derkf(n: Integer; f: Derkf_Proc; var t: Double; var y: array of Double; tout: Double; var work: array of Double; var iwork: array of Integer; var info: Integer; rtol: Double = 1.0e-10; atol: Double = 1.0e-10; mode: Integer = 0);
procedure DerkfInt(n: Integer; t: Double; var y: array of Double; var work: array of Double);

procedure Rfft1f(n: Integer; var r: array of Double; var wsave: array of Double; var info: Integer; inc: Integer = 1);
procedure Rfft1b(n: Integer; var r: array of Double; var wsave: array of Double; var info: Integer; inc: Integer = 1);
procedure Rfft1i(n: Integer; var wsave: array of Double; var info: Integer);

type Lmdif1_Proc = procedure (m, n: Integer; var x, fvec: array of Double; var iflag: Integer);
procedure Lmdif1(f: Lmdif1_Proc; m, n: Integer; var x, fvec: array of Double; var info: Integer; tol: Double = 1.0e-10);

procedure InitGenRand(s: Cardinal); cdecl; external dll name '_init_genrand';
function GenRandInt32: Cardinal; cdecl; external dll name '_genrand_int32';
function GenRandInt31: Integer; cdecl; external dll name '_genrand_int31';
function GenRandRes53: Double; cdecl; external dll name '_genrand_res53';

function Dlamch(cmach: Char): Double; cdecl; external dll name '_dlamch';

{*******************}
   implementation
{*******************}

{** External functions **}

function _chebs(var c: Double; n: Integer; x: Double; var errno: Integer): Double; cdecl; external dll name 'x_chebs';

function _besi0(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besi0';
function _besi0e(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besi0e';
function _besi1(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besi1';
function _besi1e(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besi1e';

function _besk0(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besk0';
function _besk0e(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besk0e';
function _besk1(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besi1';
function _besk1e(x: Double; var errno: Integer): Double; cdecl; external dll name 'x_besi1e';

function XL_Dlange(norm: Char; m, n, lda: Integer; var a, work: Double): Double; cdecl; external dll name '_dlange';
function XL_Dlansy(norm, uplo: Char; n, lda: Integer; var a, work: Double): Double; cdecl; external dll name '_dlansy';

procedure XL_Dgesv(n, nrhs, lda: Integer; var a: Double; var ipiv: Integer; ldb: Integer; var b: Double; var info: Integer); cdecl; external dll name '_dgesv';
procedure XL_Dgecon(norm: Char; n, lda: Integer; var a: Double; anorm: Double; var rcond: Double; var work: Double; var iwork: Integer; var info: Integer); cdecl; external dll name '_dgecon';

procedure XL_Dposv(uplo: Char; n, nrhs, lda: Integer; var a: Double; ldb: Integer; var b: Double; var info: Integer); cdecl; external dll name '_dposv';
procedure XL_Dpocon(uplo: Char; n, lda: Integer; var a: Double; anorm: Double; var rcond, work: Double; var iwork, info: Integer); cdecl; external dll name '_dpocon';

procedure XL_Dsyev(jobz, uplo: Char; n, lda: Integer; var a, w, work: Double; lwork: Integer; var info: Integer); cdecl; external dll name '_dsyev';

procedure XL_Dgels(trans: Char; m, n, nrhs, lda: Integer; var a: Double; ldb: Integer; var b, work: Double; lwork: Integer; var info: Integer); cdecl; external dll name '_dgels';
procedure XL_Dgecov(job, n, lda: Integer; var a: Double; var ci: Double; var info: Integer); cdecl; external dll name '_dgecov';

procedure XL_Pchse(n: Integer; var x, f, d: Double; incfd: Integer; var work: Double; lwork: Integer; var info: Integer); cdecl; external dll name '_pchse';
procedure XL_Pchfe(n: Integer; var x, f, d: Double; incfd, skip, ne: Integer; var xe, fe: Double; var info: Integer); cdecl; external dll name '_pchfe';
function XL_Pchia(n: Integer; var x, f, d: Double; incfd, skip: Integer; a, b: Double; var info: Integer): Double; cdecl; external dll name '_pchia';

procedure XL_Rpzero2(n: Integer; var a, zr, zi: Double; iflag, maxiter: Integer; var iter: Integer; var s, work: Double; var info: Integer); cdecl; external dll name '_rpzero2';

procedure XL_Dfzero_r(var b, c: Double; r, re, ae: Double; var info: Integer; var xx: Double; yy: Double; var irev: Integer); cdecl; external dll name '_dfzero_r';
procedure XL_Hybrd1_r(n: Integer; var x, fvec: Double; xtol: Double; var work: Double; lwork: Integer; var info: Integer; var xx, yy: Double; var irev: Integer); cdecl; external dll name '_hybrd1_r';

procedure XL_Dfmin_r(a, b, tol: Double; var xx: Double; yy: Double; var irev: Integer); cdecl; external dll name '_dfmin_r';
procedure XL_Optif0_r(n: Integer; var x, xpls, fpls, work: Double; lwork: Integer; var info: Integer; var xx: Double; yy: Double; var irev: Integer); cdecl; external dll name '_optif0_r';

procedure XL_Qk15_r(a, b: Double; var result, abserr, resabs, resasc, xx: Double; yy: Double; var irev: Integer); cdecl; external dll name '_qk15_r';
procedure XL_Qag_r(a, b, epsabs, epsrel: Double; key, limit: Integer; var result, abserr: Double; var neval, last: Integer; var work: Double; lwork: Integer; var iwork, info: Integer; var xx: Double; yy: Double; var irev: Integer); cdecl; external dll name '_qag_r';
procedure XL_Qagi_r(bound: Double; inf: Integer; epsabs, epsrel: Double; limit: Integer; var result, abserr: Double; var neval, last: Integer; var work: Double; lwork: Integer; var iwork, info: Integer; var xx: Double; yy: Double; var irev: Integer); cdecl; external dll name '_qagi_r';

procedure XL_Derkf_r(n: Integer; var t: Double; var y: Double; tout: Double; var rtol, atol: Double; itol, mode: Integer; var work: Double; lwork: Integer; var iwork: Integer; liwork: Integer; var info: Integer; var tt, yy, yyp: Double; var irev: Integer); cdecl; external dll name '_derkf_r';
procedure XL_DerkfInt(n: Integer; t: Double; var y: Double; var work: Double); cdecl; external dll name '_derkf_int';

procedure XL_Rfft1f(n, inc: Integer; var r: Double; lr: Integer; var wsave: Double; lwsave: Integer; var work: Double; lwork: Integer; var info: Integer); cdecl; external dll name '_rfft1f';
procedure XL_Rfft1b(n, inc: Integer; var r: Double; lr: Integer; var wsave: Double; lwsave: Integer; var work: Double; lwork: Integer; var info: Integer); cdecl; external dll name '_rfft1b';
procedure XL_Rfft1i(n: Integer; var wsave: Double; lwsave: Integer; var info: Integer); cdecl; external dll name '_rfft1i';

procedure XL_Lmdif1_r(m, n: Integer; var x, fvec: Double; tol: Double; var work: Double; lwork: Integer; var iwork: Integer; var info: Integer; var xx, yy: Double; var irev: Integer); cdecl; external dll name '_lmdif1_r';

{** Function definitions **}

function Chebs(var c: array of Double; n: Integer; x: Double; var errno: Integer): Double;
begin
	Chebs := _chebs(c[0], n, x, errno);
end;

function Besi0(x: Double; var errno: Integer; kode: Integer): Double;
begin
	if kode <> 2 then
		Besi0 := _besi0(x, errno)
	else
		Besi0 := _besi0e(x, errno);
end;

function Besi1(x: Double; var errno: Integer; kode: Integer): Double;
begin
	if kode <> 2 then
		Besi1 := _besi1(x, errno)
	else
		Besi1 := _besi1e(x, errno);
end;

function Besk0(x: Double; var errno: Integer; kode: Integer): Double;
begin
	if kode <> 2 then
		Besk0 := _besk0(x, errno)
	else
		Besk0 := _besk0e(x, errno);
end;

function Besk1(x: Double; var errno: Integer; kode: Integer): Double;
begin
	if kode <> 2 then
		Besk1 := _besk1(x, errno)
	else
		Besk1 := _besk1e(x, errno);
end;

function Dlange(norm: Char; m, n, lda: Integer; var a: array of Double): Double;
	var work: array of Double;
begin
	SetLength(work, m);
	Dlange := XL_Dlange(norm, m, n, lda, a[0], work[0]);
end;

function Dlansy(norm, uplo: Char; n, lda: Integer; var a: array of Double): Double;
	var work: array of Double;
begin
	SetLength(work, n);
	Dlansy := XL_Dlansy(norm, uplo, n, lda, a[0], work[0]);
end;

procedure Dgesv(n, lda: Integer; var a: array of Double; var ipiv: array of Integer; ldb: Integer; var b: array of Double; var info: Integer; nrhs: Integer);
begin
	XL_Dgesv(n, nrhs, lda, a[0], ipiv[0], ldb, b[0], info);
end;

procedure Dgecon(norm: Char; n, lda: Integer; var a: array of Double; anorm: Double; var rcond: Double; var info: Integer);
	var work: array of Double; iwork: array of Integer;
begin
	SetLength(work, 4*n);
	SetLength(iwork, n);
	XL_Dgecon(norm, n, lda, a[0], anorm, rcond, work[0], iwork[0], info);
end;

procedure Dposv(uplo: Char; n, lda: Integer; var a: array of Double; ldb: Integer; var b: array of Double; var info: Integer; nrhs: Integer);
begin
	XL_Dposv(uplo, n, nrhs, lda, a[0], ldb, b[0], info);
end;

procedure Dpocon(uplo: Char; n, lda: Integer; var a: array of Double; anorm: Double; var rcond: Double; var info: Integer);
	var work: array of Double; iwork: array of Integer;
begin
	SetLength(work, 3*n);
	SetLength(iwork, n);
	XL_Dpocon(uplo, n, lda, a[0], anorm, rcond, work[0], iwork[0], info);
end;

procedure Dsyev(jobz, uplo: Char; n, lda: Integer; var a: array of Double; var w: array of Double; var info: Integer);
	var work: array of Double; lwork: Integer;
begin
	SetLength(work, 1);
	lwork := -1;
	XL_Dsyev(jobz, uplo, n, lda, a[0], w[0], work[0], lwork, info);
	if info = 0 then
	begin
		lwork := Trunc(work[0]);
		SetLength(work, lwork);
		XL_Dsyev(jobz, uplo, n, lda, a[0], w[0], work[0], lwork, info);
	end;
end;

procedure Dgels(trans: Char; m, n, lda: Integer; var a: array of Double; ldb: Integer; var b: array of Double; var info: Integer; nrhs: Integer);
	var work: array of Double; lwork: Integer;
begin
	SetLength(work, 1);
	lwork := -1;
	XL_Dgels(trans, m, n, nrhs, lda, a[0], ldb, b[0], work[0], lwork, info);
	if info = 0 then
	begin
		lwork := Trunc(work[0]);
		SetLength(work, lwork);
		XL_Dgels(trans, m, n, nrhs, lda, a[0], ldb, b[0], work[0], lwork, info);
	end;
end;

procedure Dgecov(job, n, lda: Integer; var a: array of Double; var ci: array of Double; var info: Integer);
begin
	XL_Dgecov(job, n, lda, a[0], ci[0], info);
end;

procedure Pchse(n: Integer; var x, f, d: array of Double; var info: Integer; incfd: Integer);
	var work: array of Double; lwork: Integer;
begin
	if n < 2 then
	begin
		info := -1;
		Exit;
	end;
	lwork := 2*n;
	SetLength(work, lwork);
	XL_Pchse(n, x[0], f[0], d[0], incfd, work[0], lwork, info);
end;

procedure Pchfe(n: Integer; var x, f, d: array of Double; ne: Integer; var xe, fe: array of Double; var info: Integer; incfd, skip: Integer);
begin
	XL_Pchfe(n, x[0], f[0], d[0], incfd, skip, ne, xe[0], fe[0], info);
end;

function Pchia(n: Integer; var x, f, d: array of Double; a, b: Double; var info: Integer; incfd, skip: Integer): Double;
begin
	Pchia := XL_Pchia(n, x[0], f[0], d[0], incfd, skip, a, b, info);
end;

procedure Rpzero2(n: Integer; var a, zr, zi: array of Double; var iter: Integer; var s: array of Double; var info: Integer; iflag: Integer; maxiter: Integer);
	var work: array of Double;
begin
	if n < 1 then
	begin
		info := -1;
		Exit;
	end;
	if maxiter <= 0 then
		maxiter := 25*n;
	SetLength(work, 8*n + 6);
	XL_Rpzero2(n, a[0], zr[0], zi[0], iflag, maxiter, iter, s[0], work[0], info);
end;

procedure Dfzero(f: Dfzero_Func; var b, c: Double; r: Double; var info: Integer; re: Double; ae: Double);
	var irev: Integer; xx, yy: Double;
begin
	irev := 0; yy := 0.0;
	repeat
		XL_Dfzero_r(b, c, r, re, ae, info, xx, yy, irev);
		if irev <> 0 then yy := f(xx);
	until irev = 0;
end;

procedure Hybrd1(f: Hybrd1_Proc; n: Integer; var x, fvec: array of Double; var info: Integer; xtol: Double);
	var iflag, irev: Integer; xx, yy: array of Double;
	var work: array of Double; lwork: Integer;
begin
	if n < 1 then
	begin
		info := -2;
		Exit;
	end;
	if xtol < 0 then
	begin
		info := -5;
		Exit;
	end;
	SetLength(xx, n);
	SetLength(yy, n);
	lwork := n*(3*n + 13) div 2;
	SetLength(work, lwork);
	irev := 0;
	repeat
		XL_Hybrd1_r(n, x[0], fvec[0], xtol, work[0], lwork, info, xx[0], yy[0], irev);
		if (irev >= 1) and (irev <= 4) then
		begin
			iflag := (irev + 1) div 2;
			f(n, xx, yy, iflag);
		end;
	until irev = 0;
end;

function Dfmin(a, b: Double; f: Dfmin_Func; tol: Double): Double;
	var irev: Integer; xx, yy: Double;
begin
	irev := 0; yy := 0.0;
	repeat
		XL_Dfmin_r(a, b, tol, xx, yy, irev);
		if irev <> 0 then yy := f(xx);
	until irev = 0;
	Dfmin := xx;
end;

procedure Optif0(n: Integer; var x: array of Double; f: Optif0_Proc; var xpls: array of Double; var fpls: Double; var info: Integer);
	var irev: Integer; xx: array of Double; yy: Double;
	var work: array of Double; lwork: Integer;
begin
	if n < 2 then
	begin
		info := -1;
		Exit;
	end;
	SetLength(xx, n);
	lwork := n*(n + 10);
	SetLength(work, lwork);
	irev := 0;
	repeat
		XL_Optif0_r(n, x[0], xpls[0], fpls, work[0], lwork, info, xx[0], yy, irev);
		if (irev >= 1) and (irev <= 20) then
		begin
			f(n, xx, yy);
		end;
	until irev = 0;
end;

procedure Qk15(f: Qk15_Func; a, b: Double; var result, abserr: Double);
	var irev: Integer; resabs, resasc, xx, yy: Double;
begin
	irev := 0; yy := 0.0;
	repeat
		XL_Qk15_r(a, b, result, abserr, resabs, resasc, xx, yy, irev);
		if irev <> 0 then yy := f(xx);
	until irev = 0;
end;

procedure Qag(f: Qag_Func; a, b: Double; var result, abserr: Double; var neval, last, info: Integer; epsabs: Double; epsrel: Double; key: Integer; limit: Integer);
	var irev: Integer; xx, yy: Double;
	var work: array of Double; lwork: Integer; iwork: array of Integer;
begin
	if limit < 1 then
	begin
		info := -7;
		Exit;
	end;
	lwork := 4*limit;
	SetLength(work, lwork);
	SetLength(iwork, limit);
	irev := 0; yy := 0.0;
	repeat
		XL_Qag_r(a, b, epsabs, epsrel, key, limit, result, abserr, neval, last, work[0], lwork, iwork[0], info, xx, yy, irev);
		if irev <> 0 then yy := f(xx);
	until irev = 0;
end;

procedure Qagi(f: Qagi_Func; bound: Double; inf: Integer; var result, abserr: Double; var neval, last, info: Integer; epsabs: Double; epsrel: Double; limit: Integer);
	var irev: Integer; xx, yy: Double;
	var work: array of Double; lwork: Integer; iwork: array of Integer;
begin
	if limit < 1 then
	begin
		info := -6;
		Exit;
	end;
	lwork := 4*limit;
	SetLength(work, lwork);
	SetLength(iwork, limit);
	irev := 0; yy := 0.0;
	repeat
		XL_Qagi_r(bound, inf, epsabs, epsrel, limit, result, abserr, neval, last, work[0], lwork, iwork[0], info, xx, yy, irev);
		if irev <> 0 then yy := f(xx);
	until irev = 0;
end;

procedure Derkf(n: Integer; f: Derkf_Proc; var t: Double; var y: array of Double; tout: Double; var work: array of Double; var iwork: array of Integer; var info: Integer; rtol: Double; atol: Double; mode: Integer);
var
	itol: Integer; lwork: Integer; liwork: Integer;
	irev: Integer; tt: Double; yy, yyp: array of Double;
begin
	if n < 1 then
	begin
		info := -1;
		Exit;
	end;
	itol := 0;
	lwork := Length(work);
	liwork := Length(iwork);
	SetLength(yy, n);
	SetLength(yyp, n);
	irev := 0;
	repeat
		XL_Derkf_r(n, t, y[0], tout, rtol, atol, itol, mode, work[0], lwork, iwork[0], liwork, info, tt, yy[0], yyp[0], irev);
		if (irev >= 1) and (irev <= 11) then
		begin
			f(n, tt, yy, yyp);
		end;
	until irev = 0;
	if (info = -2) or (info = -4) then
		info := info - 1
	else if (info = -5) or (info = -6) then
		info := info - 4
	else if info = -10 then
		info := -6
	else if info < 0 then
		info := info + 5;
end;

procedure DerkfInt(n: Integer; t: Double; var y: array of Double; var work: array of Double);
begin
	if (n < 1) or (Length(y) < n) or (Length(work) < 9*n + 20) then
		Exit;
	XL_DerkfInt(n, t, y[0], work[0]);
end;

procedure Rfft1f(n: Integer; var r: array of Double; var wsave: array of Double; var info: Integer; inc: Integer);
var
	lr: Integer; lwsave: Integer;
	work: array of Double; lwork: Integer;
begin
	if n < 1 then
	begin
		info := -1;
		Exit;
	end;
	lr := Length(r); lwsave := Length(wsave);
	lwork := n;
	SetLength(work, lwork);
	XL_Rfft1f(n, inc, r[0], lr, wsave[0], lwsave, work[0], lwork, info);
end;

procedure Rfft1b(n: Integer; var r: array of Double; var wsave: array of Double; var info: Integer; inc: Integer);
var
	lr: Integer; lwsave: Integer;
	work: array of Double; lwork: Integer;
begin
	if n < 1 then
	begin
		info := -1;
		Exit;
	end;
	lr := Length(r); lwsave := Length(wsave);
	lwork := n;
	SetLength(work, lwork);
	XL_Rfft1b(n, inc, r[0], lr, wsave[0], lwsave, work[0], lwork, info);
end;

procedure Rfft1i(n: Integer; var wsave: array of Double; var info: Integer);
var
	lwsave: Integer;
begin
	lwsave := Length(wsave);
	XL_Rfft1i(n, wsave[0], lwsave, info);
end;

procedure Lmdif1(f: Lmdif1_Proc; m, n: Integer; var x, fvec: array of Double; var info: Integer; tol: Double);
	var iflag, irev: Integer; xx, yy: array of Double;
	var work: array of Double; lwork: Integer; var iwork: array of Integer; liwork: Integer;
begin
	if m < n then
	begin
		info := -2;
		Exit;
	end;
	if n < 1 then
	begin
		info := -3;
		Exit;
	end;
	if tol < 0 then
	begin
		info := -6;
		Exit;
	end;
	SetLength(xx, n);
	SetLength(yy, m);
	lwork := n*(m + 6) + 2*m;
	SetLength(work, lwork);
	liwork := n;
	SetLength(iwork, liwork);
	irev := 0;
	repeat
		XL_Lmdif1_r(m, n, x[0], fvec[0], tol, work[0], lwork, iwork[0], info, xx[0], yy[0], irev);
		if (irev >= 1) and (irev <= 3) then
		begin
			iflag := 1;
			if irev = 3 then iflag := 2;
			f(m, n, xx, yy, iflag);
			if iflag < 0 then
			begin
				info := 5;
				irev := 0;
			end;
		end;
	until irev = 0;
end;

end.
