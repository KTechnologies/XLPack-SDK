# *********************************************
# *                                           *
# *  XLPack Julia Numerical Library           *
# *  Version 6.0 (February 12, 2022)          *
# *  (C) 2014-2022  K Technologies            *
# *                                           *
# *********************************************/
module XLPack
export errno, d1num, chebs, sqrt1pm1, powm1, cos_pi, sin_pi
export li, ei, e1, digamma, besj0, besj1, besjnu, besy0, besy1
export besynu, besi0, besi1, besinu, besk0, besk1, besknu, celli1, celli2, celli3
export dconst, dlange, dlansy, dgesv, dgecon, dposv, dpocon, dsyev, dgels, dgecov
export pchse, pchfe, pchia, rpzero2, dfzero, hybrd1, dfmin, optif0, qk15, qag
export qagi, derkf, derkf_int, rfft1f, rfft1b, rfft1i, lmdif1, init_genrand
export genrand_int32, genrand_int31, genrand_res53, dlamch, get_errno, set_errno
@static if Sys.isapple()
    const lib = "/Library/Application Support/Microsoft/Office365/User Content.localized/Add-Ins.localized/XLPackLite.dylib"
elseif Sys.WORD_SIZE == 64
    const lib = "XLPack"
else # Sys.WORD_SIZE == 32
    const lib = "XLPack_32"
end
errno = Int32(0)
function d1num(i::Integer)
	ccall((:d1num, lib), Cdouble, (Cint, ), i)
end
function chebs(c::Array{Float64}, n::Integer, x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_chebs, lib), Cdouble, (Ptr{Cdouble}, Cint, Cdouble, Ptr{Cint}, ), c, n, x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function sqrt1pm1(x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_sqrt1pm1, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function powm1(x::Real, y::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_powm1, lib), Cdouble, (Cdouble, Cdouble, Ptr{Cint}, ), x, y, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function cos_pi(x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_cos_pi, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function sin_pi(x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_sin_pi, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function li(x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_li, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function ei(x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_ei, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function e1(x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_e1, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function digamma(x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_digamma, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function besj0(x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_besj0, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function besj1(x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_besj1, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	if errno1[] != 0
		errno = errno1[]
	end
	result
end
function besjnu(nu::Real, x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_besjnu, lib), Cdouble, (Cdouble, Cdouble, Ptr{Cint}, ), nu, x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function besy0(x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_besy0, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function besy1(x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_besy1, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function besynu(nu::Real, x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_besynu, lib), Cdouble, (Cdouble, Cdouble, Ptr{Cint}, ), nu, x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function besi0(x::Real, scaling::Bool = false)
	errno1 = Ref{Cint}()
	if ! scaling
		result = ccall((:_besi0, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	else
		result = ccall((:_besi0e, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	end
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function besi1(x::Real, scaling::Bool = false)
	errno1 = Ref{Cint}()
	if ! scaling
		result = ccall((:_besi1, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	else
		result = ccall((:_besi1e, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	end
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function besinu(nu::Real, x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_besinu, lib), Cdouble, (Cdouble, Cdouble, Ptr{Cint}, ), nu, x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function besk0(x::Real, scaling::Bool = false)
	errno1 = Ref{Cint}()
	if ! scaling
		result = ccall((:_besk0, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	else
		result = ccall((:_besk0e, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	end
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function besk1(x::Real, scaling::Bool = false)
	errno1 = Ref{Cint}()
	if ! scaling
		result = ccall((:_besk1, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	else
		result = ccall((:_besk1e, lib), Cdouble, (Cdouble, Ptr{Cint}, ), x, errno1)
	end
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function besknu(nu::Real, x::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_besknu, lib), Cdouble, (Cdouble, Cdouble, Ptr{Cint}, ), nu, x, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function celli1(k::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_celli1, lib), Cdouble, (Cdouble, Ptr{Cint}, ), k, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function celli2(k::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_celli2, lib), Cdouble, (Cdouble, Ptr{Cint}, ), k, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function celli3(n::Real, k::Real)
	errno1 = Ref{Cint}()
	result = ccall((:_celli3, lib), Cdouble, (Cdouble, Cdouble, Ptr{Cint}, ), n, k, errno1)
	if errno1[] != 0
		global errno = errno1[]
	end
	result
end
function dconst(i::Integer)
	ccall((:dconst, lib), Cdouble, (Cint, ), i)
end
function dlange(norm::Char, m::Integer, n::Integer, a::Array{Float64})
    anorm = 0.0
    info = 0
    if ! occursin(norm, "Mm1OoIiFfEe")
        info = -1
    elseif m < 0
        info = -2
    elseif n < 0
        info = -3
    elseif ndims(a) == 1
        lda = m
    elseif ndims(a) == 2
        lda = size(a, 1)
    else
        info = -4
    end
    if info == 0
        work = Vector{Cdouble}(undef, m)
        anorm = ccall((:dlange, lib), Cdouble, (Cchar, Cint, Cint, Cint, Ref{Cdouble}, Ptr{Cdouble}), norm, m, n,  lda, a, work)
    end
    anorm, info
end
function dlansy(norm::Char, uplo::Char, n::Integer, a::Array{Float64})
    anorm = 0.0
    info = 0
    if ! occursin(norm, "Mm1OoIiFfEe")
        info = -1
    elseif ! occursin(uplo, "UuLl")
        info = -2
    elseif n < 0
        info = -3
    elseif ndims(a) == 1
        lda = n
    elseif ndims(a) == 2
        lda = size(a, 1)
    else
        info = -4
    end
    if info == 0
        work = Vector{Cdouble}(undef, n)
        anorm = ccall((:dlansy, lib), Cdouble, (Cchar, Cchar, Cint, Cint, Ref{Cdouble}, Ptr{Cdouble}), norm, uplo, n,  lda, a, work)
    end
    anorm, info
end
function dgesv(n::Integer, a::Array{Float64}, ipiv::Array{Int32}, b::Array{Float64}, nrhs::Integer = 1)
    info = Ref{Cint}()
    if ndims(a) == 1
        lda = n
    else
        lda = size(a, 1)
    end
    if ndims(b) == 1
        ldb = n
    else
        ldb = size(b, 1)
    end
    ccall((:dgesv, lib), Nothing, (Cint, Cint, Cint, Ptr{Cdouble}, Ptr{Cint}, Cint, Ptr{Cdouble}, Ptr{Cint}), n, nrhs, lda, a, ipiv, ldb, b, info)
    if info[] == -2
        info[] = -5
    elseif info[] == -3
        info[] = -2
    elseif info[] == -6
        info[] = -4
    end
    info[]
end
function dgecon(norm::Char, n::Integer, a::Array{Float64}, anorm::Real)
    info = Ref{Cint}()
    rcond = Ref{Cdouble}()
    if ndims(a) == 1
        lda = n
        if size(a, 1) < lda*n
            return 0.0, -3
        end
    elseif ndims(a) == 2
        lda = size(a, 1)
        if size(a, 2) < n
            return 0.0, -3
        end
    else
        return 0.0, -3
    end
    work = Vector{Cdouble}(undef, 4*n)
    iwork = Vector{Cint}(undef, n)
    ccall((:dgecon, lib), Nothing, (Cchar, Cint, Cint, Ref{Cdouble}, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), norm, n, lda, a, anorm, rcond, work, iwork, info)
    if info[] == -5
        info[] = -4
    end
    rcond[], info[]
end
function dposv(uplo::Char, n::Integer, a::Array{Float64}, b::Array{Float64}, nrhs::Integer = 1)
    info = Ref{Cint}()
    if ndims(a) == 1
        lda = n
    else
        lda = size(a, 1)
    end
    if ndims(b) == 1
        ldb = n
    else
        ldb = size(b, 1)
    end
    ccall((:dposv, lib), Nothing, (Cchar, Cint, Cint, Cint, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Ptr{Cint}), uplo, n, nrhs, lda, a, ldb, b, info)
    if info[] == -6
        info[] = -5
    end
    info[]
end
function dpocon(uplo::Char, n::Integer, a::Array{Float64}, anorm::Real)
    info = Ref{Cint}()
    rcond = Ref{Cdouble}()
    if ndims(a) == 1
        lda = n
        if size(a, 1) < lda*n
            return 0.0, -3
        end
    elseif ndims(a) == 2
        lda = size(a, 1)
        if size(a, 2) < n
            return 0.0, -3
        end
    else
        return 0.0, -3
    end
    work = Vector{Cdouble}(undef, 3*n)
    iwork = Vector{Cint}(undef, n)
    ccall((:dpocon, lib), Nothing, (Cchar, Cint, Cint, Ref{Cdouble}, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}), uplo, n, lda, a, anorm, rcond, work, iwork, info)
    if info[] == -5
        info[] = -4
    end
    rcond[], info[]
end
function dsyev(jobz::Char, uplo::Char, n::Integer, a::Array{Float64}, w::Array{Float64})
    info = Ref{Cint}()
    if ndims(a) == 1
        lda = n
    else
        lda = size(a, 1)
    end
    work0 = Ref{Cdouble}()
    lwork::Cint = -1
    ccall((:dsyev, lib), Nothing, (Cchar, Cchar, Cint, Cint, Ref{Cdouble}, Ref{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}), jobz, uplo, n, lda, a, w, work0, lwork, info)
    if info[] == 0
        lwork = work0[]
        work = Vector{Cdouble}(undef, lwork)
        ccall((:dsyev, lib), Nothing, (Cchar, Cchar, Cint, Cint, Ref{Cdouble}, Ref{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}), jobz, uplo, n, lda, a, w, work, lwork, info)
    end
    info[]
end
function dgels(trans::Char, m::Integer, n::Integer, a::Array{Float64}, b::Array{Float64}, nrhs::Integer = 1)
    info = Ref{Cint}()
    if ndims(a) == 1
        lda = n
    else
        lda = size(a, 1)
    end
    if ndims(b) == 1
        ldb = m
    else
        ldb = size(b, 1)
    end
    work0 = Ref{Cdouble}()
    lwork::Cint = -1
    ccall((:dgels, lib), Nothing, (Cchar, Cint, Cint, Cint, Cint, Ref{Cdouble}, Cint, Ref{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}), trans, m, n, nrhs, lda, a, ldb, b, work0, lwork, info)
    if info[] == 0
        lwork = work0[]
        work = Vector{Cdouble}(undef, lwork)
        ccall((:dgels, lib), Nothing, (Cchar, Cint, Cint, Cint, Cint, Ref{Cdouble}, Cint, Ref{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}), trans, m, n, nrhs, lda, a, ldb, b, work, lwork, info)
    end
    if info[] == -7
        info[] = -6
    end
    info[]
end
function dgecov(job::Integer, n::Integer, a::Array{Float64}, ci::Array{Float64})
    info = Ref{Cint}()
    if ndims(a) == 1
        lda = n
    else
        lda = size(a, 1)
    end
    ccall((:dgecov, lib), Nothing, (Cint, Cint, Cint, Ref{Cdouble}, Ref{Cdouble}, Ptr{Cint}), job, n, lda, a, ci, info)
    info[]
end
function pchse(n::Integer, x::Array{Float64}, f::Array{Float64}, d::Array{Float64}, incfd::Integer = 1)
    info = Ref{Cint}()
    lwork::Cint = 2*n
    work = Vector{Cdouble}(undef, lwork)
    ccall((:pchse, lib), Nothing, (Cint, Ref{Cdouble}, Ref{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cint}), n, x, f, d, incfd, work, lwork, info)
    info[]
end
function pchfe(n::Integer, x::Array{Float64}, f::Array{Float64}, d::Array{Float64}, ne::Integer, xe::Array{Float64}, fe::Array{Float64}, incfd::Integer = 1; skip::Bool = false)
    info = Ref{Cint}()
    iskip::Cint = skip ? 1 : 0
    ccall((:pchfe, lib), Nothing, (Cint, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Cint, Cint, Cint, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cint}), n, x, f, d, incfd, iskip, ne, xe, fe, info)
    info[]
end
function pchia(n::Integer, x::Array{Float64}, f::Array{Float64}, d::Array{Float64}, a::Real, b::Real, incfd::Integer = 1; skip::Bool = false)
    info = Ref{Cint}()
    iskip::Cint = skip ? 1 : 0
    s = ccall((:pchia, lib), Cdouble, (Cint, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Cint, Cint, Cdouble, Cdouble, Ptr{Cint}), n, x, f, d, incfd, iskip, a, b, info)
    s, info[]
end
function rpzero2(n::Integer, a::Array{Float64}, rr::Array{Float64}, ri::Array{Float64}, s::Array{Float64}, iflag::Integer = 0; maxiter::Integer = 100)
    info = Ref{Cint}()
    iter = Ref{Cint}()
    lwork::Cint = 8*n + 6
    work = Vector{Cdouble}(undef, lwork)
    s = ccall((:rpzero2, lib), Nothing, (Cint, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}), n, a, rr, ri, iflag, maxiter, iter, s, work, info)
    iter[], info[]
end
function dfzero(f::Function, b::Real, c::Real, r::Real, re::Real = 1.0e-10, ae::Real = 1.0e-10)
    info = Ref{Cint}()
    irev = Ref{Cint}()
    b1 = Ref{Cdouble}()
    c1 = Ref{Cdouble}()
    b1[] = b
    c1[] = c
    xx = Ref{Cdouble}()
    xx[] = 0
    yy::Cdouble = 0
    irev[] = 0
    while true
        ccall((:dfzero_r, lib), Nothing, (Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Ptr{Cint}), b1, c1, r, re, ae, info, xx, yy, irev)
        if irev[] == 0
            break
        end
        yy = f(xx[])
    end
    b1[], info[]
end
function hybrd1(f::Function, n::Integer, x::Array{Float64}, fvec::Array{Float64}, xtol::Real = 1.0e-10)
    info = Ref{Cint}()
    irev = Ref{Cint}()
    if n < 1
        return -2
    end
    if xtol < 0
        return -5
    end
    lwork::Cint = n*(3*n + 13)/2
    work = Vector{Cdouble}(undef, lwork)
    xx = Vector{Cdouble}(undef, n)
    yy = Vector{Cdouble}(undef, n)
    irev[] = 0
    while true
        ccall((:hybrd1_r, lib), Nothing, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ref{Cdouble}, Ptr{Cint}), n, x, fvec, xtol, work, lwork, info, xx, yy, irev)
        if irev[] == 1 || irev[] == 2
            iflag = 1
        elseif irev[] == 3 || irev[] == 4
            iflag = 2
        elseif irev[] == 0
            break
        end
        f(n, xx, yy, iflag)
    end
    info[]
end
function dfmin(a::Real, b::Real, f::Function, tol::Real = 1.0e-10)
    irev = Ref{Cint}()
    xx = Ref{Cdouble}()
    yy::Cdouble = 0
    irev[] = 0
    while true
        ccall((:dfmin_r, lib), Nothing, (Cdouble, Cdouble, Cdouble, Ptr{Cdouble}, Cdouble, Ptr{Cint}), a, b, tol, xx, yy, irev)
        if irev[] == 0
            break
        end
        yy = f(xx[])
    end
    xx[]
end
function optif0(n::Integer, x::Array{Float64}, f::Function, xpls::Array{Float64})
    info = Ref{Cint}()
    irev = Ref{Cint}()
    if n < 1
        return 0.0, -1
    end
    fpls = Ref{Cdouble}()
    lwork = n*(n + 10)
    work = Vector{Cdouble}(undef, lwork)
    xx = Vector{Cdouble}(undef, n)
    yy::Cdouble = 0
    irev[] = 0
    while true
        ccall((:optif0_r, lib), Nothing, (Cint, Ref{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Ptr{Cint}), n, x, xpls, fpls, work, lwork, info, xx, yy, irev)
        if irev[] >= 1 && irev[] <= 20
            yy = f(n, xx)
        elseif irev[] == 0
            break
        end
    end
    fpls[], info[]
end
function qk15(f::Function, a::Real, b::Real)
    info = Ref{Cint}()
    irev = Ref{Cint}()
    result = Ref{Cdouble}()
    abserr = Ref{Cdouble}()
    resabs = Ref{Cdouble}()
    resasc = Ref{Cdouble}()
    xx = Ref{Cdouble}()
    yy::Cdouble = 0
    irev[] = 0
    while true
        ccall((:qk15_r, lib), Nothing, (Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cint}), a, b, result, abserr, resabs, resasc, xx, yy, irev)
        if irev[] >= 1 && irev[] <= 5
            yy = f(xx[])
        elseif irev[] == 0
            break
        end
    end
    result[], abserr[]
end
function qag(f::Function, a::Real, b::Real, epsabs::Real = 1.0e-10, epsrel::Real = 1.0e-10, key::Integer = 1, limit::Integer = 100)
    info = Ref{Cint}()
    irev = Ref{Cint}()
    result = Ref{Cdouble}()
    abserr = Ref{Cdouble}()
    neval = Ref{Cint}()
    last = Ref{Cint}()
    lwork = 4*limit
    work = Vector{Cdouble}(undef, lwork)
    iwork = Vector{Cint}(undef, limit)
    xx = Ref{Cdouble}()
    yy::Cdouble = 0
    irev[] = 0
    while true
        ccall((:qag_r, lib), Nothing, (Cdouble, Cdouble, Cdouble, Cdouble, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Ptr{Cint}), a, b, epsabs, epsrel, key, limit, result, abserr, neval, last, work, lwork, iwork, info, xx, yy, irev)
        if irev[] >= 1 && irev[] <= 15
            yy = f(xx[])
        elseif irev[] == 0
            break
        end
    end
    if info[] == -6
        info[] = -7
    end
    result[], abserr[], info[]
end
function qagi(f::Function, bound::Real, inf::Integer, epsabs::Real = 1.0e-10, epsrel::Real = 1.0e-10, limit::Integer = 100)
    info = Ref{Cint}()
    irev = Ref{Cint}()
    result = Ref{Cdouble}()
    abserr = Ref{Cdouble}()
    neval = Ref{Cint}()
    last = Ref{Cint}()
    lwork = 4*limit
    work = Vector{Cdouble}(undef, lwork)
    iwork = Vector{Cint}(undef, limit)
    xx = Ref{Cdouble}()
    yy::Cdouble = 0
    irev[] = 0
    while true
        ccall((:qagi_r, lib), Nothing, (Cdouble, Cint, Cdouble, Cdouble, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Ptr{Cint}), bound, inf, epsabs, epsrel, limit, result, abserr, neval, last, work, lwork, iwork, info, xx, yy, irev)
        if irev[] >= 1 && irev[] <= 18
            yy = f(xx[])
        elseif irev[] == 0
            break
        end
    end
    if info[] == -5
        info[] = -6
    end
    result[], abserr[], info[]
end
function derkf(info::Integer, n::Integer, f::Function, t::Real, y::Array{Float64}, tout::Real, wsave::Array{Float64}, iwsave::Array{Int32}; rtol::Real = 1.0e-10, atol::Real = 1.0e-10, mode::Integer = 0)
    info1 = Ref{Cint}()
    irev = Ref{Cint}()
    t1 = Ref{Cdouble}()
    rtol1 = Ref{Cdouble}()
    atol1 = Ref{Cdouble}()
    if n < 1
        return t, -2
    elseif ndims(y) != 1 || size(y, 1) < n
        return t, -5
    elseif ndims(wsave) != 1
        return t, -7
    elseif ndims(iwsave) != 1
        return t, -8
    end
    t1[] = t
    rtol1[] = rtol
    atol1[] = atol
    itol::Cint = 0
    lwsave = size(wsave, 1)
    liwsave = size(iwsave, 1)
    info1[] = info
    tt = Ref{Cdouble}()
    yy = Vector{Cdouble}(undef, n)
    yyp = Vector{Cdouble}(undef, n)
    irev[] = 0
    while true
        ccall((:derkf_r, lib), Nothing, (Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Cint, Ptr{Cint}, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ref{Cdouble}, Ptr{Cint}), n, t1, y, tout, rtol1, atol1, itol, mode, wsave, lwsave, iwsave, liwsave, info1, tt, yy, yyp, irev)
        if irev[] >= 1 && irev[] <= 11
            f(n, tt[], yy, yyp)
        elseif irev[] == 0
            break
        end
    end
    if info1[] == -12
        info1[] = -1
    elseif info1[] == -10
        info1[] = -7
    elseif info1[] == -12
        info1[] = -8
    elseif info1[] == -5
        info1[] = -9
    elseif info1[] == -6
        info1[] = -10
    elseif info1[] < 0
        info1[] = info1[] - 2
    end
    t1[], info1[]
end
function derkf_int(n::Integer, t::Real, y::Array{Float64}, wsave::Array{Float64})
    if n < 1
        return -1
    elseif ndims(y) != 1
        return -3
    elseif ndims(wsave) != 1
        return -4
    end
    ccall((:derkf_int, lib), Nothing, (Cint, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), n, t, y, wsave)
    return 0
end
function rfft1f(n::Integer, r::Array{Float64}, wsave::Array{Float64}, inc::Integer = 1)
    info = Ref{Cint}()
    if ndims(r) != 1
        return -2
    elseif ndims(wsave) != 1
        return -3
    end
    lr = size(r, 1)
    lwsave = size(wsave, 1)
    lwork = n
    work = Vector{Cdouble}(undef, lwork)
    ccall((:rfft1f, lib), Nothing, (Cint, Cint, Ptr{Cdouble}, Cint, Ref{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cint}), n, inc, r, lr, wsave, lwsave, work, lwork, info)
    if info[] == -2
        info[] = -4
    elseif info[] == -4
        info[] = -2
    elseif info[] == -6
        info[] = -3
    end
    info[]
end
function rfft1b(n::Integer, r::Array{Float64}, wsave::Array{Float64}, inc::Integer = 1)
    info = Ref{Cint}()
    if ndims(r) != 1
        return -2
    elseif ndims(wsave) != 1
        return -3
    end
    lr = size(r, 1)
    lwsave = size(wsave, 1)
    lwork = n
    work = Vector{Cdouble}(undef, lwork)
    ccall((:rfft1b, lib), Nothing, (Cint, Cint, Ptr{Cdouble}, Cint, Ref{Cdouble}, Cint, Ptr{Cdouble}, Cint, Ptr{Cint}), n, inc, r, lr, wsave, lwsave, work, lwork, info)
    if info[] == -2
        info[] = -4
    elseif info[] == -4
        info[] = -2
    elseif info[] == -6
        info[] = -3
    end
    info[]
end
function rfft1i(n::Integer, wsave::Array{Float64})
    info = Ref{Cint}()
    if ndims(wsave) != 1
        return -3
    end
    lwsave = size(wsave, 1)
    ccall((:rfft1i, lib), Nothing, (Cint, Ptr{Cdouble}, Cint, Ptr{Cint}), n, wsave, lwsave, info)
    if info[] == -3
        info[] = -2
    end
    info[]
end
function lmdif1(f::Function, m::Integer, n::Integer, x::Array{Float64}, fvec::Array{Float64}, tol::Real = 1.0e-10)
    info = Ref{Cint}()
    irev = Ref{Cint}()
    if m < n
        return -2
    elseif n < 1
        return -3
    elseif tol < 0
        return -6
    end
    lwork = n*(m + 5) + m
    work = Vector{Cdouble}(undef, lwork)
    liwork = n
    iwork = Vector{Cint}(undef, lwork)
    xx = Vector{Cdouble}(undef, n)
    yy = Vector{Cdouble}(undef, m)
    irev[] = 0
    while true
        ccall((:lmdif1_r, lib), Nothing, (Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ref{Cdouble}, Ptr{Cint}), m, n, x, fvec, tol, work, lwork, iwork, info, xx, yy, irev)
        if irev[] == 1 || irev[] == 2
            iflag = 1
        elseif irev[] == 3
            iflag = 2
        elseif irev[] == 0
            break
        end
        f(m, n, xx, yy, iflag)
    end
    info[]
end
function init_genrand(s::Integer)
	ccall((:init_genrand, lib), Nothing, (Culong, ), s)
end
function genrand_int32()
	ccall((:genrand_int32, lib), Cint, ())
end
function genrand_int31()
	ccall((:genrand_int31, lib), Cint, ())
end
function genrand_res53()
	ccall((:genrand_res53, lib), Cdouble, ())
end
function dlamch(c::Char)
	ccall((:dlamch, lib), Cdouble, (Cchar, ), c)
end
function get_errno()
    errno
end
function set_errno(n::Integer)
    global errno = Int32(n)
end
end
