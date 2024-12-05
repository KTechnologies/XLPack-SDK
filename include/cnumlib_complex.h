/****************************************
 *                                      *
 *  XLPack Numerical Library            *
 *  Version 7.0 (January 31, 2023)      *
 *  (C) 2014-2023  K Technologies       *
 *                                      *
 ****************************************/
#pragma once

#if defined(__cplusplus)

#include <complex>
#define doublecomplex	std::complex<double>
#define floatcomplex	std::complex<float>

#else

#ifndef _CRT_USE_C_COMPLEX_H
	#define _CRT_USE_C_COMPLEX_H
	#include <complex.h>
	#undef _CRT_USE_C_COMPLEX_H
#else
	#include <complex.h>
#endif
#if defined(__INTEL_LLVM_COMPILER) || (defined(__INTEL_COMPILER) && defined(__STDC_VERSION__)) || defined(__clang__) || defined(__GNUC__) || defined(__APPLE__)
	#define doublecomplex	double _Complex
	#define floatcomplex	float _Complex
	#if !defined(CMPLX)
		#define CMPLX(r, i) ((double _Complex)((double)(r) + _Complex_I * (double)(i)))
	#endif
	#if !defined(CMPLXF)
		#define CMPLXF(r, i) ((float _Complex)((float)(r) + _Complex_I * (float)(i)))
	#endif
#elif !defined(__INTEL_COMPILER) && defined(_MSC_VER)
	#define doublecomplex	_Dcomplex
	#define floatcomplex	_Fcomplex
	#if !defined(CMPLX)
		#define CMPLX(r,i) _Cbuild(r, i)
	#endif
	#if !defined(CMPLXF)
		#define CMPLXF(r,i) _FCbuild(r, i)
	#endif
#else
	typedef struct { double r, i; } doublecomplex;
	typedef struct { float r, i; } floatcomplex;
	#if !defined(CMPLX)
		#define CMPLX(r,i) _cmplx(r, i)
	#endif
	#if !defined(CMPLXF)
		#define CMPLXF(r,i) _cmplxf(r, i)
	#endif
#endif

#endif
