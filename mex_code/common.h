// Maximum entropy toolbox
// This is definitions to help compile with or without the Intel compiler and MKL/IPP libraries
// Ori Maoz, August 2016


#pragma once
#include <cstdlib>



// comment these out if you don't want MKL and IPP
#include <mkl.h>
#include <ipps.h>
#include <ippvm.h>


#if defined(_MKL_H_) || defined(__IPPDEFS_H__)


#ifdef _WINDOWS
#define DECLARE_ALIGNED _declspec(align(64))
#else
#define DECLARE_ALIGNED
#endif

#define malloc_aligned ippsMalloc_64f
#define free_aligned ippsFree

#else
// not using the intel compiler

#define malloc_aligned malloc
#define free_aligned free
#define MKL_INT int


#define DECLARE_ALIGNED


inline void ippsAdd_64f_I(double * src, double * dst, int len)
{
	for (int j = 0; j < len; j++)
	{
		dst[j] += src[j];
	}
}


inline void ippsSub_64f_I(double * src, double * dst, int len)
{
	for (int j = 0; j < len; j++)
	{
		dst[j] -= src[j];
	}
}

inline void ippsDotProd_64f(double *a, double *b, int len, double * sum)
{
	for (int i = 0; i < len; i++)
	{
		*sum += a[i] * b[i];
	}
}

inline void vdAdd(int n, double *a, double *b, double *y)
{
	for (uint32_t j = 0; j < n; j++)
	{
		y[j] = a[j] + b[j];
	}

}


inline void vdSub(int n, double *a, double *b, double *y)
{
	for (uint32_t j = 0; j < n; j++)
	{
		y[j] = a[j] - b[j];
	}
}


inline void cblas_dgthr(MKL_INT nz, double * y, double * x, MKL_INT * idx)
{
	for (uint32_t j = 0; j < nz; j++)
	{
		x[j] = y[idx[j]];
	}

}

inline void cblas_daxpyi(MKL_INT nz, double a, double * x, MKL_INT *idx, double * y)
{
	for (uint32_t j = 0; j < nz; j++)
	{
		y[j] += a*x[idx[j]];
	}
}

#endif // NO_INTEL_COMPILER