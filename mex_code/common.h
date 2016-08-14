// Maximum entropy toolbox
// This is definitions to help compile with or without the Intel compiler and MKL/IPP libraries
// Ori Maoz, August 2016


#pragma once

// comment these out if you don't want MKL and IPP
//#include <mkl.h>
//#include <ipps.h>
//#include <ippvm.h>


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

//
//void ippsAdd_64f_I(double * src, double * dst, int len)
//{
//	for (int j = 0; j < len; j++)
//	{
//		dst[j] += src[j];
//	}
//}
//
//
//void ippsSub_64f_I(double * src, double * dst, int len)
//{
//	for (int j = 0; j < len; j++)
//	{
//		dst[j] -= src[j];
//	}
//}



#endif // NO_INTEL_COMPILER