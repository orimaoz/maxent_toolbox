#!/bin/sh
# Builds the maximum entropy toolbox for MATLAB on unix hosts
# Usage:
#
# "source build_unix.sh host" - builds with maximum optimizations for current host
# "source build_unix.sh compatible" - builds backwards-compatible code for old computers


source compilervars.sh intel64
export COMPILERLIBDIR=/apps/RH6U4/intel/2016/lib/intel64
export MKLLIB=/apps/RH6U4/intel/2016/mkl/lib/intel64
export IPPLIB=/apps/RH6U4/intel/2016/ipp/lib/intel64


if [ "$1" == "" ]; then
	echo
	echo Usage: source build_unix.sh [optimization_type]
	echo where [optimization_type] is either:
	echo       optimized 	  - builds with maximum optimizations for current host 
	echo	   compatible - builds backwards-compatible code for old computers
	return
fi

if [ "$1" == "optimized" ]; then
# The following options file uses the maximum optimizations available for the computer
# used to build the code:
	echo Optimizating for current host
	export OPTIONSFILE="mex_g++_glnxa64_optimized.xml"
fi

if [ "$1" == "compatible" ]; then
# The following options file uses optimizations only up to SSE3, this will result in slower
# code that can run on older computers (circa 2005 and onwards). Uncomment it if this is what you need.
	echo Using backwards-compatible optimizations
	export OPTIONSFILE="mex_g++_glnxa64_compatible.xml"
fi


sh ./buildmex.sh
