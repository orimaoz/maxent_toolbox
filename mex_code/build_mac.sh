#!/bin/sh
# Builds the maximum entropy toolbox for MATLAB on MacOS hosts
# Usage:
#
# "source build_mac.sh host" - builds with maximum optimizations for current host
# "source build_mac.sh compatible" - builds backwards-compatible code for old computers


source /opt/intel/compilers_and_libraries_2017/mac/bin/compilervars.sh intel64
export COMPILERLIBDIR=/opt/intel/compilers_and_libraries_2017/mac/lib

if [ "$1" == "" ]; then
	echo
	echo Usage: source build_mac.sh [optimization_type]
	echo where [optimization_type] is either:
	echo       optimized 	  - builds with maximum optimizations for current host 
	echo	   compatible - builds backwards-compatible code for old computers
	return
fi

if [ "$1" == "optimized" ]; then
# The following options file uses the maximum optimizations available for the computer
# used to build the code:
	echo Optimizating for current host
	export OPTIONSFILE="mex_C++_maci64_optimized.xml"
fi

if [ "$1" == "compatible" ]; then
# The following options file uses optimizations only up to SSE3, this will result in slower
# code that can run on older computers (circa 2005 and onwards). Uncomment it if this is what you need.
	echo Using backwards-compatible optimizations
	export OPTIONSFILE="mex_C++_maci64_compatible.xml"
fi

export MKLLIB=$MKLROOT/lib
export IPPLIB=$IPPROOT/lib


sh ./buildmex.sh

