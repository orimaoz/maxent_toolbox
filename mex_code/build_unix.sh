#!/bin/sh
# Builds the maximum entropy toolbox for MATLAB on unix hosts
# Usage:
#
# "source build_unix.sh" - builds the MEX functions


source compilervars.sh intel64
export COMPILERLIBDIR=/apps/RH6U4/intel/2016/lib/intel64
export MKLLIB=/apps/RH6U4/intel/2016/mkl/lib/intel64
export IPPLIB=/apps/RH6U4/intel/2016/ipp/lib/intel64


# The following options file uses two execution paths, one with AVX2 for newer CPUs and one with SSE4.2 for older CPUs
export OPTIONSFILE="mex_g++_glnxa64.xml"

sh ./buildmex.sh
