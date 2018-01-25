#!/bin/sh
# Builds the maximum entropy toolbox for MATLAB on MacOS hosts
# Usage:
#
# "source build_mac.sh" - builds with the mex code


source /opt/intel/compilers_and_libraries/mac/bin/compilervars.sh intel64


# The following options file creates two execution branches, one with AVX2 and one with SSE4.2 for older processors
export OPTIONSFILE="mex_C++_maci64.xml"


sh ./buildmex.sh

