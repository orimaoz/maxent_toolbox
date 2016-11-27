
DIRECTIVES="CXX=icpc LD=icpc"
SUPPORT_FILES=$(find *.cpp -not -name 'mex*.cpp')
INCLUDES="-I$CPATH"
LIBRARIES="$MKLLIB/libmkl_intel_lp64.a $MKLLIB/libmkl_core.a $MKLLIB/libmkl_sequential.a $IPPLIB/libipps.a $IPPLIB/libippcore.a $COMPILERLIBDIR/libirc.a $COMPILERLIBDIR/libimf.a $COMPILERLIBDIR/libirng.a"
FLAGS="-D_IPP_SEQUENTIAL_STATIC -v -f $OPTIONSFILE"

echo Compiling gibbsSampler
mex mexGibbsSampler.cpp $DIRECTIVES $SUPPORT_FILES $LIBRARIES $INCLUDES $FLAGS
echo Compiling mexLogProbability
mex mexLogProbability.cpp $DIRECTIVES $SUPPORT_FILES $INCLUDES $LIBRARIES $FLAGS
echo Compiling mexGetMarginals
mex mexEmpiricalMarginals.cpp $DIRECTIVES $SUPPORT_FILES $INCLUDES $LIBRARIES $FLAGS
echo Compiling mexGetExhaustiveMarginals
mex mexExhaustiveMarginals.cpp $DIRECTIVES $SUPPORT_FILES $INCLUDES $LIBRARIES $FLAGS
echo Compiling wangLandauSampler
mex mexWangLandauSampler.cpp $DIRECTIVES $SUPPORT_FILES $INCLUDES $LIBRARIES $FLAGS

