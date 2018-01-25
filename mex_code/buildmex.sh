
DIRECTIVES="CXX=icpc LD=icpc"
SUPPORT_FILES=$(find *.cpp -not -name 'mex*.cpp')
INCLUDES="-I$CPATH"
LIBRARIES=""
FLAGS="-v -f $OPTIONSFILE"

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

