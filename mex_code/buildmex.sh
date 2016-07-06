

DIRECTIVES="CXX=icpc LD=icpc"
SUPPORT_FILES="mtrand.cpp IsingEnergy.cpp  MerpSparseEnergy.cpp MerpFastEnergy.cpp KSyncEnergy.cpp KIsingEnergy.cpp IndependentEnergy.cpp EnergyFunctionFactory.cpp"
INCLUDES="-I$MKLINC -I$IPPINC"
LIBRARIES="$MKLLIB/libmkl_intel_lp64.a $MKLLIB/libmkl_core.a $MKLLIB/libmkl_sequential.a $IPPLIB/libipps.a $IPPLIB/libippcore.a $COMPILERDIR/libirc.a $COMPILERDIR/libimf.a"
FLAGS="-D_IPP_SEQUENTIAL_STATIC -v"

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

