
export MKLDIR=/opt/intel/composer_xe_2015/mkl
export IPPDIR=/opt/intel/composer_xe_2015/ipp
export COMPILERDIR=/opt/intel/composer_xe_2015


export MKLINC=/opt/intel/composer_xe_2015/mkl/include
export IPPINC=/opt/intel/composer_xe_2015/ipp/include
export MKLLIB=/opt/intel/composer_xe_2015/mkl/lib
export IPPLIB=/opt/intel/composer_xe_2015/ipp/lib
export COMPILERDIR=/opt/intel/composer_xe_2015/lib
export OPTIONSFILE="mex_C++_maci64.xml"


sh ./buildmex.sh
