// Computes the marginals of a set of samples according to the inputed model.
//

#pragma warning(disable:4996)


#include <cmath>
#include <strstream>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include "mex.h"
#include "mtrand.h"
#include <mkl.h>

#include "IsingEnergy.h"
#include "EnergyFunctionFactory.h"
#include "maxent_functions.h"

#include <algorithm>
#include <iterator>
#include <vector>

#define MAX_STRLEN_MODELTYPE 20
#define DEFAULT_SEPARATION 1

//#define DEBUG_PRINTS

void printVector(std::strstream & str, char* name, double vec[], size_t len);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	EnergyFunction* pModel;

	if(nrhs<1) {
	    mexErrMsgIdAndTxt("mexExhaustiveMarginals:init",
                      "Usage: mexExhaustiveMarginals(model)");
	}


	// Process structure containing information about the model
	const mxArray * model_struct = prhs[0];

	EnergyFunctionFactory factory;
	pModel = factory.createEnergyFunction(model_struct);



	// debug prints of the arguments
#ifdef DEBUG_PRINTS
	std::strstream strOut;
	strOut << "Number of samples: " << (uint32_t)nsamples << "\n";	
	strOut << "Number of dimensions: " << (uint32_t)ndims << "\n";	
	strOut.clear();

	// print lambdas
	//printVector(strOut,"x",px,ndims);

	mexPrintf(strOut.str());
	mexEvalString("drawnow;"); // force it to update the screen so we can see the text
#endif


	// Allocate a result array
	uint32_t nFactors = pModel->getNumFactors();
	mxArray * mxMarginals= mxCreateNumericMatrix(1,nFactors, mxDOUBLE_CLASS, mxREAL);
	double * pMarginals = (double*)mxGetData(mxMarginals);
	plhs[0] = mxMarginals;

	getMarginals(pModel, pMarginals);
	
	// Delete the model that we had previously allocated
	delete pModel;	
}


// Prints the beginning and end of a vector to the display
void printVector(std::strstream & str, char* name, double vec[], size_t len)
{

	str  << name << ": [ ";
	for (size_t i=0; i < std::min<size_t>(len,4); i++)
	{
		str << vec[i] << " ";
	}

	if (len>1)
	{
		str  << "... " << vec[len-1] << " ]\n";
	}
}
