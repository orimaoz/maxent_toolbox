// Computes the marginals of a set of samples according to the inputed model.
//

#pragma warning(disable:4996)


#include <cmath>
#include <strstream>
#include <string.h>
#include <time.h>
#include <stdint.h>
#include "mex.h"
#include "matlab_utils.h"
#include "maxent_functions.h"

#include "EnergyFunctionFactory.h"
#include <mkl.h>

#include <algorithm>
#include <iterator>
#include <vector>

#define DEFAULT_SEPARATION 1

//#define DEBUG_PRINTS

void printVector(std::strstream & str, char* name, double vec[], size_t len);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	EnergyFunction* pModel;

	if(nrhs<2) {
	    mexErrMsgIdAndTxt("mexEmpiricalMarginals:init",
                      "Usage: mexEmpiricalMarginals(x,model,[probabilities])");
	}

	mkl_disable_fast_mm(); // make MKL use simple and safe memory management


	// get the set of inputs
	void * px = (uint32_t*)mxGetData(prhs[0]);
	size_t nsamples = mxGetN(prhs[0]);
	size_t ndims = mxGetM(prhs[0]);
	mxClassID sampleClassid = mxGetClassID(prhs[0]);


	// Process structure containing information about the model
	const mxArray * model_struct = prhs[1];

	EnergyFunctionFactory factory;
	pModel = factory.createEnergyFunction(model_struct);


	if (pModel->getDim() != ndims)
	{
		mexErrMsgIdAndTxt("mexEmpiricalMarginals:init",
			"model dimension does not match sample dimension");	
	}


	// Check if the probabilities of patterns were inputed.
	// If they were not, treat them as all ones.
	bool bInputedProbabilities = false;
	double * pInputProbabilities = NULL;
	double uniform_prob = (double)1/nsamples;
	if(nrhs>=3) {
		bInputedProbabilities = true;

		const mxArray * mxProbabilities = prhs[2];

		// Check the datatype - probabilities must be zero
		mxClassID probsClassid = mxGetClassID(mxProbabilities);
		if (probsClassid != mxDOUBLE_CLASS)	
		{
			mexErrMsgIdAndTxt("mexEmpiricalMarginals:init",
							  "pattern probabilities must be of type double");
		}

		// Check that it is the same length as the number of samples. We will allow it to be
		// either a row or column vector (because it is transparent to us anyway)
		if (!(((mxGetN(mxProbabilities) == nsamples) && (mxGetM(mxProbabilities) == 1)) || 
			((mxGetM(mxProbabilities) == nsamples) && (mxGetN(mxProbabilities) == 1))))
		{
			mexErrMsgIdAndTxt("mexEmpiricalMarginals:init",
							  "number of probabilities not equal to number of samples");
		}


		pInputProbabilities = (double*)mxGetData(mxProbabilities);
	}



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
	
	// make sure that the datatype is correct
	uint32_t * pInputArray = reallocate_uint32(px, sampleClassid, nsamples*ndims);

	// build a temp output array that is aligned to 64-bit boundaries
	double * sum_factors = (double*)mkl_malloc(sizeof(double)*nFactors, 64);

	// run the main function and get the marginals
	getEmpiricalMarginals(pModel, nsamples, pInputArray, pInputProbabilities, sum_factors);


	// copy the result to the output array and release the temp array
	memcpy(pMarginals, sum_factors, sizeof(double) * nFactors);
	mkl_free(sum_factors);

	// Delete the model that we had previously allocated
	delete pModel;

	// if we needed to translate the input datatype, delete that temporary buffer too
	if (pInputArray != px)
	{
		delete[] pInputArray;
	}

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
