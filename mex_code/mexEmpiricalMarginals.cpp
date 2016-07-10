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

	// Get the log probability for each sample
	double z = pModel->getLogZ();
	std::vector<uint32_t> x;
	double * sum_factors = (double*)mkl_malloc(sizeof(double)*nFactors, 64);
	memset(sum_factors, 0, sizeof(double)*nFactors);
	uint32_t offset = 0;
	for (size_t i=0; i < nsamples; i++)
	{
		if ((sampleClassid == mxUINT32_CLASS) || (sampleClassid == mxINT32_CLASS))	
		{
			// fetch the sample
			x.assign((uint32_t*)px + offset,(uint32_t*)px + offset +ndims);
		}
		else if ((sampleClassid == mxCHAR_CLASS) || (sampleClassid == mxLOGICAL_CLASS) || (sampleClassid == mxUINT8_CLASS)|| (sampleClassid == mxINT8_CLASS))
		{
			x.assign((unsigned char*)px + offset,(unsigned char*)px + offset +ndims);
		}
		else if (sampleClassid == mxDOUBLE_CLASS)	
		{
			x.assign((double*)px + offset,(double*)px + offset + ndims);		
		}
		else
		{
			mexErrMsgIdAndTxt("mexEmpiricalMarginals:init",
							  "x is of an unsupported type");
		}

		// get the probability for this factor (if it was given)
		double prob;
		if (bInputedProbabilities)
		{
			prob = pInputProbabilities[i];
		}
		else
		{
			// default is uniform probability. We will perform the division later.
			prob = uniform_prob;
		}
			

		// Compute the factor and and sum it
		pModel->sumSampleFactor(x, sum_factors, prob);

		// advance to next sample
		offset += ndims;
	}

	// copy it to the output array
	for (size_t i=0; i < nFactors;i++)
	{
		pMarginals[i] = sum_factors[i];
	}
	mkl_free(sum_factors);

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
