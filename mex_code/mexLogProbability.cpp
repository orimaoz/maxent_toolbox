// Computes the log probabilities of a set of samples on the inputed model.
// Setting up the model may take a while so this is mainly effective when a set of samples
// is inputed and not just a single sample
//
// Generic Wang-Landau sampler for estimating the energy density of a Boltzmann distribution.
// Performs a certain amount of random-walk steps on a weighted energy function and returns
// ending state, accumulated histogram and updated density function.
//
// This can work with any generic Boltzmann distribution which implements the proper (TBD) interfaces.
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

	if(nrhs<2) {
	    mexErrMsgIdAndTxt("mexLogProbability:init",
                      "Usage: mexLogProbability(x,model)");
	}


	// get the set of inputs
	uint32_t * px = (uint32_t*)mxGetData(prhs[0]);
	size_t nsamples = mxGetN(prhs[0]);
	size_t ndims = mxGetM(prhs[0]);
	mxClassID sampleClassid = mxGetClassID(prhs[0]);

	// Process structure containing information about the model
	const mxArray * model_struct = prhs[1];

	EnergyFunctionFactory factory;
	pModel = factory.createEnergyFunction(model_struct);


	if (pModel->getDim() != ndims) {
		mexErrMsgIdAndTxt("mexLogProbability:init",
			"data differs in dimension from model");
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
	mxArray * mxLogprobs= mxCreateNumericMatrix(1,nsamples, mxDOUBLE_CLASS, mxREAL);
	double * pLogprobs = (double*)mxGetData(mxLogprobs);
	plhs[0] = mxLogprobs;

	// Get the log probability for each sample
	double z = pModel->getLogZ();
	std::vector<uint32_t> x;
	uint32_t offset = 0;
	for (size_t i = 0; i < nsamples; i++)
	{
		if ((sampleClassid == mxUINT32_CLASS) || (sampleClassid == mxINT32_CLASS))
		{
			// fetch the sample
			x.assign((uint32_t*)px + offset, (uint32_t*)px + offset + ndims);
		}
		else if ((sampleClassid == mxCHAR_CLASS) || (sampleClassid == mxLOGICAL_CLASS) || (sampleClassid == mxUINT8_CLASS) || (sampleClassid == mxINT8_CLASS))
		{
			x.assign((unsigned char*)px + offset, (unsigned char*)px + offset + ndims);
		}
		else if (sampleClassid == mxDOUBLE_CLASS)
		{
			x.assign((double*)px + offset, (double*)px + offset + ndims);
		}
		else
		{
			mexErrMsgIdAndTxt("mexLogProbability:init",
				"x is of an unsupported type");
		}


		// compute the log probability (negative energy plus log partition function)
		pLogprobs[i] = -pModel->getEnergy(x) + z;


		// advance to next sample
		offset += ndims;
	}


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
