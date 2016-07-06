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

#include <algorithm>
#include <iterator>
#include <vector>

#define MAX_STRLEN_MODELTYPE 20
#define DEFAULT_SEPARATION 1

//#define DEBUG_PRINTS

void printVector(std::strstream & str, char* name, double vec[], size_t len);

void recursiveComputeMarginals(EnergyFunction * pModel, unsigned int curr_bit, double * pMarginals, double & z);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	EnergyFunction* pModel;

	if(nrhs<1) {
	    mexErrMsgIdAndTxt("mexGetMarginals:init",
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
	double * pAlignedMarginals= (double*)mkl_malloc(sizeof(double)*nFactors, 64);
	plhs[0] = mxMarginals;




	memset(pAlignedMarginals, 0, nFactors * sizeof(double));

	// Initialize the marginals to zero
	for (uint32_t i=0; i < nFactors; i++)
	{
		pMarginals[i]=0;
	}
	

	// partition function initially zero
	double z=0;

	// Initial the as zero
	std::vector<uint32_t> x(pModel->getDim());
	std::fill(x.begin(), x.end(), 0);

	pModel->getEnergy(x);

	// Sum up the un-normalized marginals
	//recursiveComputeMarginals(pModel,0,pMarginals,z);
	recursiveComputeMarginals(pModel, 0, pAlignedMarginals, z);

	// Use the partition function to normalize the results, and also divide by the number
	// of samples
	uint32_t nDims = pModel->getDim();
	double nSamples = pow(2,(double)nDims);
	for (uint32_t i=0; i < nFactors; i++)
	{
		pMarginals[i] = pAlignedMarginals[i] / z;
	}	


	// also return the partition function
	plhs[1] = mxCreateDoubleScalar(z);

	// Delete the model that we had previously allocated
	delete pModel;
	mkl_free(pAlignedMarginals);
	
}

// recursively computes the weighted sum of factors for all the subpatterns from (curr_bit) and
// towards the LSB. Also sums up the probabilities (i.e. partition function) for this
// subset of probabilities.
// This function assumes that the state of the model it receives has zero in all the bits
// from curr_bit towards the LSB, and is responsible for returning in the same state.
void recursiveComputeMarginals(EnergyFunction * pModel, unsigned int curr_bit, double * pMarginals, double & z)
{
	if (curr_bit < (pModel->getDim()-1))
	{

		// recurse when the current bit is zero (and the rest is zero)
		recursiveComputeMarginals(pModel,curr_bit+1,pMarginals,z);

		// switch the current bit to one
		pModel->propose(curr_bit);
		pModel->accept();

		// recurse when the current bit is one (and the rest is zero too)
		recursiveComputeMarginals(pModel,curr_bit+1,pMarginals,z);

		// switch back to zero
		pModel->propose(curr_bit);
		pModel->accept();
	}
	else
	{
		double energy, prob;
		// we are at the bottom, switch and sum

		// switch the current bit to one
		energy = pModel->propose(curr_bit);
		prob = exp(-energy);
		z += prob;
		pModel->accept(pMarginals,prob);

		// switch back to zero
		energy = pModel->propose(curr_bit);
		prob = exp(-energy);
		z += prob;
		pModel->accept(pMarginals,prob);


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
