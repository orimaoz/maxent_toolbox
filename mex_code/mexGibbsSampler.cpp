//
// Generic monte-carlo sampler based on the Metropolis-Hasting algorithm (unlike the misleading name, it does not actually do Gibbs sampling).

// This can work with any generic Boltzmann distribution which implements the EnergyFunction interface.
// 
// Last update: 
// Ori Maoz, July 2016
//

#pragma warning(disable:4996)


#include <vector>
#include <cmath>
#include <strstream>
#include <string.h>
#include <time.h>
#include <mex.h>
#include "mtrand.h"

#include "EnergyFunction.h"
#include "EnergyFunctionFactory.h"

#include "maxent_functions.h"

static bool bAlreadySeededRand = false;
bool bSequentialBits = true;

#define MAX_STRLEN_MODELTYPE 20
#define DEFAULT_BURNIN 10000
//#define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
void printVector(std::strstream & str, char* name, double vec[], size_t len);
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	EnergyFunction* pModel;
	uint32_t nSeparation;
	uint32_t nBurnin = DEFAULT_BURNIN;
	bool b_supplied_x0;
	std::vector<uint32_t> initial_x;
	uint32_t randSeed = 0;

	if(nrhs<3) {
	    mexErrMsgIdAndTxt("mexGibbsSampler:init",
					  "Usage: mexGibbsSampler(model,nsamples,x0,[,params])\n"
					  "possible fields for params structure:\n"
					  "burnin - length of initial burn-in (default: 10000)\n"
					  "separation - sample separation (default: 1)\n"
					  "randseed - seed for random number generator (default: from clock)\n"					  
					  "randbits - true if bits to flip are chosen randomly (not sequentially)\n"
					  );
	}



	// Process structure containing information about the model
	const mxArray * model_struct = prhs[0];

	// Initialize the model according to the string that we got
	EnergyFunctionFactory factory;
	pModel = factory.createEnergyFunction(model_struct);

	// get number of samples that we need to generate
	double * npSamples = mxGetPr(prhs[1]);
	uint32_t nSamples = (uint32_t)*npSamples;


	// get initial state x0
	size_t nDims = mxGetN(prhs[2]);
	mxClassID x0_id = mxGetClassID(prhs[2]);
	if (nDims>1)
	{
		// we were given a starting point for the random walk, extract it
		b_supplied_x0 = true;

		initial_x.resize(nDims,0);

		if (x0_id == mxUINT32_CLASS)
		{
			uint32_t * ptr_initial_x  = (uint32_t*)mxGetData(prhs[2]);
			//copy(initial_x.begin(),initial_x.end(),ptr_initial_x);
			std::copy(ptr_initial_x, ptr_initial_x + nDims, initial_x.begin());

		}
		else if (x0_id == mxDOUBLE_CLASS)
		{
			// Data was supplied as double, copy it to the internal uint32 class
			double * ptr_initial_x  = (double*)mxGetData(prhs[2]);
			for (unsigned int i=0; i < nDims; i++)
			{
				initial_x[i] = (ptr_initial_x[i] > 0);
			}
		}
		else
		{
			mexErrMsgIdAndTxt("mexGibbsSampler:init",
							  "x must be of type uint32_t or double");
		}
	}
	else
	{
		// Nothing was supplied, use a starting point of all-zeros
		//  (this is actually a moderately-ok starting point for most neural data)
		b_supplied_x0 = false;
	}


	if (b_supplied_x0)
	{
		// Check that the starting point has the same dimension as supported by the model
		if (nDims != pModel->getDim())
		{
			mexErrMsgIdAndTxt("mexGibbsSampler:init",
							  "x0 and model are of different dimensions");
		}
	}
	else
	{
		// Starting point was not supplied, we will create it now that we know the
		// dimension that we are working in. We will start with an initial point
		// of all-zeros which is an appropriate guess for much of the neural spiking data.
		nDims = pModel->getDim();
		initial_x.resize(nDims,0);
		std::fill(initial_x.begin(),(initial_x.begin() + nDims),0);
	}

	// Default separation (if nothing was supplied) is one step per bit
	nSeparation = pModel->getDim();


	// Start parsing the structure containing optional parameters
	if (nrhs >= 4)
	{
		const mxArray * params_struct = prhs[3];
		if (params_struct )
		{
			// Get number of burnin steps
			mxArray * mxBurnin = mxGetField(params_struct,0,"burnin");
			if (mxBurnin)
			{
				double * npBurnin= mxGetPr(mxBurnin);
				nBurnin = (uint32_t) * npBurnin;
			}

			// Get length of separation
			mxArray * mxSeparation = mxGetField(params_struct,0,"separation");
			if (mxSeparation)
			{
				double * npSeparation= mxGetPr(mxSeparation);
				nSeparation = (uint32_t) * npSeparation;
			}


			// Get optional RNG seed
			mxArray * mxRandSeed = mxGetField(params_struct, 0, "randseed");
			if (mxRandSeed)
			{
				// if user supplied a seed, save it and force re-seed
				randSeed = (uint32_t)mxGetScalar(mxRandSeed);
				bAlreadySeededRand = false;
			}

			// Get optional flag to choose changed bits randomly (rather than sequentially)
			mxArray * mxRandBits = mxGetField(params_struct, 0, "randbits");
			if (mxRandBits)
			{
				// if user supplied a seed, save it and force re-seed
				unsigned char randbits_flag = (unsigned char)mxGetScalar(mxRandBits);
				if (randbits_flag)
				{
					bSequentialBits = false;
				}
				
			}


		} // params_struct
	} // nrhs>=4


	// debug prints of the arguments
#ifdef DEBUG_PRINTS
	std::strstream strOut;
	strOut << "Number of steps: " << (uint32_t)nSteps << "\n";	
	strOut << "Number of dimensions (input): " << (uint32_t)nDims << "\n";	
	strOut << "Number of dimensions (model): " << (uint32_t)pModel->getDim() << "\n";	
	strOut << "Burnin: " << nBurnin << "\n";
	strOut << "Separation: " << nSeparation << "\n";
	strOut << "Output arguments: " << nlhs << "\n";
	strOut.clear();

	// print lambdas
	//strOut  << "Lambdas:";
//	printVector(strOut,"x",ptr_initial_x,nDims);
//	printVector(strOut,"lambdas",pLambdas,nLambdas);

	mexPrintf(strOut.str());
	mexEvalString("drawnow;"); // force it to update the screen so we can see the text
#endif

	// Seed the random number generator
	if (!bAlreadySeededRand)
	{
		// This flag is in the data area of the DLL and should bge persistent across function calls.
		// If somehow it is not persistent then it should reset to false and then we will simply re-seed the RNG.
		bAlreadySeededRand = true;
		if (randSeed != 0)
		{
			// Use user-supplied seed (mostly used for debugging)
			MTRand engine(randSeed);
		}
		else
		{
			// Use system time for seed
			MTRand engine((uint32_t)time(NULL));
		}
	}

	bool bReturnResults = (nlhs>0);  // does the user want us to return anything

	// Allocate an output array
	unsigned char * pOutputSamples = NULL;
	if (bReturnResults)
	{
		mxArray * samples = mxCreateNumericMatrix(nDims,nSamples, mxUINT8_CLASS, mxREAL);
		pOutputSamples = (unsigned char*)mxGetData(samples);
		plhs[0] = samples;
	} // bReturnResults 

	// If there is a burnin phase, start by cycling through some samples without returning them
	if (nBurnin>0)
	{
		runGibbsSampler(pModel, nBurnin, initial_x.data(), NULL, nSeparation, bSequentialBits);
	}

	// Now work on the samples we do want
	runGibbsSampler(pModel, nSamples, initial_x.data(), pOutputSamples, nSeparation, bSequentialBits);


	if (b_supplied_x0)
	{
		// return the last x in the chain back to matlab as the current state
		if (x0_id == mxUINT32_CLASS)
		{
			uint32_t * ptr_initial_x  = (uint32_t*)mxGetData(prhs[2]);
			std::copy(initial_x.begin(),initial_x.end(),ptr_initial_x);

		}
		else if (x0_id == mxDOUBLE_CLASS)
		{
			// Data was supplied as double, copy it to the internal uint32 class
			double * ptr_initial_x  = (double*)mxGetData(prhs[2]);
			for (unsigned int i=0; i < nDims; i++)
			{
				ptr_initial_x[i] = initial_x[i]>0;
			}
		}
	} // b_supplied_x0

	// Delete the model that we had previously allocated
	delete pModel;
	
}

#ifdef DEBUG_PRINTS
// Prints the beginning and end of a vector to the display
void printVector(std::strstream & str, char* name, double vec[], size_t len)
{

	str  << name << ": [ ";
	for (size_t i=0; i < std::min<size_t>(len,4); i++)
	{
		str << vec[i] << " ";
	}
	str  << "... " << vec[len-1] << " ]\n";
}
#endif
