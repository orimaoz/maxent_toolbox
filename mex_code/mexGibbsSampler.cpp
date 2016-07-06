//
// Generic monte-carlo sampler based on the Metropolis-Hasting algorithm (unlike the misleading name, it does not actually do Gibbs sampling).

// This can work with any generic Boltzmann distribution which implements the EnergyFunction interface.
//

#pragma warning(disable:4996)


#include <vector>
#include <cmath>
#include <strstream>
#include <string.h>
#include <time.h>
#include <mex.h>
#include "mtrand.h"

#include <mkl.h>
#include "EnergyFunction.h"
#include "EnergyFunctionFactory.h"

static bool bAlreadySeededRand = false;
bool bSequentialBits = true;

#define MAX_STRLEN_MODELTYPE 20
#define DEFAULT_SEPARATION 1
#define DEFAULT_BURNIN 10000
//#define DEBUG_PRINTS

#ifdef DEBUG_PRINTS
void printVector(std::strstream & str, char* name, double vec[], size_t len);
#endif

void runGibbsSampler(uint32_t nSteps, std::vector<uint32_t> & initial_x, unsigned char * pOutputSamples, EnergyFunction & model, uint32_t nSeparation,uint32_t bReturnResults);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	EnergyFunction* pModel;
	uint32_t nSeparation = DEFAULT_SEPARATION;
	uint32_t nBurnin = DEFAULT_BURNIN;
	bool b_supplied_x0;
	std::vector<uint32_t> initial_x;
	uint32_t randSeed = 0;

	if(nrhs<3) {
	    mexErrMsgIdAndTxt("GibbsSampler:init",
                      "Usage: gibbsSampler(x0,nsteps,model[,params])\n"
					  "possible fields for params:\n"
					  "burnin - length of initial burn-in (default: 10000)\n"
					  "separation - sample separation (default: 1)\n"
					  "randseed - seed for random number generator (default: from clock)\n"					  
					  "randbits - true if bits to flip are chosen randomly (not sequentially)\n"
					  );
	}


	// get initial state x0
	size_t nDims = mxGetN(prhs[0]);
	mxClassID x0_id = mxGetClassID(prhs[0]);
	if (nDims>1)
	{
		// we were given a starting point for the random walk, extract it
		b_supplied_x0 = true;

		initial_x.resize(nDims,0);

		if (x0_id == mxUINT32_CLASS)
		{
			uint32_t * ptr_initial_x  = (uint32_t*)mxGetData(prhs[0]);
			//copy(initial_x.begin(),initial_x.end(),ptr_initial_x);
			std::copy(ptr_initial_x, ptr_initial_x + nDims, initial_x.begin());

		}
		else if (x0_id == mxDOUBLE_CLASS)
		{
			// Data was supplied as double, copy it to the internal uint32 class
			double * ptr_initial_x  = (double*)mxGetData(prhs[0]);
			for (unsigned int i=0; i < nDims; i++)
			{
				initial_x[i] = (ptr_initial_x[i] > 0);
			}
		}
		else
		{
			mexErrMsgIdAndTxt("GibbsSampler:init",
							  "x must be of type uint32_t or double");
		}
	}
	else
	{
		// Nothing was supplied, use a starting point of all-zeros
		//  (this is actually a moderately-ok starting point for most neural data)
		b_supplied_x0 = false;
	}


	// get number of steps
	double * npSteps = mxGetPr(prhs[1]);
	uint32_t nSteps = (uint32_t) *npSteps;


	// Process structure containing information about the model
	const mxArray * model_struct = prhs[2];

	EnergyFunctionFactory factory;
	pModel = factory.createEnergyFunction(model_struct);


	if (b_supplied_x0)
	{
		// Check that the starting point has the same dimension as supported by the model
		if (nDims != pModel->getDim())
		{
			mexErrMsgIdAndTxt("GibbsSampler:init",
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
		mxArray * samples = mxCreateNumericMatrix(nDims,nSteps, mxUINT8_CLASS, mxREAL);
		pOutputSamples = (unsigned char*)mxGetData(samples);
		plhs[0] = samples;
	} // bReturnResults 

	// If there is a burnin phase, start by cycling through some samples without returning them
	if (nBurnin>0)
	{
		runGibbsSampler(nBurnin, initial_x, pOutputSamples, *pModel, nSeparation,false);
	}

	// Now work on the samples we do want
	runGibbsSampler(nSteps, initial_x, pOutputSamples, *pModel, nSeparation,bReturnResults);


	if (b_supplied_x0)
	{
		// return the last x in the chain back to matlab as the current state
		if (x0_id == mxUINT32_CLASS)
		{
			uint32_t * ptr_initial_x  = (uint32_t*)mxGetData(prhs[0]);
			std::copy(initial_x.begin(),initial_x.end(),ptr_initial_x);

		}
		else if (x0_id == mxDOUBLE_CLASS)
		{
			// Data was supplied as double, copy it to the internal uint32 class
			double * ptr_initial_x  = (double*)mxGetData(prhs[0]);
			for (unsigned int i=0; i < nDims; i++)
			{
				ptr_initial_x[i] = initial_x[i]>0;
			}
		}
	} // b_supplied_x0

	// Delete the model that we had previously allocated
	delete pModel;
	
}


// Performs a Metropolis-Hasting type MCMC random walk on the state space and returns samples.
// For an n-bit word, n bit flips are performed between each returned sample. This can be scaled up with the "nSeparation" argument.
// The bits are flipped in a sequential order unless specified otherwise by the global bSequentialBits argument.
//
// Input:
//		nsteps (in)			- number of samples to generate
//		x (in/out)			- accepts initial state and returns the final state of the random walk.
//		pOutputSamples (out)- pointer to preallocated array that will contain the output samples
//		model (in)			- model to generate the samples from
//		nSeparation (in)	- how many samples between every returned sample. Is used to generate less-correlated data.
//		bReturnResults (in) - If true, pOutputSamples is used to return the results. Otherwise it is ignored and all the samples are thrown away,
//							in this case the function is used to perform burn-in.
//
// Returns:  
//		The energy (un-normalized log probability) of the new state.
void runGibbsSampler(uint32_t nsteps, std::vector<uint32_t> & x, unsigned char * pOutputSamples, EnergyFunction & model, uint32_t nSeparation, uint32_t bReturnResults)
{
	std::vector<uint32_t> current_x(x);  // inputed x is the st
	std::vector<uint32_t> proposed_x;
	double current_energy=0;  // energy of the current state
	double proposed_energy=0; // energy of the proposed state
	double transition_probability;
	double rand_double;
	uint32_t bit_to_flip;
	MTRand engine_double; // Random number generator. Class is a singleton so it's
					// ok that it's just sitting on the stack like this
	MTRand_int32 engine_integer; // Same but for integers. It shares the internal state
				// of the other engine so does not need to be initialized anywhere.

	uint32_t n = x.size(); // dimension of the data

	// Initial energy for x
	current_energy = model.getEnergy(x);

	for (uint32_t outputIdx = 0; outputIdx < nsteps; outputIdx++)
	{
		// For each sample make as many steps as needed		
		for (uint32_t iteration = 0; iteration < nSeparation; iteration++)
		{
			for (uint32_t current_bit = 0; current_bit < n; current_bit++)
			{
				if (bSequentialBits)
				{
					// Bits are flipped one by one in a sequential order
					bit_to_flip = current_bit;
				}
				else
				{
					// Bits to flip are chosen randomly
					bit_to_flip = engine_integer() % x.size();
				}

				// Find the energy and bin of the proposed state
				proposed_energy = model.propose(bit_to_flip);

				// Transition probability is exponential in the difference in densities
				transition_probability = exp(current_energy - proposed_energy);

				// Randomly choose if to accept or reject the proposal
				rand_double = engine_double();

				uint32_t bAccepted = rand_double < transition_probability;
				current_energy = model.accept(bAccepted);
			}
		}

		if (bReturnResults)
		{
			std::vector<uint32_t> * px = model.getX();
			for (unsigned int i = 0; i < n; i++)
			{
				pOutputSamples[outputIdx*n + i] = ((*px)[i] == 1);
			}
		}
	}

	// Return x as the last state
	for (unsigned int i=0; i < n; i++)
	{
		std::vector<uint32_t> * px = model.getX();
		x[i] = ((*px)[i] == 1);
	}

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
