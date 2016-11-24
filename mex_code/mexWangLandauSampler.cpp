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
#include "mex.h"
#include "mtrand.h"

#include "EnergyFunctionFactory.h"
#include "maxent_functions.h"

#include <algorithm>
#include <iterator>
#include <vector>

static bool bAlreadySeededRand = false;

#define MAX_STRLEN_MODELTYPE 20
#define DEFAULT_SEPARATION 1

//#define DEBUG_PRINTS

void printVector(std::strstream & str, char* name, double vec[], size_t len);
//void runWangLandauStep(uint32_t nsteps, uint32_t * x, EnergyFunction *  pModel, uint32_t nbins, double bin_limits[], double g[], double h[], double update_size, uint32_t nSeparation);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	EnergyFunction* pModel;
	uint32_t nSeparation;

	if(nrhs<4) {
	    mexErrMsgIdAndTxt("mexWangLandauSampler:mexWangLandauSampler",
                      "Usage: mexWangLandauSampler(x0,nsteps,energy_params,model [,separation])");
	}


	// get initial state x0
	uint32_t * ptr_initial_x  = (uint32_t*)mxGetData(prhs[0]);
	size_t nDims = mxGetN(prhs[0]);
	mxClassID id = mxGetClassID(prhs[0]);
	if (id != mxUINT32_CLASS)	
	{
		mexErrMsgIdAndTxt("mexWangLandauSampler:init",
						  "separation must be of type uint32_t");
	}


	// get number of steps
	double * npSteps = mxGetPr(prhs[1]);
	uint32_t nSteps = (uint32_t) *npSteps;


	// Start parsing the structure containing information about the energy
	const mxArray * energy_struct = prhs[2];

	// get energy histogram limits
	mxArray * mxBinLimits = mxGetField(energy_struct,0,"bins");
	if (!mxBinLimits) mexErrMsgIdAndTxt("mexWangLandauSampler:init","cannot find energies.bins");
	double * bin_limits = (double*)mxGetData(mxBinLimits);
	size_t nBins = mxGetN(mxBinLimits);

	// get energy densities
	mxArray * mxG = mxGetField(energy_struct,0,"densities");
	if (!mxG) mexErrMsgIdAndTxt("mexWangLandauSampler:init","cannot find energies.densities");
	double * g = (double*)mxGetData(mxG);
	size_t nEnergies = mxGetN(mxG);

	// get energy histogram
	mxArray * mxH = mxGetField(energy_struct,0,"histogram");
	if (!mxH) mexErrMsgIdAndTxt("mexWangLandauSampler:init","cannot find energies.histogram");
	double * h = (double*)mxGetData(mxH);
	size_t nHistograms = mxGetN(mxH);

	// Check that all three arrays are of the same length
	if (!((nBins == nEnergies) && (nBins == nHistograms )))
	{
		mexErrMsgIdAndTxt("mexWangLandauSampler:init",
						  "all histogram arrays must be of the same length");
	}

	// Get update size (for energy density function)
	mxArray * mxUpdateSize = mxGetField(energy_struct,0,"update_factor");
	if (!mxH) mexErrMsgIdAndTxt("mexWangLandauSampler:init","cannot find energies.update_factor");
	double * pUpdateSize = mxGetPr(mxUpdateSize);
	double updateSize = (double) *pUpdateSize;

	// Process structure containing information about the model
	const mxArray * model_struct = prhs[3];


	EnergyFunctionFactory factory;
	pModel = factory.createEnergyFunction(model_struct);

	// Default separation (if nothing was supplied) is one step per bit
	nSeparation = DEFAULT_SEPARATION;

	// The separation argument is optional (if it's not given, we assume it is equal to the number of bits)
	if (nrhs>=5) {
		id = mxGetClassID(prhs[4]);
		if (id != mxDOUBLE_CLASS)
		{
			mexErrMsgIdAndTxt("mexWangLandauSampler:wangLandauStep",
							  "h must be of type double");
		}	

		double * npSeparation = mxGetPr(prhs[4]);
		nSeparation = (uint32_t) *npSeparation;
	}


	// debug prints of the arguments
#ifdef DEBUG_PRINTS
	std::strstream strOut;
	strOut << "Number of steps: " << (uint32_t)nSteps << "\n";	
	strOut << "Number of dimensions: " << (uint32_t)nDims << "\n";	
	strOut << "Update size: " << updateSize << "\n";
	strOut << "Separation: " << nSeparation << "\n";
	strOut.clear();

	// print lambdas
	//strOut  << "Lambdas:";
	//printVector(strOut,"x",ptr_initial_x,nDims);
	printVector(strOut,"lambdas",pLambdas,nLambdas);
	printVector(strOut,"E",bin_limits,nBins);
	printVector(strOut,"g",g,nBins);
	printVector(strOut,"h",h,nBins);

	mexPrintf(strOut.str());
	mexEvalString("drawnow;"); // force it to update the screen so we can see the text
#endif

	// Seed the random number generator
	if (!bAlreadySeededRand)
	{
		// This flag is in the data area of the DLL and should bge persistent across function calls.
		// If somehow it is not persistent then it should reset to false and then we will simply re-seed the RNG.
		bAlreadySeededRand = true;
		MTRand engine((uint32_t)time(NULL));
	}

	runWangLandauStep(nSteps, ptr_initial_x,pModel,nBins,bin_limits,g,h,updateSize,nSeparation);


	// Delete the model that we had previously allocated
	delete pModel;
	
}




