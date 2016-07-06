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

#include <algorithm>
#include <iterator>
#include <vector>

static bool bAlreadySeededRand = false;

#define MAX_STRLEN_MODELTYPE 20
#define DEFAULT_SEPARATION 1

//#define DEBUG_PRINTS

void printVector(std::strstream & str, char* name, double vec[], size_t len);
void runWangLandauStep(uint32_t nsteps, std::vector<uint32_t> & x, EnergyFunction & model, uint32_t nbins, double bin_limits[], double g[], double h[], double update_size, uint32_t nSeparation);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	EnergyFunction* pModel;
	uint32_t nSeparation;

	if(nrhs<5) {
		nSeparation = DEFAULT_SEPARATION;
	}

	if(nrhs<4) {
	    mexErrMsgIdAndTxt("WangLandauSampler:WangLandauSampler",
                      "Usage: WangLandauSampler(x0,nsteps,energy_params,model [,separation])");
	}


	// get initial state x0
	uint32_t * ptr_initial_x  = (uint32_t*)mxGetData(prhs[0]);
	size_t nDims = mxGetN(prhs[0]);
	mxClassID id = mxGetClassID(prhs[0]);
	if (id != mxUINT32_CLASS)	
	{
		mexErrMsgIdAndTxt("Ori:wangLandauStep",
						  "separation must be of type uint32_t");
	}
	std::vector<uint32_t> initial_x(ptr_initial_x,ptr_initial_x+nDims);


	// get number of steps
	double * npSteps = mxGetPr(prhs[1]);
	uint32_t nSteps = (uint32_t) *npSteps;


	// Start parsing the structure containing information about the energy
	const mxArray * energy_struct = prhs[2];

	// get energy histogram limits
	mxArray * mxBinLimits = mxGetField(energy_struct,0,"bins");
	if (!mxBinLimits) mexErrMsgIdAndTxt("Ori:wangLandauGeneric","cannot find energies.bins");
	double * bin_limits = (double*)mxGetData(mxBinLimits);
	size_t nBins = mxGetN(mxBinLimits);

	// get energy densities
	mxArray * mxG = mxGetField(energy_struct,0,"densities");
	if (!mxG) mexErrMsgIdAndTxt("Ori:wangLandauGeneric","cannot find energies.densities");
	double * g = (double*)mxGetData(mxG);
	size_t nEnergies = mxGetN(mxG);

	// get energy histogram
	mxArray * mxH = mxGetField(energy_struct,0,"histogram");
	if (!mxH) mexErrMsgIdAndTxt("Ori:wangLandauGeneric","cannot find energies.histogram");
	double * h = (double*)mxGetData(mxH);
	size_t nHistograms = mxGetN(mxH);

	// Check that all three arrays are of the same length
	if (!((nBins == nEnergies) && (nBins == nHistograms )))
	{
		mexErrMsgIdAndTxt("Ori:wangLandauStep",
						  "all histogram arrays must be of the same length");
	}

	// Get update size (for energy density function)
	mxArray * mxUpdateSize = mxGetField(energy_struct,0,"update_factor");
	if (!mxH) mexErrMsgIdAndTxt("Ori:wangLandauGeneric","cannot find energies.update_factor");
	double * pUpdateSize = mxGetPr(mxUpdateSize);
	double updateSize = (double) *pUpdateSize;

	// Process structure containing information about the model
	const mxArray * model_struct = prhs[3];


	EnergyFunctionFactory factory;
	pModel = factory.createEnergyFunction(model_struct);

	// The separation argument is optional (if it's not given, we assume it is equal 1)
	if (nrhs>=5) {
		id = mxGetClassID(prhs[4]);
		if (id != mxDOUBLE_CLASS)
		{
			mexErrMsgIdAndTxt("Ori:wangLandauStep",
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

	runWangLandauStep(nSteps,initial_x,*pModel,nBins,bin_limits,g,h,updateSize,nSeparation);

	// return the last x in the chain back to matlab
	// TODO: is this call safe?
	copy(initial_x.begin(),initial_x.end(),ptr_initial_x);

	// Delete the model that we had previously allocated
	delete pModel;
	
}


// Actually runs the Wang-Landau steps of random walk on the energy function
// nsteps - number of steps
// x - starting point (in/out - returns ending state)
// model - model that we use to compute energy
// bin_limits - bin limits for the energy function discretization
// g - energy density function (discretized)
// h - energy function histogram
void runWangLandauStep(uint32_t nsteps, std::vector<uint32_t> & x, EnergyFunction& model, uint32_t nbins, double bin_limits[], double g[], double h[], double update_size, uint32_t nSeparation)
{
	std::vector<uint32_t> current_x(x);  // inputed x is the st
	std::vector<uint32_t> proposed_x;
	double current_energy;  // energy of the current state
	double proposed_energy; // energy of the proposed state
	uint32_t current_bin;		// bin of the current energy
	uint32_t proposed_bin;		// bin of the proposed energy
	double transition_probability;
	double rand_double;
	uint32_t current_bit;
	MTRand engine_double; // Random number generator. Class is a singleton so it's
					// ok that it's just sitting on the stack like this
	MTRand_int32 engine_integer; // Same but for integers. It shares the internal state
				// of the other engine so does not need to be initialized anywhere.


	// compute energy and parameters for initial x
	current_energy = model.getEnergy(x);
	current_bin = std::lower_bound(bin_limits,bin_limits+nbins-1,current_energy) - bin_limits;

	// iterate...
	for (uint32_t iteration=0; iteration < nsteps; iteration++)
	{
		
		// Choose a random bit to flip
		current_bit = engine_integer() % x.size();


		// Find the energy and bin of the proposed state
		proposed_energy = model.propose(current_bit);
		proposed_bin = std::lower_bound(bin_limits,bin_limits+nbins-1,proposed_energy) - bin_limits;

		// Transition probability is exponential in the difference in densities
		transition_probability = exp(g[current_bin] - g[proposed_bin]);

		// Randomly choose if to accept or reject the proposal
		rand_double = engine_double();
		if (rand_double < transition_probability)
		{

			// Accept the proposed x and update everything						
			model.accept();
			current_energy = proposed_energy;
			current_bin = proposed_bin;
		}

		// update density and histogram
		g[current_bin] += update_size;
		
		// If we have a separation value, we only advance the histogram every several steps
		// in order to decorrelate the samples
		if (nSeparation>1)
		{
			if (iteration % nSeparation == 0)
				h[current_bin]++;
		}
		else
		{
			h[current_bin]++;
		}


	}

	// Return x as the last state
	x = (*model.getX());
}



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
