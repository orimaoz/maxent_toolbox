// Main functions for maximum entropy toolbox.
// The platform-specific parts are abstracted from these functions so that they can work in more than one environment.
// Ori Maoz, August 2016

#include "EnergyFunction.h"
#include <vector>
#include <mkl.h>
#include "mtrand.h"

void getLogProbability(EnergyFunction * energy, uint32_t npatterns, uint32_t * patterns_in, double * logprobs_out, bool bNormalize = true)
{




}


// Performs a Metropolis-Hasting type MCMC random walk on the state space and returns samples.
// For an n-bit word, n bit flips are performed between each returned sample. This can be scaled up with the "nSeparation" argument.
// The bits are flipped in a sequential order unless specified otherwise by the global bSequentialBits argument.
//
// Input:
//		pModel (in)			- model to generate the samples from
//		nsteps (in)			- number of samples to generate
//		x0 (in/out)			- accepts initial state and returns the final state of the random walk.
//		pOutputSamples (out)- pointer to preallocated array that will contain the output samples, or NULL if we don't want to return samples (for burn-in)
//		nSeparation (in)	- how many samples between every returned sample. Is used to generate less-correlated data.
//		bSequentialBits (in)- true if we want bits to be flipped in a sequential order, false for random order
//
// Returns:  
//		The energy (un-normalized log probability) of the new state.
void runGibbsSampler(EnergyFunction * pModel, uint32_t nsteps, uint32_t * x0, unsigned char * pOutputSamples, uint32_t nSeparation, bool bSequentialBits)
{
	std::vector<uint32_t> current_x;;  // inputed x is the st
	std::vector<uint32_t> proposed_x;
	double current_energy = 0;  // energy of the current state
	double proposed_energy = 0; // energy of the proposed state
	double transition_probability;
	double rand_double;
	uint32_t bit_to_flip = 0;
	MTRand engine_double; // Random number generator. Class is a singleton so it's
						  // ok that it's just sitting on the stack like this
	MTRand_int32 engine_integer; // Same but for integers. It shares the internal state
								 // of the other engine so does not need to be initialized anywhere.


	uint32_t n = pModel->getDim(); // dimension of the data

	// starting point for MCMC walk
	current_x.assign(x0, x0 + n);



	if (!bSequentialBits)
	{
		// Bits to flip are chosen randomly - randomly choose the last bit (because usually we do this
		// only at the end of the loop, and we don't want the first bit to be fixed at zero)
		bit_to_flip = engine_integer() % n;
	}


	// Initial energy for x
	current_energy = pModel->getEnergy(current_x.data());

	for (uint32_t outputIdx = 0; outputIdx < nsteps; outputIdx++)
	{
		// For each sample make as many steps as needed		
		for (uint32_t iteration = 0; iteration < nSeparation; iteration++)
		{

			// Find the energy and bin of the proposed state
			proposed_energy = pModel->propose(bit_to_flip);

			// Transition probability is exponential in the difference in densities
			transition_probability = exp(current_energy - proposed_energy);

			// Randomly choose if to accept or reject the proposal
			rand_double = engine_double();

			uint32_t bAccepted = rand_double < transition_probability;
			current_energy = pModel->accept(bAccepted);

			if (bSequentialBits)
			{
				// Bits are flipped one by one in a sequential order
				bit_to_flip = (bit_to_flip + 1) % n;
			}
			else
			{
				// Bits to flip are chosen randomly
				bit_to_flip = engine_integer() % n;
			}


		}

		if (pOutputSamples) // we return the results only if the output buffer is not NULL
		{
			uint32_t * px = pModel->getX();
			for (unsigned int i = 0; i < n; i++)
			{
				pOutputSamples[outputIdx*n + i] = px[i];
			}
		}
	}

	// Return x as the last state
	uint32_t * px = pModel->getX();
	for (unsigned int i = 0; i < n; i++)
	{
		x0[i] = px[i];
	}

}

void recursiveComputeMarginals(EnergyFunction * pModel, unsigned int curr_bit, double * pMarginals, double & z);


// Returns the marginals of a model (exhaustively computed)
// Input:
//		pModel (in)			- model to generate the samples from
//		pMarginals (out		- preallocated buffer in which the marginals are returned
void getMarginals(EnergyFunction * pModel, double * pMarginals)
{
	// Fix the starting state as all zeros
	int n = pModel->getDim();
	int nfactors = pModel->getNumFactors();
	double z = 0;	// partition function
	
	std::vector<uint32_t> x(pModel->getDim());
	std::fill(x.begin(), x.end(), 0);
	pModel->getEnergy(x.data());

	// make a memory-aligned version of the marginals (on 64-bit boundry for faster operation)
	double * pAlignedMarginals = (double*)mkl_malloc(sizeof(double)*nfactors, 64);
	memset(pAlignedMarginals, 0, nfactors * sizeof(double));

	// Call the recursive part to walk over all the pattern combinations
	recursiveComputeMarginals(pModel, 0, pAlignedMarginals, z);

	// Use the partition function to normalize the results, and also divide by the number
	// of samples
	uint32_t nDims = pModel->getDim();
	double nSamples = pow(2, (double)nDims);
	for (uint32_t i = 0; i < nfactors; i++)
	{
		pMarginals[i] = pAlignedMarginals[i] / z;
	}

	mkl_free(pAlignedMarginals);	// release the aligned buffer
}





// recursively computes the weighted sum of factors for all the subpatterns from (curr_bit) and
// towards the LSB. Also sums up the probabilities (i.e. partition function) for this
// subset of probabilities.
// This function assumes that the state of the model it receives has zero in all the bits
// from curr_bit towards the LSB, and is responsible for returning in the same state.
void recursiveComputeMarginals(EnergyFunction * pModel, unsigned int curr_bit, double * pMarginals, double & z)
{
	if (curr_bit < (pModel->getDim() - 1))
	{

		// recurse when the current bit is zero (and the rest is zero)
		recursiveComputeMarginals(pModel, curr_bit + 1, pMarginals, z);

		// switch the current bit to one
		pModel->propose(curr_bit);
		pModel->accept();

		// recurse when the current bit is one (and the rest is zero too)
		recursiveComputeMarginals(pModel, curr_bit + 1, pMarginals, z);

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
		pModel->accept(pMarginals, prob);

		// switch back to zero
		energy = pModel->propose(curr_bit);
		prob = exp(-energy);
		z += prob;
		pModel->accept(pMarginals, prob);


	}
}
