// Main functions for maximum entropy toolbox.
// The platform-specific parts are abstracted from these functions so that they can work in more than one environment.
// Ori Maoz, August 2016

#include "EnergyFunction.h"

// Returns the log probabilities of samples according to the model
// Input:
//		pModel (in)			- model to generate the samples from
//	    npatterns (in)		- number of patterns
//		patterns_in (in)	- samples in UINT32 form (32 bit integer)
//		logprobs_out (out)	- preallocated buffer in which the log probabilities are returned
void getLogProbability(EnergyFunction * pModel, uint32_t npatterns, uint32_t * patterns_in, double * logprobs_out);


// Returns the empirical marginals for a set of samples
// Input:
//		pModel (in)			- model to generate the samples from
//	    npatterns (in)		- number of patterns
//		patterns_in (in)	- samples in UINT32 form (32 bit integer)
//		weights (in)		- probabilites used to reweight the patterns or NULL to treat them as uniform
//		pMarginals (out)	- preallocated buffer in which the marginals are returned
void getEmpiricalMarginals(EnergyFunction * pModel, uint32_t npatterns, uint32_t * patterns_in, double * weights, double * pMarginals);


// Performs a Metropolis-Hasting type MCMC random walk on the state space and returns samples.
// For an n-bit word, n bit flips are performed between each returned sample. This can be changed with the "nSeparation" argument.
// The bits are flipped in a sequential order unless specified otherwise by the global bSequentialBits argument.
//
// Input:
//		pModel (in)			- model to generate the samples from
//		nsteps (in)			- number of samples to generate
//		x0 (in/out)			- initial state for the random walk
//		x0_out (out)		- final state of the random walk.
//		pOutputSamples (out)- pointer to preallocated array that will contain the output samples, or NULL if we don't want to return samples (for burn-in)
//		nSeparation (in)	- how many samples between every returned sample. Is used to generate less-correlated data.
//		bSequentialBits (in)- true if we want bits to be flipped in a sequential order, false for random order
//	    indices_to_change   - array of indices to be changed by the random walk. Default, if NULL,  is equivalent to an array containing 0..(ncells-1)
//						      which denotes that all the bits are to be changed.
//		num_indices_to_change - number of elements in the array indices_to_change.
//
// Returns:  
//		Generated samples in array pointed by pOutputSamples
void runGibbsSampler(EnergyFunction * pModel, uint32_t nsteps, uint32_t * x0, uint32_t * x0_out, uint32_t * pOutputSamples, uint32_t nSeparation, bool bSequentialBits, uint32_t * indices_to_change = NULL, uint32_t num_indices_to_change=0);

// Returns the marginals of a model (exhaustively computed)
// Input:
//		pModel (in)			- model to generate the samples from
//		pMarginals (out		- preallocated buffer in which the marginals are returned
//
// Returns: 
//		The partition function of the model.
double getMarginals(EnergyFunction * pModel, double * pMarginals);


// Runs the Wang-Landau steps of random walk on the energy function
// nsteps - number of steps
// x - starting point (in/out - returns ending state)
// pModel - model that we use to compute energy
// bin_limits - bin limits for the energy function discretization
// g - energy density function (discretized)
// h - energy function histogram
// separation - how many samples to skip in each MCMC operation
void runWangLandauStep(uint32_t nsteps, uint32_t * x, EnergyFunction *  pModel, uint32_t nbins, double bin_limits[], double g[], double h[], double update_size, uint32_t nSeparation);


// Seeds the random number generator using the system time
void seedRNG();

// Seeds the random number generator with a user-specified seed
// seed - seed for the RNG
void seedRNG(uint32_t seed);