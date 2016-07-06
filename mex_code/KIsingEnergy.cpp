#include "KIsingEnergy.h"



// Constructor - constructs it with a set of factors
KIsingEnergy::KIsingEnergy(uint32_t ncells, uint32_t nfactors, double * factors) :
	m_KSyncEnergy(ncells + 1, factors),
	m_IsingEnergy(ncells, factors + (ncells + 1)),
	m_logz(0)
{
	// The first (n+1) factors are k-sync factors, the rest are ising			

	m_ndims = ncells;
	m_nfactors = nfactors;
}


// Accepts a vector x and returns its energy. A class implementing this interface is
// also expected to store x as the current state of the random walk which is used
// when proposing a new state.
//
// Input:
//		x - state as a vector of boolean entries (0/1)
//
// Returns:  
//		The energy (un-normalized log probability) of the inputed state
double KIsingEnergy::getEnergy(std::vector<uint32_t> & x)
{
	return m_KSyncEnergy.getEnergy(x) + m_IsingEnergy.getEnergy(x);
}

// Proposes a new state obtained by a single bit-flip from the current state,
// and returns the new energy level. This implementation of this function may assume that getEnergy() has been called
// at some point in the past.
//
// Input:
//		nbit - bit to flip
//
// Returns:  
//		The energy (un-normalized log probability) of the new state after the bit flip
double KIsingEnergy::propose(uint32_t nbit)
{
	m_proposed_energy = m_KSyncEnergy.propose(nbit) + m_IsingEnergy.propose(nbit);
	return m_proposed_energy;
}

// Accept/reject the proposed change and updates x.
// This function may assume that propose() has been called prior to calling this function.
//
// Input:
//		bAccept - if true, then fixes the proposed state as the current state. Otherwise dumps the proposed state and 
//				  reverts to the pre-proposal state.
//
// Returns:  
//		The energy (un-normalized log probability) of the new state.
double KIsingEnergy::accept(uint32_t bAccept)
{
	return m_KSyncEnergy.accept(bAccept) + m_IsingEnergy.accept(bAccept);
}

// Accepts the proposed changes and adds factors of current state to a running sum (marginal).
// This function may assume that propose() has been called prior to calling this function.
//
// Input:
//		factor_sum - vector of marginals (one for each factor). The function should sum the factors of the proposed state
//					 into this vector after multipying them by the argument 'p'.
//		p	-		 multiply factors by this number before summing them into factor_sum.
//
// Returns:  
//		(none)
void KIsingEnergy::accept(double * factor_sum, double prob)
{
	m_KSyncEnergy.accept(factor_sum, prob);
	m_IsingEnergy.accept(factor_sum + (m_ndims + 1), prob);

}


// Returns the current state of the system
std::vector<uint32_t> * KIsingEnergy::getX()
{
	return m_IsingEnergy.getX();
}

// Returns the dimensions of this energy functions' inputs
uint32_t KIsingEnergy::getDim()
{
	return m_ndims;
}

// Returns the number of factors (parameters) of the model
uint32_t KIsingEnergy::getNumFactors()
{
	return m_nfactors;
}

// Adds the factors of the inputed sample to a vector of factor sums, multiplied by the inputed probabilities. This is mainly used to compute marginals.
//
// Input:
//		x	-		 state to compute the factors of.
//		factor_sum - vector of marginals (one for each factor). The function should sum the factors of the inputed state
//					 into this vector after multipying them by the argument 'p'.
//		p	-		 multiply factors by this number before summing them into factor_sum.
//
// Returns:  
//		(none)
void KIsingEnergy::sumSampleFactor(std::vector<uint32_t> & x, double * factor_sum, double p)
{
	// Each of the energy functions only sees part of the factors, so we transfer
	// them corresponding halves of the factor sum to each of the sub-functions
	m_KSyncEnergy.sumSampleFactor(x, factor_sum, p);
	m_IsingEnergy.sumSampleFactor(x, factor_sum + (m_ndims + 1), p);
}

// Returns the partition function (it it is known, otherwise just returns 0)
double KIsingEnergy::getLogZ()
{
	return m_logz;
}

// Sets the partition function
void KIsingEnergy::setLogZ(double logz)
{
	m_logz = logz;
}

