#include "KSyncEnergy.h"

#define MAX_LENGTH 500


// Constructor - constructs it with a set of factors
KSyncEnergy::KSyncEnergy(uint32_t nfactors, double * factors) :
	m_factors(factors, factors + nfactors), m_logz(0)
{
	// There is 1 factor more than there are dimensions (because they are for
	// all the cases from 0 cells firing up to n cells firing)
	m_ndims = nfactors - 1;
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
double KSyncEnergy::getEnergy(std::vector<uint32_t> & x)
{
	uint32_t bitsum = 0;

	// Just count how many ones we have
	for (unsigned int i = 0; i < m_ndims; i++)
	{
		if (x[i]) {
			bitsum++;
		}
	}

	// This the value of x for later so we can make incremental changes
	m_x = x;

	m_nbitson = bitsum;
	m_energy = m_factors[bitsum];

	return m_energy;
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
double KSyncEnergy::propose(uint32_t nbit)
{
	// Save this proposed bit for later
	m_proposed_bit = nbit;

	// Check the bit that changes
	if (m_x[nbit])
	{
		// it is changing 1->0
		m_proposed_energy = m_factors[m_nbitson - 1];
	}
	else
	{
		// it is changing 0->1
		m_proposed_energy = m_factors[m_nbitson + 1];
	}


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
double KSyncEnergy::accept(uint32_t bAccept)
{
	if (bAccept)
	{
		m_x[m_proposed_bit] = !m_x[m_proposed_bit];

		// keep track of how many bits are 1
		if (m_x[m_proposed_bit])
		{
			m_nbitson++;
		}
		else
		{
			m_nbitson--;
		}

		m_energy = m_proposed_energy;
	}
	return m_energy;
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
void KSyncEnergy::accept(double * factor_sum, double p)
{

	// change the state to the proposed state
	accept();

	sumSampleFactor(m_x, factor_sum, p);
}


// Returns the current state of the system
std::vector<uint32_t> * KSyncEnergy::getX()
{
	return &m_x;
}

// Returns the dimensions of this energy functions' inputs
uint32_t KSyncEnergy::getDim()
{
	return m_ndims;
}


// Returns the number of factors (parameters) of the model
uint32_t KSyncEnergy::getNumFactors()
{
	return m_factors.size();
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
void KSyncEnergy::sumSampleFactor(std::vector<uint32_t> & x, double * factor_sum, double p)
{
	uint32_t bitsum = 0;

	// Just count how many ones we have
	for (unsigned int i = 0; i < m_ndims; i++)
	{
		if (x[i]) {
			bitsum++;
		}
	}

	factor_sum[bitsum] += p;

}


// Returns the partition function (it it is known, otherwise just returns 0)
double KSyncEnergy::getLogZ()
{
	return m_logz;
}


// Sets the partition function
void KSyncEnergy::setLogZ(double logz)
{
	m_logz = logz;
}