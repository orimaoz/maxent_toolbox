#include "IndependentEnergy.h"
#include <cstring>


// constructor
IndependentEnergy::IndependentEnergy(uint32_t ndims, double * factors) :
m_logz(0)
{
	m_ndims = ndims;
	m_factors = factors;
	m_x = new uint32_t[ndims]; // allocate memory for current state
}

// destructor
IndependentEnergy::~IndependentEnergy()
{
	// remove memory allocated for current state
	delete m_x;
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
double IndependentEnergy::getEnergy(uint32_t * x)
{
	double energy = 0;
	unsigned long active_bits = 0;

	std::memcpy(m_x, x, m_ndims * sizeof(uint32_t));

	for (unsigned long i = 0; i < m_ndims; i++)
	{
		// add active lagrange multipliers
		energy += (double)x[i] * m_factors[i];
	}


	m_energy = energy;
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
double IndependentEnergy::propose(uint32_t nbit)
{
	m_proposed_bit = nbit;
	if (m_x[nbit])
	{
		m_proposed_energy = m_energy - m_factors[nbit];
	}
	else
	{
		m_proposed_energy = m_energy + m_factors[nbit];
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
double IndependentEnergy::accept(uint32_t bAccept)
{
	if (bAccept) {
		m_x[m_proposed_bit] = !m_x[m_proposed_bit];
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
void IndependentEnergy::accept(double * factor_sum, double p)
{
	// change the state to the proposed state
	accept();

	sumSampleFactor(m_x, factor_sum, p);
}


// Returns the current state of the system
uint32_t * IndependentEnergy::getX()
{
	return m_x;
}


// Returns the dimensions of this energy functions' inputs
uint32_t IndependentEnergy::getDim()
{
	return m_ndims;
}

// Returns the number of factors (parameters) of the model
uint32_t IndependentEnergy::getNumFactors()
{
	// number off factors in this model is the same as the number of dimensions
	return m_ndims;
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
void IndependentEnergy::sumSampleFactor(uint32_t * x, double * factor_sum, double p)
{
	double bitsum = 0;

	// activate the active bits multiplied by the population activity
	for (unsigned int i=0; i < m_ndims; i++)
	{
		factor_sum[i] += p * (double)x[i];
	}

}
// Returns the partition function (it it is known, otherwise just returns 0)
double IndependentEnergy::getLogZ()
{
	return m_logz;
}

// Sets the partition function
void IndependentEnergy::setLogZ(double logz)
{
	m_logz = logz;
}