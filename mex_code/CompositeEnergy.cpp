// composite energy function
#include "CompositeEnergy.h"



// Constructor - constructs an empty composite
CompositeEnergy::CompositeEnergy(uint32_t ncells) :
	m_logz(0)
{
	m_ndims = ncells;
}


// destructor - removes all internally kept energy function (the destructor deletes them)
CompositeEnergy::~CompositeEnergy()
{
	for (int i=0; i < m_internalEnergies.size(); i++)
	{
		delete m_internalEnergies[i];
	}
}

void CompositeEnergy::addEnergyFunction(EnergyFunction * energy)
{
	// store the energy
	m_internalEnergies.push_back(energy);

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
double CompositeEnergy::getEnergy(uint32_t * x)
{
	// the energy is just the sum of all the internal energies
	double energy_sum = 0;
	for (int i=0; i < m_internalEnergies.size(); i++)
	{
		energy_sum += m_internalEnergies[i]->getEnergy(x);
	}

	return energy_sum;
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
double CompositeEnergy::propose(uint32_t nbit)
{
	// the energy is just the sum of all the internal energies
	m_proposed_energy = 0;
	for (int i=0; i < m_internalEnergies.size(); i++)
	{
		m_proposed_energy += m_internalEnergies[i]->propose(nbit);
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
double CompositeEnergy::accept(uint32_t bAccept)
{
	// delegate to the children
	double energy_sum = 0;
	for (int i=0; i < m_internalEnergies.size(); i++)
	{
		energy_sum += m_internalEnergies[i]->accept(bAccept);
	}

	return energy_sum;

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
void CompositeEnergy::accept(double * factor_sum, double prob)
{
	// delegate to the children, keep track of the location in the factor sum vector
	uint32_t factor_sum_location = 0;
	for (int i=0; i < m_internalEnergies.size(); i++)
	{
		m_internalEnergies[i]->accept(factor_sum + factor_sum_location,prob);
		factor_sum_location += m_internalEnergies[i]->getNumFactors();
	}

}


// Returns the current state of the system
uint32_t * CompositeEnergy::getX()
{
	return m_internalEnergies[0]->getX();
}

// Returns the dimensions of this energy functions' inputs
uint32_t CompositeEnergy::getDim()
{
	return m_ndims;
}

// Returns the number of factors (parameters) of the model
uint32_t CompositeEnergy::getNumFactors()
{
	uint32_t factor_sum = 0;
	for (int i=0; i < m_internalEnergies.size(); i++)
	{
		factor_sum += m_internalEnergies[i]->getNumFactors();
	}
	return factor_sum;
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
void CompositeEnergy::sumSampleFactor(uint32_t * x, double * factor_sum, double p)
{
	// delegate to the children, keep track of the location in the factor sum vector
	uint32_t factor_sum_location = 0;
	for (int i=0; i < m_internalEnergies.size(); i++)
	{
		m_internalEnergies[i]->sumSampleFactor(x,factor_sum + factor_sum_location,p);
		factor_sum_location += m_internalEnergies[i]->getNumFactors();
	}

}

// Returns the partition function (it it is known, otherwise just returns 0)
double CompositeEnergy::getLogZ()
{
	return m_logz;
}

// Sets the partition function
void CompositeEnergy::setLogZ(double logz)
{
	m_logz = logz;
}

