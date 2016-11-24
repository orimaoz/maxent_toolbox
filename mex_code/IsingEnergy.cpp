#include "IsingEnergy.h"

#define MAX_LENGTH 500


// Constructor - constructs it from a  vector of lagrange multipliers
IsingEnergy::IsingEnergy(uint32_t ndims, double * lambdas) :
	m_logz(0)
{
	m_ndims = ndims;
	m_lambdaMatrix = (double*)malloc_aligned(sizeof(double)*ndims*ndims);;
	m_double_x = (double*)malloc_aligned(sizeof(double)*ndims);
	m_lambdaVector = lambdas;

	// fill in the matrix of lambdas - first the single factors (matrix diagonal)...
	double * location = lambdas;
	for (unsigned long i = 0; i < ndims; i++)
	{
		m_lambdaMatrix[i + i*ndims] = *location;
		location++;
	}

	// and now the pairs (both sides of the off-diagonal)
	for (unsigned long i = 0; i < ndims; i++)
	{
		for (unsigned long j = i + 1; j < ndims; j++)
		{
			m_lambdaMatrix[j + i*ndims] = *location;
			m_lambdaMatrix[i + j*ndims] = *location;
			location++;
		}
	}

}


// destructor
IsingEnergy::~IsingEnergy()
{
	// release allocated memory
	free_aligned(m_double_x);
	free_aligned(m_lambdaMatrix);
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
double IsingEnergy::getEnergy(uint32_t *  x)
{
	double logprob = 0;

	// First we will find which values of x are non-zero so we can parse it more efficiently
	// (this way we avoid going over all pair combinations and just go over the ones with value of 1)
	unsigned int ones_vec[MAX_LENGTH];
	unsigned int total_ones = 0;
	for (unsigned int i=0; i < m_ndims; i++)
	{
		if (x[i]) {
			// remember this nonzero value
			ones_vec[total_ones++] = i;

			// while we are it, sum up its probability
			logprob += m_lambdaMatrix[i + i*m_ndims];
		}		

		m_double_x[i] = x[i];
	}

	// Now we go over the pairs
	for (unsigned int i=0; i < total_ones; i++)
	{
		for (unsigned int j=i+1; j < total_ones; j++)
		{
			unsigned int first = ones_vec[i];
			unsigned int second = ones_vec[j];

			logprob += m_lambdaMatrix[first + second*m_ndims];

		}
	}



	// This the value of x for later so we can make incremental changes
	m_x.assign(x, x + m_ndims);

	m_energy = logprob;
	m_proposed_energy[0] = m_energy;
	
	return m_energy;
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
void IsingEnergy::sumSampleFactor(uint32_t *  x, double* factor_sum, double p)
{

	// First we will find which values of x are non-zero so we can parse it more efficiently
	// (this way we avoid going over all pair combinations and just go over the ones with value of 1)
	unsigned int ones_vec[MAX_LENGTH];
	unsigned int total_ones = 0;
	for (unsigned int i=0; i < m_ndims; i++)
	{
		if (x[i]) {
			// remember this nonzero value
			ones_vec[total_ones++] = i;

			// while we are it, sum it up
			factor_sum[i]+=p;
		}			
	}

	// Now we go over the pairs
	for (unsigned int i=0; i < total_ones; i++)
	{
		for (unsigned int j=i+1; j < total_ones; j++)
		{
			unsigned int first = ones_vec[i];
			unsigned int second = ones_vec[j];

			// Compute the offset in the lambda vector
			unsigned int idx = (first+1)*(m_ndims + (m_ndims - first) ) / 2 + second - first - 1;
			factor_sum[idx]+=p;

		}
	}
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
double IsingEnergy::propose(uint32_t nbit)
{
	double dEnergy=0;   // energy difference from original sample

	m_proposed_bit = nbit;

	// go over only the relevant row that was changed
	double * factor_row = m_lambdaMatrix + nbit * m_ndims;

	ippsDotProd_64f(factor_row, m_double_x, m_ndims, &dEnergy);


	// The single factors (diagonal element) will always be changed, so add it if it wasn't added during the previous loop
	dEnergy += factor_row[nbit] * !m_x[nbit];


	// Increase the energy if we change from 1 to 0 and decrease it elsewise.
	// this is equivalent to an "if" statement but hopefully a bit faster.
	m_proposed_energy[1]= m_energy - dEnergy * ((((double)m_x[nbit]) * 2) - 1);

	return m_proposed_energy[1];
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
double IsingEnergy::accept(uint32_t bAccept)
{
	// switch the bit only if bAccept is true
	m_x[m_proposed_bit] = m_x[m_proposed_bit] ^ bAccept;
	m_double_x[m_proposed_bit] = (bool)m_double_x[m_proposed_bit] ^ bAccept;
	m_energy = m_proposed_energy[bAccept];
	m_proposed_energy[0] = m_proposed_energy[bAccept];

	return m_energy;
}

// accepts the proposed changes and adds current factor state to a running sum (marginal)
void IsingEnergy::accept(double * factor_sum, double prob)
{
	// change the state to the proposed state
	accept();

	sumSampleFactor(m_x.data(),factor_sum,prob);
}


// Returns the current state of the system
uint32_t *  IsingEnergy::getX()
{
	return m_x.data();
}



// Returns the dimensions of this energy functions' inputs
uint32_t IsingEnergy::getDim()
{
	return m_ndims;
}

// Returns the number of factors (parameters) of the model
uint32_t IsingEnergy::getNumFactors()
{
	return m_ndims * (m_ndims + 1) / 2;

}


// Returns the partition function (it it is known, otherwise just returns 0)
double IsingEnergy::getLogZ()
{
	return m_logz;
}

// Sets the partition function
void IsingEnergy::setLogZ(double logz)
{
	m_logz = logz;
}
