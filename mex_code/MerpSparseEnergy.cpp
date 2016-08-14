#include "MerpSparseEnergy.h"
#include <cstring>

#ifdef _WINDOWS
#define DECLARE_ALIGNED _declspec(align(64))
#else
#define DECLARE_ALIGNED
#endif


// Constructor - constructs it from a random projection and a single threshold
MerpSparseEnergy::MerpSparseEnergy(double* in_W, double * in_lambda, uint32_t ncells, uint32_t nfactors, double * in_thresholds, double bStochasticSlope)
	: m_ndims(ncells), m_nfactors(nfactors), m_logz(0), m_bProposed(false)
{


	m_ndims = ncells;
	m_nfactors = nfactors;

	// Check whether our model has stochastic firing
	m_stochasticSlope = bStochasticSlope;
	m_bStochastic = (m_stochasticSlope != 0);


	// allocate required arrays
	m_W_val = (double**)ippsMalloc_64f(sizeof(double*)*m_ndims);	// values of each column of W
	m_W_idx = (MKL_INT**)ippsMalloc_64f(sizeof(MKL_INT*)*m_ndims);	// indices for each column of W
	m_W_nz = (uint32_t*)ippsMalloc_64f(sizeof(uint32_t)*m_ndims);	// # of elements for each column of W

	m_sparse_lambda = (double**)ippsMalloc_64f(sizeof(double*)*m_ndims);	// values of lambda for each column of W

	m_lambda = (double*)ippsMalloc_64f(sizeof(double)*m_nfactors);
	m_threshold = (double*)ippsMalloc_64f(sizeof(double)*m_nfactors);
	m_y = (double*)ippsMalloc_64f(sizeof(double)*m_nfactors);
	m_tmp = (double*)ippsMalloc_64f(sizeof(double)*m_nfactors);


	// copy data to the arrays
	std::memcpy(m_lambda, in_lambda, sizeof(double) * m_nfactors);
	std::memset(m_y, 0, sizeof(double) * m_nfactors);

	// The matrix W is stored in a sparse form. We will allocate and initialize each of its row separately.
	double * ptrW = in_W;
	for (uint32_t idxCol = 0; idxCol < ncells; idxCol++)
	{
		// Allocate a long enough column so that we can hold even matrix with all elements nonzero.
		// This doesn't take so much memory so we don't really need to conserve (which would mean to pre-count)
		m_W_val[idxCol] = (double*)ippsMalloc_64f(sizeof(double)*m_nfactors);	// values of current column
		m_W_idx[idxCol] = (MKL_INT*)ippsMalloc_64f(sizeof(MKL_INT)*m_nfactors);	// indices of current column

		double * p_sparse_W_val = m_W_val[idxCol];
		MKL_INT * p_sparse_W_idx = m_W_idx[idxCol];

		// Count the nonzero elements 
		unsigned int nz = 0;
		for (MKL_INT idxRow = 0; idxRow < m_nfactors; idxRow++)
		{
			// Copy the nonzero element to the sparse array
			if (ptrW[idxRow])
			{
				m_W_val[idxCol][nz] = ptrW[idxRow];
				m_W_idx[idxCol][nz] = idxRow;
				nz++;

			}
		}
		// save the number of elements
		m_W_nz[idxCol] = nz;

		// advance to the next row
		ptrW += m_nfactors;
	}

	for (uint32_t idxLambda = 0; idxLambda < ncells; idxLambda++)
	{
		// allocate memory for this sparse bit of lambda and gather the relevant entries into the array
		m_sparse_lambda[idxLambda] = (double*)ippsMalloc_64f(sizeof(double)*m_nfactors);
		cblas_dgthr(m_W_nz[idxLambda], m_lambda, m_sparse_lambda[idxLambda], m_W_idx[idxLambda]);
	}


	// copy the thresholds
	std::memcpy(m_threshold, in_thresholds, sizeof(double) * m_nfactors);
}


// destructor
MerpSparseEnergy::~MerpSparseEnergy()
{
	// free preallocated memory

	// clear sparse lambdas - first all the columns and then the base
	for (uint32_t idxLambda = 0; idxLambda < m_ndims; idxLambda++)
	{
		// allocate memory for this sparse bit of lambda and gather the relevant entries into the array
		ippsFree(m_sparse_lambda[idxLambda]);
	}
	ippsFree(m_sparse_lambda);

	// Free the sparse matrix W - first all the columns and then the base
	for (uint32_t idxCol = 0; idxCol < m_ndims; idxCol++)
	{
		ippsFree(m_W_val[idxCol]);
		ippsFree(m_W_idx[idxCol]);
	}

	ippsFree(m_W_val);
	ippsFree(m_W_idx);
	ippsFree(m_W_nz);

	ippsFree(m_lambda);
	ippsFree(m_threshold);
	ippsFree(m_y);
	ippsFree(m_tmp);
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
double MerpSparseEnergy::getEnergy(uint32_t * x)
{
	m_x.assign(x, x + m_ndims);

	// Implement random projection - we sum up the columns of W according to x

	std::memset(m_y,0, sizeof(double) * m_nfactors);
#pragma vector aligned
	for (uint32_t i = 0; i < m_ndims; i++)
	{
		// Check for each column if we are to sum it
		if (x[i])
		{
			// sparse addition into m_y
			cblas_daxpyi(m_W_nz[i], 1, m_W_val[i], m_W_idx[i], m_y);
		}
	}

	// subtract the tresholds so we can use check if we are above or below 0
	ippsSub_64f_I(m_threshold, m_y, m_nfactors);

	m_energy = applyThreshold(m_y);


	return m_energy;
}


// Applies threshold on the projection and returns the summed results
double MerpSparseEnergy::applyThreshold(double * y)
{
	DECLARE_ALIGNED double outprod = 0;
	DECLARE_ALIGNED unsigned int i;

	#pragma vector aligned 
	for (i = 0; i < m_nfactors; i++)
	{
		// deterministic - fire only if we passed the threshold
		if (y[i] > 0)
		{
			outprod += m_lambda[i];
		}
	}

	return outprod;
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
void MerpSparseEnergy::sumSampleFactor(uint32_t * x, double* factor_sum,double p)
{
	// Implement random projection - we sum up the columns of W according to x

	std::memset(m_y, 0, sizeof(double) * m_nfactors);

	for (uint32_t i = 0; i < m_ndims; i++)
	{
		// Check for each column if we are to sum it
		if (x[i])
		{
			// sparse addition into m_y
			cblas_daxpyi(m_W_nz[i], 1, m_W_val[i], m_W_idx[i], m_y);

		}
	}

	// subtract the tresholds so we can use check if we are above or below 0
	ippsSub_64f_I(m_threshold, m_y, m_nfactors);

	// set spiking/nonspiking according to the cutoff threshold
	#pragma vector aligned 
	for (unsigned int i = 0; i < m_nfactors; i++)
	{
		if (m_y[i] > 0)
		{
			factor_sum[i] += p;
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
double MerpSparseEnergy::propose(uint32_t nbit)
{
	m_proposed_bit = nbit;
	if (m_x[nbit])
	{
		m_proposed_energy = m_energy + getSumDiffBitOff(nbit);
	}
	else
	{
		m_proposed_energy = m_energy + getSumDiffBitOn(nbit);
	}

	return m_proposed_energy;
}



double MerpSparseEnergy::getSumDiffBitOn(uint32_t nbit)
{
	MKL_INT * idx = m_W_idx[nbit];
	double * val = m_W_val[nbit];
	double energydiff=0;

	// gather the relevant m_y entries into a compressed aray
	cblas_dgthr(m_W_nz[nbit], m_y, m_tmp, m_W_idx[nbit]);
	double * sparse_lambda = m_sparse_lambda[nbit];

	// now we can have a tight loop on the gathered data
#pragma vector aligned
	for (uint32_t i = 0; i < m_W_nz[nbit]; i++)
	{
		// sum the factor contribution before the bit change
		if (m_tmp[i] > 0)
		{
            energydiff -= sparse_lambda[i];
		}
	}

	// now add the relevant elements of m_y and see what passed the threshold
	ippsAdd_64f_I(val, m_tmp, m_W_nz[nbit]);

#pragma vector aligned
	for (uint32_t i = 0; i < m_W_nz[nbit]; i++)
	{
		if (m_tmp[i] > 0)
		{
			energydiff += sparse_lambda[i];
		}

	}

	return energydiff;
}




double MerpSparseEnergy::getSumDiffBitOff(uint32_t nbit)
{
	MKL_INT * idx = m_W_idx[nbit];
	double * val = m_W_val[nbit]; 
	double energydiff = 0;

	// gather the relevant m_y entries into a compressed aray
	cblas_dgthr(m_W_nz[nbit], m_y, m_tmp, m_W_idx[nbit]);
	double * sparse_lambda = m_sparse_lambda[nbit];

	// now we can have a tight loop on the gathered data
#pragma vector aligned
	for (uint32_t i = 0; i < m_W_nz[nbit]; i++)
	{
		// sum the factor contribution before the bit change
		if (m_tmp[i] > 0)
		{
            energydiff -= sparse_lambda[i];
		}
	}

	// now subtract the relevant elements of m_y and see what passed the threshold
	ippsSub_64f_I(val, m_tmp, m_W_nz[nbit]);

#pragma vector aligned
	for (uint32_t i = 0; i < m_W_nz[nbit]; i++)
	{
		if (m_tmp[i] > 0)
		{
			energydiff += sparse_lambda[i];
		}

	}

	return energydiff;
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
double MerpSparseEnergy::accept(uint32_t bAccept)
{
	if (bAccept)
	{
		// Change m_y to reflect the updated sums before we flip the bit
		if (m_x[m_proposed_bit])
		{
			// bit changing 1->0
			cblas_daxpyi(m_W_nz[m_proposed_bit], -1, m_W_val[m_proposed_bit], m_W_idx[m_proposed_bit], m_y);
	
		}
		else
		{
			// bit changing 0->1
			// sparse addition into m_y
			cblas_daxpyi(m_W_nz[m_proposed_bit], 1, m_W_val[m_proposed_bit], m_W_idx[m_proposed_bit], m_y);
		}

		// now flip the proposed bit
		m_x[m_proposed_bit] = !m_x[m_proposed_bit];

		m_energy = m_proposed_energy;
		m_bProposed  = false;
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
void MerpSparseEnergy::accept(double * factor_sum, double prob)
{
	// change the state to the proposed state
	accept();

	// after accepting, m_y contains the sub-threshold activity so we just need to sum up
	// the activity that passes the threshold
	#pragma vector aligned 
	for (size_t i = 0; i < m_nfactors; i++)
	{
		if (m_y[i]>0)
		{
			factor_sum[i] += prob;
		}
	}
}




// Returns the current state of the system
uint32_t * MerpSparseEnergy::getX()
{
	return m_x.data();
}


// Returns the dimensions of this energy functions' inputs
uint32_t MerpSparseEnergy::getDim()
{
	return m_ndims;
}

// Returns the number of factors (parameters) of the model
uint32_t MerpSparseEnergy::getNumFactors()
{
	return m_nfactors;
}


// Returns the partition function (it it is known, otherwise just returns 0)
double MerpSparseEnergy::getLogZ()
{
	return m_logz;
}


// Sets the partition function
void MerpSparseEnergy::setLogZ(double logz)
{
	m_logz = logz;
}
