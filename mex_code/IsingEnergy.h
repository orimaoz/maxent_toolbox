// Class that implements an Ising model energy function
// This an experimental accelerated version which uses a matrix representation of the lagrange multipliers
#include "EnergyFunction.h"
#include <vector>
#pragma once

#include <mkl.h>
#include <ipps.h>
#include <ippvm.h>

#ifdef _WINDOWS
#define DECLARE_ALIGNED _declspec(align(64))
#else
#define DECLARE_ALIGNED
#endif

class IsingEnergy : public EnergyFunction
{


public:



	// Constructor - constructs it from a ("squareform") vector of lagrange multipliers
	IsingEnergy(uint32_t ndims, double * lambdas) :
		m_logz(0)
	{
		m_ndims = ndims;
		m_lambdaMatrix = (double*)ippsMalloc_64f(sizeof(double)*ndims*ndims);;
		m_double_x = (double*)ippsMalloc_64f(sizeof(double)*ndims);
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
	~IsingEnergy()
	{
		// release allocated memory
		ippsFree(m_double_x);
		ippsFree(m_lambdaMatrix);
	}


public:

	// Accepts a vector x and returns its energy. A class implementing this interface is
	// also expected to store x as the current state of the random walk which is used
	// when proposing a new state.
	//
	// Input:
	//		x - state as a vector of boolean entries (0/1)
	//
	// Returns:  
	//		The energy (un-normalized log probability) of the inputed state
	virtual double getEnergy(std::vector<uint32_t> & x);

	// Proposes a new state obtained by a single bit-flip from the current state,
	// and returns the new energy level. This implementation of this function may assume that getEnergy() has been called
	// at some point in the past.
	//
	// Input:
	//		nbit - bit to flip
	//
	// Returns:  
	//		The energy (un-normalized log probability) of the new state after the bit flip
	virtual double propose(uint32_t nbit);

	// Accept/reject the proposed change and updates x.
	// This function may assume that propose() has been called prior to calling this function.
	//
	// Input:
	//		bAccept - if true, then fixes the proposed state as the current state. Otherwise dumps the proposed state and 
	//				  reverts to the pre-proposal state.
	//
	// Returns:  
	//		The energy (un-normalized log probability) of the new state.
	virtual double accept(uint32_t bAccept = 1);

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
	virtual void accept(double * factor_sum, double p);

	// Returns the current state of the system
	virtual std::vector<uint32_t> * getX();

	// Returns the dimensions of this energy functions' inputs
	virtual uint32_t getDim();

	// Returns the number of factors (parameters) of the model
	virtual uint32_t getNumFactors();

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
	virtual void sumSampleFactor(std::vector<uint32_t> & x, double * factor_sum, double p);

	// Returns the partition function (it it is known, otherwise just returns 0)
	virtual double getLogZ();

	// Sets the partition function
	virtual void setLogZ(double logz);


private:

	std::vector<uint32_t> m_x;
	DECLARE_ALIGNED double * m_double_x; // double version of m_x for accelerated computation
	uint32_t m_proposed_bit;
	double m_energy;
	double m_logz;
	double m_proposed_energy[2];
	DECLARE_ALIGNED uint32_t m_ndims;
	double * m_lambdaMatrix;  // lambdas in (symmetric) matrix format
	double * m_lambdaVector;  // lambdas in vector format
};