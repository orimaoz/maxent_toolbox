// Creates and initializes classes that implement energy functions of Boltzmann distributions.
// Ori Maoz 07/2014

#include "mex.h"
#include "EnergyFunction.h"


class EnergyFunctionFactory
{

public:
	// Reads the model parameters received by matlab and creates an appropriate energy function.
	// The returned energy function is dynamically allocate dan dshould be deleted when it is no 
	// longer in use.
	EnergyFunction * createEnergyFunction(const mxArray * model_params);

private:
	EnergyFunction * createMerpEnergy(const mxArray * model_params, bool is_sparse);
	EnergyFunction * createKSyncEnergy(const mxArray * model_params);
	EnergyFunction * createKIsingEnergy(const mxArray * model_params);
	EnergyFunction * createIsingEnergy(const mxArray * model_params);
	EnergyFunction * createIndependentEnergy(const mxArray * model_params);
	EnergyFunction * createHighOrderEnergy(const mxArray * model_params);

	// helper functions for extracting model params
	uint32_t getNcells(const mxArray * model_params);  // return number of cells
	double getZ(const mxArray * model_params);	   // return partition function
	double * getFactors(const mxArray * model_params);  // return factors
	uint32_t getNFactors(const mxArray * model_params); // return number of factors

};