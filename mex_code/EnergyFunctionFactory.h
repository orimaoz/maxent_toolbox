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

	EnergyFunction * createMerpEnergy(const mxArray * model_params);
	EnergyFunction * createFastMerpEnergy(const mxArray * model_params);
	EnergyFunction * createMerpEnergyOld(const mxArray * model_params);
	EnergyFunction * createKSyncEnergy(const mxArray * model_params);
	EnergyFunction * createKIsingEnergy(const mxArray * model_params);
	EnergyFunction * createSparseMerpEnergy(const mxArray * model_params);
	EnergyFunction * createIsingEnergy(const mxArray * model_params);
	EnergyFunction * createPopulationCouplingEnergy(const mxArray * model_params);
	EnergyFunction * createIndependentEnergy(const mxArray * model_params);
	EnergyFunction * createMushEnergy(const mxArray * model_params);
	EnergyFunction * createDendriticMerpEnergy(const mxArray * model_params);		
	
};