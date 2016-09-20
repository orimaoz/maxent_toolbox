#include "EnergyFunctionFactory.h"
#include "KSyncEnergy.h"
#include "KIsingEnergy.h"
#include "MerpFastEnergy.h"
#include "MerpSparseEnergy.h"
#include "IsingEnergy.h"
#include "IndependentEnergy.h"
#include <math.h>
#include <cstring>

using namespace std;

#define MAX_STRLEN_MODELTYPE 100


// Reads the model parameters received by matlab and creates an appropriate energy function.
// The returned energy function is dynamically allocate and should be deleted when it is no 
// longer in use.
EnergyFunction * EnergyFunctionFactory::createEnergyFunction(const mxArray * model_params)
{
	// Get model type to see what model we should initialize
	mxArray * mxModelType = mxGetField(model_params,0,"type");
	if (!mxModelType) mexErrMsgIdAndTxt("EnergyFunctionFactory","cannot find model.type");
	char strModelType[MAX_STRLEN_MODELTYPE];
	mxGetString(mxModelType,strModelType,MAX_STRLEN_MODELTYPE);
	if (strcmp(strModelType,"ising") == 0)
	{
		// collect parameters and initialize for Ising model
		return createIsingEnergy(model_params);

	}
	else if (strcmp(strModelType, "merp") == 0)
	{
		// collect parameters and initialize for MERP model
		return createFastMerpEnergy(model_params);
	}
	else if (strcmp(strModelType, "merpsparse") == 0)
	{
		// collect parameters and initialize for MERP model (sparse version)
		return createSparseMerpEnergy(model_params);
	}
	else if (strcmp(strModelType, "ksync") == 0)
	{
		// collect parameters and initialize for ksynchrony model
		return createKSyncEnergy(model_params);
	}
	else if (strcmp(strModelType,"kising") == 0)
	{
		// collect parameters and initialize for kising model
		return createKIsingEnergy(model_params);
	}
	else if (strcmp(strModelType, "indep") == 0)
	{
		// collect parameters and initialize for independent model
		return createIndependentEnergy(model_params);
	}
	else
	{
		// We didn't recognize the model string
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
						  "Unrecognized model inputed in model.type");
		return NULL;
	}
}




// Reads the parameters for an Ising model supplied by matlab and initializes an Ising model
// energy function
EnergyFunction * EnergyFunctionFactory::createIsingEnergy(const mxArray * model_params)
{
	EnergyFunction * pModel;
	mxClassID id;
	// get model factors
	mxArray * mxFactors = mxGetField(model_params, 0, "factors");
	if (!mxFactors) mexErrMsgIdAndTxt("EnergyFunctionFactory", "cannot find model.factors");
	double * pLambdas = (double*)mxGetData(mxFactors);
	size_t nLambdas = mxGetN(mxFactors);

	size_t ndims = (uint32_t)((sqrt(1 + (double)nLambdas * 8) - 1) / 2);

	// Check what type the model factors are
	id = mxGetClassID(mxFactors);
	if (id != mxDOUBLE_CLASS)
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
			"model factors must be of type double");
	}

	// Initialize the Ising model with the parameters we collected
	pModel = new IsingEnergy(ndims, pLambdas);

	// Get the partition function (if it exists)
	mxArray * mxZ = mxGetField(model_params, 0, "z");
	if (mxZ) {
		double z = mxGetScalar(mxZ);
		pModel->setLogZ(z);
	}

	return pModel;
}


EnergyFunction * EnergyFunctionFactory::createKIsingEnergy(const mxArray * model_params)
{
	EnergyFunction * pModel;
	mxClassID id;
	// get model factors
	mxArray * mxFactors = mxGetField(model_params,0,"factors");
	if (!mxFactors) mexErrMsgIdAndTxt("EnergyFunctionFactory","cannot find model.factors");
	double * pLambdas = (double*)mxGetData(mxFactors);
	size_t nLambdas = mxGetN(mxFactors);

	// get number of cells
	mxArray * mxNcells = mxGetField(model_params,0,"ncells");
	if (!mxNcells) mexErrMsgIdAndTxt("EnergyFunctionFactory","cannot find model.ncells");
	uint32_t ncells;
	id = mxGetClassID(mxNcells);
	if (id == mxDOUBLE_CLASS)
	{
		ncells = (uint32_t)mxGetScalar(mxNcells);
	}	
	else if (id == mxUINT32_CLASS)
	{
		uint32_t * pdata = (uint32_t *)mxGetData(mxNcells);
		ncells = *pdata;
	} 
	else
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
							"model.ncells must be of type double or uint32");

	}
	
	// todo - check type of ncells



	// Check what type the model factors are
	id = mxGetClassID(mxFactors);
	if (id != mxDOUBLE_CLASS)
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
							"model factors must be of type double");
	}

	// Initialize the Ising model with the parameters we collected
	pModel = new KIsingEnergy(ncells,pLambdas);

	// Get the partition function (if it exists)
	mxArray * mxZ = mxGetField(model_params,0,"z");
	if (mxZ) {
		double z = mxGetScalar(mxZ);
		pModel->setLogZ(z);
	}

	return pModel;

}


// Reads the parameters for a K-synchrony model supplied by matlab and initializes an K-Synchrony
// model energy function
EnergyFunction * EnergyFunctionFactory::createKSyncEnergy(const mxArray * model_params)
{
	EnergyFunction * pModel;
	mxClassID id;
	// get model factors
	mxArray * mxFactors = mxGetField(model_params,0,"factors");
	if (!mxFactors) mexErrMsgIdAndTxt("EnergyFunctionFactory","cannot find model.factors");
	double * pLambdas = (double*)mxGetData(mxFactors);
	size_t nLambdas = mxGetN(mxFactors);

	// Check what type the model factors are
	id = mxGetClassID(mxFactors);
	if (id != mxDOUBLE_CLASS)
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
							"model factors must be of type double");
	}

	// Initialize the model with the parameters we collected
	pModel = new KSyncEnergy(nLambdas,pLambdas);

	// Get the partition function (if it exists)
	mxArray * mxZ = mxGetField(model_params,0,"z");
	if (mxZ) {
		double z = mxGetScalar(mxZ);
		pModel->setLogZ(z);
	}

	return pModel;
}


// Reads the parameters for a independent model supplied by matlab and initializes an independent
// model energy function
EnergyFunction * EnergyFunctionFactory::createIndependentEnergy(const mxArray * model_params)
{
	EnergyFunction * pModel;
	mxClassID id;

	// get model factors
	mxArray * mxFactors = mxGetField(model_params, 0, "factors");
	if (!mxFactors) mexErrMsgIdAndTxt("EnergyFunctionFactory", "cannot find model.factors");
	double * pLambdas = (double*)mxGetData(mxFactors);
	size_t nLambdas = mxGetN(mxFactors);

	// Check what type the model factors are
	id = mxGetClassID(mxFactors);
	if (id != mxDOUBLE_CLASS)
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
			"model factors must be of type double");
	}

	// get number of cells
	mxArray * mxNcells = mxGetField(model_params, 0, "ncells");
	if (!mxNcells) mexErrMsgIdAndTxt("EnergyFunctionFactory", "cannot find model.ncells");
	uint32_t ncells;
	id = mxGetClassID(mxNcells);
	if (id == mxDOUBLE_CLASS)
	{
		ncells = (uint32_t)mxGetScalar(mxNcells);
	}
	else if (id == mxUINT32_CLASS)
	{
		uint32_t * pdata = (uint32_t *)mxGetData(mxNcells);
		ncells = *pdata;
	}
	else
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
			"model.ncells must be of type double or uint32");

	}

	if (nLambdas != ncells)
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
			"population coupling models must have an equal number of cells and factors");
	}

	// Initialize the Ising model with the parameters we collected
	pModel = new IndependentEnergy(ncells, pLambdas);

	// Get the partition function (if it exists)
	mxArray * mxZ = mxGetField(model_params, 0, "z");
	if (mxZ) {
		double z = mxGetScalar(mxZ);
		pModel->setLogZ(z);
	}

	return pModel;
}



// Initializes using the "fast" bare-bones implementation of MERP
EnergyFunction * EnergyFunctionFactory::createFastMerpEnergy(const mxArray * model_params)
{

	EnergyFunction * pModel;
	mxClassID id;

	// Get number of cells
	mxArray * mxNcells= mxGetField(model_params,0,"ncells");
	if (!mxNcells) mexErrMsgIdAndTxt("EnergyFunctionFactory","cannot find model.ncells");
	uint32_t ncells = mxGetScalar(mxNcells);

	// get model factors
	mxArray * mxFactors = mxGetField(model_params,0,"factors");
	if (!mxFactors) mexErrMsgIdAndTxt("EnergyFunctionFactory","cannot find model.factors");
	double * pLambdas = (double*)mxGetData(mxFactors);
	size_t nfactors = mxGetN(mxFactors);
	id = mxGetClassID(mxFactors);
	if (id != mxDOUBLE_CLASS)
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
							"model factors must be of type double");
	}

	// get the MERP connections and check that sizes match
	mxArray * mxConnections = mxGetField(model_params,0,"connections");
	if (!mxConnections) mexErrMsgIdAndTxt("EnergyFunctionFactory","cannot find model.connections");
	if (mxGetN(mxConnections) != ncells)
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
							"Size of MERP connection matrix size does not match number of cells");
	}
	if (mxGetM(mxConnections) != nfactors)
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
							"Size of MERP connection matrix size does not match number of factors");
	}


	// get the connection matrix
	double * pW = (double*)mxGetPr(mxConnections);

	// get the lagrange factors
	double * pLambda = (double*)mxGetPr(mxFactors);

	// Check if there is a stochastic slope
	double stochastic_slope;
	mxArray * mxSlope = mxGetField(model_params,0,"stochastic_slope");
	if (!mxSlope)
	{
		stochastic_slope = 0;		
	}
	else
	{
		stochastic_slope = mxGetScalar(mxSlope);	
	}


	// Get the threshold
	mxArray * mxThreshold = mxGetField(model_params,0,"threshold");
	if (!mxThreshold) mexErrMsgIdAndTxt("EnergyFunctionFactory","cannot find model.threshold");

	size_t nThresh = mxGetNumberOfElements(mxThreshold);

	if (nThresh == 1)
	{
		// Initialize the MERP model with a single threshold
		double threshold = mxGetScalar(mxThreshold);

		pModel = new MerpFastEnergy(pW,pLambda,ncells,nfactors,threshold,stochastic_slope);

	}
	else if (nThresh == nfactors)
	{

		// Initialize the MERP model with multiple thresholds
		double * pThreshold = (double*)mxGetPr(mxThreshold);

		pModel = new MerpFastEnergy(pW,pLambda,ncells,nfactors,pThreshold,stochastic_slope);

	}
	else
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory","model.threshold must be a scalar or a vector of the same length as factors");
	}



	// Get the partition function (if it exists)
	mxArray * mxZ = mxGetField(model_params,0,"z");
	if (mxZ) {
		double z = mxGetScalar(mxZ);
		pModel->setLogZ(z);
	}

	return pModel;
}




EnergyFunction * EnergyFunctionFactory::createSparseMerpEnergy(const mxArray * model_params)
{

	EnergyFunction * pModel;
	mxClassID id;

	// Get number of cells
	mxArray * mxNcells = mxGetField(model_params, 0, "ncells");
	if (!mxNcells) mexErrMsgIdAndTxt("EnergyFunctionFactory", "cannot find model.ncells");
	uint32_t ncells = mxGetScalar(mxNcells);

	// get model factors
	mxArray * mxFactors = mxGetField(model_params, 0, "factors");
	if (!mxFactors) mexErrMsgIdAndTxt("EnergyFunctionFactory", "cannot find model.factors");
	double * pLambdas = (double*)mxGetData(mxFactors);
	size_t nfactors = mxGetN(mxFactors);
	id = mxGetClassID(mxFactors);
	if (id != mxDOUBLE_CLASS)
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
			"model factors must be of type double");
	}

	// get the MERP connections and check that sizes match
	mxArray * mxConnections = mxGetField(model_params, 0, "connections");
	if (!mxConnections) mexErrMsgIdAndTxt("EnergyFunctionFactory", "cannot find model.connections");
	if (mxGetN(mxConnections) != ncells)
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
			"Size of MERP connection matrix size does not match number of cells");
	}
	if (mxGetM(mxConnections) != nfactors)
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory",
			"Size of MERP connection matrix size does not match number of factors");
	}


	// get the connection matrix
	double * pW = (double*)mxGetPr(mxConnections);

	// get the lagrange factors
	double * pLambda = (double*)mxGetPr(mxFactors);

	// Check if there is a stochastic slope
	double stochastic_slope;
	mxArray * mxSlope = mxGetField(model_params, 0, "stochastic_slope");
	if (!mxSlope)
	{
		stochastic_slope = 0;
	}
	else
	{
		stochastic_slope = mxGetScalar(mxSlope);
	}


	// Get the threshold
	mxArray * mxThreshold = mxGetField(model_params, 0, "threshold");
	if (!mxThreshold) mexErrMsgIdAndTxt("EnergyFunctionFactory", "cannot find model.threshold");

	size_t nThresh = mxGetNumberOfElements(mxThreshold);

	if (nThresh == 1)
	{
		// we got a single threshold, copy it multiple times and initialize with a full vector
		double threshold = mxGetScalar(mxThreshold);
		double * threshold_vector = new double[nfactors];
		for (int i = 0; i < nfactors; i++)
		{
			threshold_vector[i] = threshold;
		}

		pModel = new MerpSparseEnergy(pW, pLambda, ncells, nfactors, threshold_vector, stochastic_slope);

		delete threshold_vector;

	}
	else if (nThresh == nfactors)
	{

		// Initialize the MERP model with multiple thresholds
		double * pThreshold = (double*)mxGetPr(mxThreshold);

		pModel = new MerpSparseEnergy(pW, pLambda, ncells, nfactors, pThreshold, stochastic_slope);

	}
	else
	{
		mexErrMsgIdAndTxt("EnergyFunctionFactory", "model.threshold must be a scalar or a vector of the same length as factors");
	}



	// Get the partition function (if it exists)
	mxArray * mxZ = mxGetField(model_params, 0, "z");
	if (mxZ) {
		double z = mxGetScalar(mxZ);
		pModel->setLogZ(z);
	}

	return pModel;
}
