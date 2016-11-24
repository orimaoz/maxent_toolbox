#include "EnergyFunctionFactory.h"
#include "KSyncEnergy.h"
#include "KIsingEnergy.h"
#include "IsingEnergy.h"
#include "IndependentEnergy.h"
#include "HighOrderEnergy.h"
#include "CompositeEnergy.h"
#include <math.h>
#include <cstring>

using namespace std;

#define MAX_STRLEN_MODELTYPE 100


// Reads the model parameters received by matlab and creates an appropriate energy function.
// The returned energy function is dynamically allocate and should be deleted when it is no 
// longer in use.
EnergyFunction * EnergyFunctionFactory::createEnergyFunction(const mxArray * model_params)
{
	EnergyFunction * outEnergy = NULL;

	// Get model type to see what model we should initialize
	mxArray * mxModelType = mxGetField(model_params,0,"type");
	if (!mxModelType) mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory","cannot find model.type");

	// check what type it is
	mxClassID typeClass = mxGetClassID(mxModelType);
	if (typeClass == mxCHAR_CLASS)
	{
		// use the string to check what class of model it is
		char strModelType[MAX_STRLEN_MODELTYPE];
		mxGetString(mxModelType, strModelType, MAX_STRLEN_MODELTYPE);
		if (strcmp(strModelType, "ising") == 0)
		{
			// collect parameters and initialize for Ising model
			outEnergy = createIsingEnergy(model_params);

		}
		else if (strcmp(strModelType, "ksync") == 0)
		{
			// collect parameters and initialize for ksynchrony model
			outEnergy = createKSyncEnergy(model_params);
		}
		else if (strcmp(strModelType, "kising") == 0)
		{
			// collect parameters and initialize for kising model
			outEnergy = createKIsingEnergy(model_params);
		}
		else if (strcmp(strModelType, "indep") == 0)
		{
			// collect parameters and initialize for independent model
			outEnergy = createIndependentEnergy(model_params);
		}
		else if (strcmp(strModelType, "highorder") == 0)
		{
			// collect parameters and initialize for high-order model
			outEnergy = createHighOrderEnergy(model_params);
		}
		else if (strcmp(strModelType, "composite") == 0)
		{
			// it's a composite model made up of several smaller models
			
			// It should contain a list of inner models, get that list
			mxArray * mxInnerModels = mxGetField(model_params, 0, "innermodels");
			if (!mxInnerModels) mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory", "cannot find model.innerModels for composite model");

			// check what type it is
			mxClassID typeClass = mxGetClassID(mxInnerModels);

			if (typeClass != mxCELL_CLASS) mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory", "model.innerModels must be a cell array of models");

			// type is a cell array of composite models - recursively call on all the sub classes and build a composite class
			size_t num_nested_models = mxGetNumberOfElements(mxInnerModels);

			// check that the nested models cell array is not empty
			if (num_nested_models == 0)
			{
				mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory",
					"nested model cell array is empty");
			}

			// create a nested model
			uint32_t ncells = getNcells(model_params);
			CompositeEnergy * pCompositeEnergy = new CompositeEnergy(ncells);

			// go through all of them and create a list of nested models
			for (int i = 0; i < num_nested_models; i++)
			{
				// get a pointer to the nested model params
				mxArray * mxNestedModel = mxGetCell(mxInnerModels, i);

				// recursively create an energy function for it
				EnergyFunction * pNestedModel = createEnergyFunction(mxNestedModel);

				// add it to the composite model
				pCompositeEnergy->addEnergyFunction(pNestedModel);
			}

			// set partition function for the composite energy
			pCompositeEnergy->setLogZ(getZ(model_params));
			outEnergy = pCompositeEnergy;

		}
		else
		{
			// We didn't recognize the model string
			mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory",
				"Unrecognized model string in model.type");
			outEnergy = NULL;
		}

	}
	else
	{
		mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory", "model.type must be a string or a cell array of models");
	}



	// We should never actually reach this point
	return outEnergy;
}




// Reads the parameters for an Ising model supplied by matlab and initializes an Ising model
// energy function
EnergyFunction * EnergyFunctionFactory::createIsingEnergy(const mxArray * model_params)
{
	EnergyFunction * pModel;
	mxClassID id;

	// get model factors
	double * pLambdas = getFactors(model_params);
	size_t nLambdas = getNFactors(model_params);

	size_t ndims = getNcells(model_params);

	// Initialize the Ising model with the parameters we collected
	pModel = new IsingEnergy(ndims, pLambdas);

	// Get the partition function (if it exists)
	pModel->setLogZ(getZ(model_params));

	return pModel;
}


EnergyFunction * EnergyFunctionFactory::createKIsingEnergy(const mxArray * model_params)
{
	EnergyFunction * pModel;
	mxClassID id;

	// get model factors
	double * pLambdas = getFactors(model_params);
	size_t nLambdas = getNFactors(model_params);

	// get number of cells
	uint32_t ncells = getNcells(model_params);

	pModel = new KIsingEnergy(ncells,pLambdas);

	// Get the partition function (if it exists)
	pModel->setLogZ(getZ(model_params));

	return pModel;

}


// Reads the parameters for a K-synchrony model supplied by matlab and initializes an K-Synchrony
// model energy function
EnergyFunction * EnergyFunctionFactory::createKSyncEnergy(const mxArray * model_params)
{
	EnergyFunction * pModel;
	mxClassID id;

	// get model factors
	double * pLambdas = getFactors(model_params);
	size_t nLambdas = getNFactors(model_params);

	// Initialize the model with the parameters we collected
	pModel = new KSyncEnergy(nLambdas,pLambdas);

	// Get the partition function (if it exists)
	pModel->setLogZ(getZ(model_params));

	return pModel;
}


// Reads the parameters for a independent model supplied by matlab and initializes an independent
// model energy function
EnergyFunction * EnergyFunctionFactory::createIndependentEnergy(const mxArray * model_params)
{
	EnergyFunction * pModel;
	mxClassID id;

	// get model factors
	double * pLambdas = getFactors(model_params);
	size_t nLambdas = getNFactors(model_params);

	// get number of cells
	uint32_t ncells = getNcells(model_params);

	if (nLambdas != ncells)
	{
		mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory",
			"population coupling models must have an equal number of cells and factors");
	}

	// Initialize the Ising model with the parameters we collected
	pModel = new IndependentEnergy(ncells, pLambdas);

	// Get the partition function (if it exists)
	pModel->setLogZ(getZ(model_params));

	return pModel;
}






// Initializes a model using an arbitrary list of high order correlations
EnergyFunction * EnergyFunctionFactory::createHighOrderEnergy(const mxArray * model_params)
{
	EnergyFunction * pModel;
	mxClassID id;

	// get number of cells
	uint32_t ncells = getNcells(model_params);


	// get model factors
	double * pLambda = getFactors(model_params);
	size_t nfactors = getNFactors(model_params);

	// get the factors and check that sizes match
	mxArray * mxConnections = mxGetField(model_params, 0, "factorMatrix");
	if (!mxConnections) mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory", "cannot find model.factorMatrix");
	if (mxGetN(mxConnections) != ncells)
	{
		mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory",
			"Size of high order factor matrix does not match number of cells");
	}
	if (mxGetM(mxConnections) != nfactors)
	{
		mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory",
			"Size of high order factor matrix does not match number of factors");
	}


	// get the connection matrix
	float * pW = (float*)mxGetPr(mxConnections);


	// Get the threshold
	mxArray * mxThreshold = mxGetField(model_params, 0, "threshold");
	if (!mxThreshold) mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory", "cannot find model.threshold");

	size_t nThresh = mxGetNumberOfElements(mxThreshold);

	if (nThresh == nfactors)
	{
		float * pThreshold = (float*)mxGetPr(mxThreshold);

		pModel = new HighOrderEnergy(pW, pLambda, ncells, nfactors, pThreshold);

	}
	else
	{
		mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory", "model.threshold must be a vector of the same length as factors");
	}



	// Get the partition function (if it exists)
	pModel->setLogZ(getZ(model_params));



	return pModel;
}



// returns the number of cells in a model params structure
uint32_t EnergyFunctionFactory::getNcells(const mxArray * model_params)
{

	mxArray * mxNcells = mxGetField(model_params, 0, "ncells");
	if (!mxNcells) mexErrMsgIdAndTxt("EnergyFunctionFactory", "cannot find model.ncells");
	uint32_t ncells;
	mxClassID id = mxGetClassID(mxNcells);
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
		mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory",
			"model.ncells must be of type double or uint32");

	}

	return ncells;
}


// returns a pointer to the model factors in a models param structure
double * EnergyFunctionFactory::getFactors(const mxArray * model_params)
{
	// get model factors
	mxArray * mxFactors = mxGetField(model_params, 0, "factors");
	if (!mxFactors) mexErrMsgIdAndTxt("EnergyFunctionFactory", "cannot find model.factors");
	double * pLambdas = (double*)mxGetData(mxFactors);
	mxClassID id = mxGetClassID(mxFactors);
	if (id != mxDOUBLE_CLASS)
	{
		mexErrMsgIdAndTxt("maxent:EnergyFunctionFactory",
			"model factors must be of type double");
	}
	return pLambdas;
}



// returns number of factors in a models param structure
uint32_t EnergyFunctionFactory::getNFactors(const mxArray * model_params)
{
	// get model factors
	mxArray * mxFactors = mxGetField(model_params, 0, "factors");
	size_t nfactors = mxGetN(mxFactors);
	return nfactors;
}



// returns the partition function model params structure
double EnergyFunctionFactory::getZ(const mxArray * model_params)	   // return partition function
{
	double z;
	mxArray * mxZ = mxGetField(model_params, 0, "z");
	if (mxZ) {
		z = mxGetScalar(mxZ);
	}
	else
	{
		z = 0;
	}
	return z;
}
