// MATLAB-specific utility functions for the maxent toolbox
// Ori Maoz, August 2016

#include "matlab_utils.h"


// allocates a new buffer and copies the data elementwise into it while converting the data
uint32_t * reallocate_uint32(void * buffer_in, mxClassID sampleClassid, uint32_t nElements)
{
	uint32_t * buffer_out;

	// If the input is a different type, we will need to convert it to 32-bit integers
	if ((sampleClassid == mxUINT32_CLASS) || (sampleClassid == mxINT32_CLASS))
	{
		buffer_out = (uint32_t *)buffer_in;
	}
	else if ((sampleClassid == mxCHAR_CLASS) || (sampleClassid == mxLOGICAL_CLASS) || (sampleClassid == mxUINT8_CLASS) || (sampleClassid == mxINT8_CLASS))
	{
		// buffer_out a new array
		buffer_out = new uint32_t[nElements];

		// copy the data to the new array while converting to the new datatype
		char * p_old = (char*)buffer_in;
		for (uint32_t i = 0; i < nElements; i++)
		{
			buffer_out[i] = p_old[i] != 0;
		}
	}
	else if (sampleClassid == mxDOUBLE_CLASS)
	{
		// allocate a new array
		buffer_out = new uint32_t[nElements];

		// copy the data to the new array while converting to the new datatype
		double * p_old = (double*)buffer_in;
		for (uint32_t i = 0; i < nElements; i++)
		{
			buffer_out[i] = p_old[i] != 0;
		}
	}
	else
	{
		mexErrMsgIdAndTxt("maxent_toolbox:init",
			"unsupported type");
	}

	return buffer_out;
}