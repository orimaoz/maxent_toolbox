// MATLAB-specific utility functions for the maxent toolbox
// Ori Maoz, August 2016


#include <stdint.h>
#include "mex.h"

// allocates a new buffer and copies the data elementwise into it while converting the data
uint32_t * reallocate_uint32(void * buffer_in, mxClassID sampleClassid, uint32_t nElements);