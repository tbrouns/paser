Compile with:
mex spc_mex.c SPC.c MergeSort.c L.c FAIRSPLIT.c ALGRAPH.c
 
Change notes:

In "SPC.c":
Added: #include <math.h>
 
In "spc_mex.c":
const mwSize *dims;
dims = mxGetDimensions(prhs[0]); //dimensions of the data array
(other code doesn't work on 64-bit machines)
