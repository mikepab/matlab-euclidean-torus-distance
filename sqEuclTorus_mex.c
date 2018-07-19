/*
 * sqEuclTorus_mex.c - MEX implementation of a distance function for a
 * square Euclidean Torus (i.e. square periodic domain). 
 *
 * Inputs:
 * Apos: Nx2 matrix of (x,y) coordinates.
 * Bpos: Mx2 matrix of (x,y) coordinates.
 * L   : Window length along one axis. Assumed square bounds.
 *
 * Outputs:
 * D   : NxM double matrix of distances
 *
 * The calling syntax is:
 * D = sqEuclTorus_mex(Apos,Bpos,L)
 *
 * This is a MEX file for MATLAB.
 * Written by: Mike Pablo (Nov. 1, 2016).
 */

#include "mex.h"
#include "math.h"

// Gatefway fn
// lhs: output args, #, and ptr to array
// rhs: input  args, #, and pts to array
void mexFunction(mwSize nlhs, mxArray *plhs[],
                 mwSize nrhs, const mxArray *prhs[])
{
    // ARGUMENT CHECKS =================================================  
    // Check # of input arguments.
    if (nrhs!=3)
        mexErrMsgIdAndTxt("myMEX:sqEuclTorus:nrhs","Three inputs only.");
    
    if (nlhs!=1)
        mexErrMsgIdAndTxt("myMEX:sqEuclTorus:nlhs","One output only.");

    mwSize ncolsA = mxGetN(prhs[0]);
    mwSize nrowsA = mxGetM(prhs[0]);
    mwSize ncolsB = mxGetN(prhs[1]);
    mwSize nrowsB = mxGetM(prhs[1]);
    
    // First two arguments should be matrices of type double.
    if ( ncolsA != 2 )
        mexErrMsgIdAndTxt("myMEX:sqEuclTorus:nlhs","Matrix must be Nx2");
    if ( ncolsB != 2 )
        mexErrMsgIdAndTxt("myMEX:sqEuclTorus:nlhs","Matrix must be Nx2");
    
    
    // Second two arguments should be scalars.
    if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) 
                              || mxGetNumberOfElements(prhs[2]) != 1 ) {
        mexErrMsgIdAndTxt("myMEX:sqEuclTorus:nlhs",
                          "Domain length must be a scalar, type double.");
    }
    // END ARGUMENT CHECKS ==============================================                
    double *posA = mxGetPr(prhs[0]);
    double *posB = mxGetPr(prhs[1]);
/*        
    // Sanity check
    mexPrintf("A; Cols: %d\nRows: %d\n",ncolsA,nrowsA);
    mexPrintf("B; Cols: %d\nRows: %d\n",ncolsB,nrowsB);
*/  
    // Create memory for output matrix & assign pDists argument to pointer
    plhs[0] = mxCreateDoubleMatrix(nrowsA,nrowsB,mxREAL);
    double *pDists = mxGetPr(plhs[0]);
    
    double L = mxGetScalar(prhs[2]); // Window length.
    int i; 
    int j;
    // Compute square euclidean dists on torus.
    for (i=0;i<nrowsB;i++) {
        for (j=0;j<nrowsA;j++) {
            double absxdiff = fabs(*(posA+j) - *(posB+i));
            absxdiff = fmin(absxdiff,L-absxdiff);
            double absydiff = fabs(*(posA+j+nrowsA) - *(posB+i+nrowsB));
            absydiff = fmin(absydiff,L-absydiff);
                        
            *(pDists+j+i*nrowsA) = absxdiff*absxdiff + absydiff*absydiff;
        }
    }
}