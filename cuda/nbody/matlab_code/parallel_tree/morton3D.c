/*=================================================================
 *
 * MORTON3D.C	.MEX file
 *	        Compute morton code for 3-d coordinates 
 *
 * The calling syntax is:
 *
 *		[code] = morton3D(x, y, z)
 *
 *  Author: Y. Liu
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2011 The MathWorks, Inc.
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	IN_ARG1	prhs[0]
#define	IN_ARG2	prhs[1]
#define	IN_ARG3	prhs[2]

/* Output Arguments */

#define	OUT_ARG1	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

// Expands a 10-bit integer into 30 bits
// by inserting 2 zeros after each bit.
unsigned int expandBits(unsigned int v)
{
    v = (v * 0x00010001u) & 0xFF0000FFu;
    v = (v * 0x00000101u) & 0x0F00F00Fu;
    v = (v * 0x00000011u) & 0xC30C30C3u;
    v = (v * 0x00000005u) & 0x49249249u;
    return v;
}

// Calculates a 30-bit Morton code for the
// given 3D point located within the unit cube [0,1].
unsigned int morton3D(float x, float y, float z)
{
    x = MIN(MAX(x * 1024.0f, 0.0f), 1023.0f);
    y = MIN(MAX(y * 1024.0f, 0.0f), 1023.0f);
    z = MIN(MAX(z * 1024.0f, 0.0f), 1023.0f);
    unsigned int xx = expandBits((unsigned int)x);
    unsigned int yy = expandBits((unsigned int)y);
    unsigned int zz = expandBits((unsigned int)z);
    return xx * 4 + yy * 2 + zz;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    double *outputArg1; 
    double *inputArg1, *inputArg2, *inputArg3; 
    size_t m1,n1,m2,n2,m3,n3; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 3) { 
	    mexErrMsgIdAndTxt( "Morton:invalidNumInputs",
                "Three input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "Morton:maxlhs",
                "Too many output arguments."); 
    } 
    
    /* Check the dimensions of input arrays. input2 and 3 can be m*1. */ 
    
    m1 = mxGetM(IN_ARG1);n1 = mxGetN(IN_ARG1); 
    m2 = mxGetM(IN_ARG2);n2 = mxGetN(IN_ARG2);
    m3 = mxGetM(IN_ARG3);n3 = mxGetN(IN_ARG3);
    //m4 = mxGetM(IN_ARG4);n4 = mxGetN(IN_ARG4);
    
    //printf("columns of inputArg3: %d\n", n3);
    
    if ((m1 != m2) || (m2 != m3) || (n1 != 1) || (n2 != 1) || (n3 != 1))
    {
        mexErrMsgIdAndTxt( "Morton:invalidinputArray",
                "Invalid size of inputArray.");
    }
    
    if (!mxIsDouble(IN_ARG1) || mxIsComplex(IN_ARG1) || mxIsSparse(IN_ARG1)||
        !mxIsDouble(IN_ARG2) || mxIsComplex(IN_ARG2) || mxIsSparse(IN_ARG2)||
        !mxIsDouble(IN_ARG3) || mxIsComplex(IN_ARG3) || mxIsSparse(IN_ARG3))
    {
        mexErrMsgIdAndTxt( "Morton:invalidinputArray",
                "Invalid type of inputArray.");
    }
    
    /* Create a matrix for the return argument */ 
    OUT_ARG1 = mxCreateDoubleMatrix( (mwSize)m1,1, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    outputArg1 = mxGetPr(OUT_ARG1);
    
    inputArg1 = mxGetPr(IN_ARG1); 
    inputArg2 = mxGetPr(IN_ARG2);
    inputArg3 = mxGetPr(IN_ARG3);
    
    /* Do the actual computations in a subroutine */
    for(mwSize i = 0; i < m1; i++)
    {
        outputArg1[i] = (double)morton3D((float)inputArg1[i],(float)inputArg2[i],(float)inputArg3[i]);
    }
    
    return;
    
}