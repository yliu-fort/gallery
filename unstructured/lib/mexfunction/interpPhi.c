/*=================================================================
 *
 * YPRIME.C	Sample .MEX file corresponding to YPRIME.M
 *	        Solves simple 3 body orbit problem 
 *
 * The calling syntax is:
 *
 *		[yp] = yprime(t, y)
 *
 *  You may also want to look at the corresponding M-code, yprime.m.
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
#define	IN_ARG4	prhs[3]

/* Output Arguments */

#define	OUT_ARG	plhs[0]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

static void interpPhi2d(
		   double	outputArg[],
		   double	(*inputArg1)[2],
 		   double	*inputArg2,
           double	*inputArg3,
           double	(*inputArg4)[2],
           int      nfaces
		   )
{
            
    for(int ifc = 0; ifc < nfaces ; ifc = ifc + 1){
        
        outputArg[ifc] =
                0.5*
                (
                  (
                    *(*(inputArg1 + (int)*(inputArg2 + ifc) - 1))
                    +
                    *(*(inputArg1 + (int)*(inputArg3 + ifc) - 1))
                    
                  )* (*(*(inputArg4 + ifc)))
                +
                  (
                    *(*(inputArg1 + (int)*(inputArg2 + ifc) - 1) + 1)
                    +
                    *(*(inputArg1 + (int)*(inputArg3 + ifc) - 1) + 1)
                    
                  )* (*(*(inputArg4 + ifc) + 1))
                );
        
    }
               
    return;
}

static void interpPhi3d(
		   double	outputArg[],
		   double	(*inputArg1)[3],
 		   double	*inputArg2,
           double	*inputArg3,
           double	(*inputArg4)[3],
           int      nfaces
		   )
{
    
    int owner, neighbour;
            
    for(int ifc = 0; ifc < nfaces ; ifc = ifc + 1){
        
        owner = (int)*(inputArg2 + ifc) - 1;
        neighbour = (int)*(inputArg3 + ifc) - 1;
    
        outputArg[ifc] =
                0.5*
                (
                  (
                    *(*(inputArg1 + owner))
                    +
                    *(*(inputArg1 + neighbour))
                  )* (*(*(inputArg4 + ifc)))
                +
                  (
                    *(*(inputArg1 + owner) + 1)
                    +
                    *(*(inputArg1 + neighbour) + 1)
                  )* (*(*(inputArg4 + ifc) + 1))
                +
                  (
                    *(*(inputArg1 + owner) + 2)
                    +
                    *(*(inputArg1 + neighbour) + 2)
                  )* (*(*(inputArg4 + ifc) + 2))
                );

   }
               
    return;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *outputArg; 
    double *inputArg1, *inputArg2, *inputArg3, *inputArg4; 
    size_t m1,n1,m2,n2,m3,n3,m4,n4; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 4) { 
	    mexErrMsgIdAndTxt( "scalarInterp:invalidNumInputs",
                "Three input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgIdAndTxt( "scalarInterp:maxlhs",
                "Too many output arguments."); 
    } 
    
    /* Check the dimensions of input arrays. input2 and 3 can be m*1. */ 
    
    m1 = mxGetM(IN_ARG1);n1 = mxGetN(IN_ARG1); 
    m2 = mxGetM(IN_ARG2);n2 = mxGetN(IN_ARG2);
    m3 = mxGetM(IN_ARG3);n3 = mxGetN(IN_ARG3);
    m4 = mxGetM(IN_ARG4);n4 = mxGetN(IN_ARG4);
    
    //printf("columns of inputArg3: %d\n", n3);
    
    if ((m2 != m3) || (n2 != 1) || (n3 != 1)||(m1 > 3)||(m1 != m4)||
            (m3 != n4)||(m1 < 2)){
        mexErrMsgIdAndTxt( "scalarInterp:invalidinputArray",
                "Invalid size of inputArray.");
    }
    
    if (!mxIsDouble(IN_ARG1) || mxIsComplex(IN_ARG1) || mxIsSparse(IN_ARG1)||
            !mxIsDouble(IN_ARG2) || mxIsComplex(IN_ARG2) || mxIsSparse(IN_ARG2)||
            !mxIsDouble(IN_ARG3) || mxIsComplex(IN_ARG3) || mxIsSparse(IN_ARG3)||
            !mxIsDouble(IN_ARG4) || mxIsComplex(IN_ARG4) || mxIsSparse(IN_ARG4)
            ){
        mexErrMsgIdAndTxt( "scalarInterp:invalidinputArray",
                "Invalid type of inputArray.");
    }
    
    /* Create a matrix for the return argument */ 
    OUT_ARG = mxCreateDoubleMatrix( (mwSize)m3, (mwSize)n3, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    outputArg = mxGetPr(OUT_ARG);
    
    inputArg1 = mxGetPr(IN_ARG1); 
    inputArg2 = mxGetPr(IN_ARG2); 
    inputArg3 = mxGetPr(IN_ARG3); 
    inputArg4 = mxGetPr(IN_ARG4);
    
    /* Do the actual computations in a subroutine */
    
    if(m1 == 2){
    interpPhi2d(outputArg,inputArg1,inputArg2,inputArg3,inputArg4,m3);
    }
    
    if(m1 == 3){
    interpPhi3d(outputArg,inputArg1,inputArg2,inputArg3,inputArg4,m3);
    }
    
    return;
    
}
