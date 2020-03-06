/*=================================================================
 *
 * TRAVERSE.C	.MEX file
 *	        Traverse a BVH tree. 
 *
 * The calling syntax is:
 *
 *		    [ax, ay, az] = traverse(obj_id,x,y,z,w,...
 *          connectivity,w_stat, boundingbox, is_octree_node, ...
 *          c1, c2);
 *
 *  Author: Y. Liu
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2011 The MathWorks, Inc.
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"

#include "defines.h"
/* Input Arguments */

#define	IN_ARG1	prhs[0]
#define	IN_ARG2	prhs[1]
#define	IN_ARG3	prhs[2]
#define	IN_ARG4	prhs[3]
#define	IN_ARG5	prhs[4]
#define	IN_ARG6	prhs[5]
#define	IN_ARG7	prhs[6]
#define	IN_ARG8	prhs[7]
#define	IN_ARG9	prhs[8]
#define	IN_ARG10	prhs[9]
#define	IN_ARG11	prhs[10]
#define	IN_ARG12	prhs[11]
#define	IN_ARG13	prhs[12]

/* Output Arguments */

#define	OUT_ARG1	plhs[0]
#define	OUT_ARG2	plhs[1]
#define	OUT_ARG3	plhs[2]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

#ifdef DEBUG
static int max_pos = 0;
#endif

static int stack[32];

double4 
computeAccel(int idx, double4 my,uint32_t* restrict sortedID,
		double *posList_x, double *posList_y, double *posList_z, double *posList_w,
		double* childA, double* childB, 
        double *mass_x, double *mass_y, double *mass_z, double *mass_w,
        double *bbox_w,
        double *isOctreeNode, double *segmentA, double *segmentB, double *level,
        double softeningSquared, double graviConstant, const int nObjs)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.
	//int stack[32];
	int* stackPtr = &stack[0];

	*stackPtr++ = -1; // push a stop sign

    double4 ai = {0.0f, 0.0f, 0.0f, 0.0f}, r;

    // Traverse nodes starting from the root.
    int node = 0;
    
#ifdef DEBUG
    int stack_pos = 1;
    int iter = 0;
#endif

    while (node != -1)
    {
#ifdef DEBUG
        if(stack_pos > max_pos) {max_pos = stack_pos;}
#endif
        
        // Check each child node for overlap.

        int childL = (int)childA[node];
        int childR = (int)childB[node];

    	bool overlapL = true;
    	bool overlapR = true;

    	if(childL < nObjs-1 && isOctreeNode[childL] > 0) // is octree node
    	{
    		double4 bj = {mass_x[childL],mass_y[childL],mass_z[childL],mass_w[childL]};
    		double aabb = bbox_w[childL];
    		r.x = bj.x - my.x;
    		r.y = bj.y - my.y;
    		r.z = bj.z - my.z;
	        double distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
	        distSqr += softeningSquared;

	        double invDist = 1.0f/sqrt(distSqr);

    		overlapL = (aabb * invDist) > 0.5; // default parameter for barnes-hut algorithm.
    		if(!overlapL)
    		{
    			double invDistCube = invDist * invDist * invDist;
                double s = graviConstant * bj.w * invDistCube; // mass * invdistcube

    			if(my.w < 0)
                {
                    s = -s;// mass * invdistcube
                    ai.x += r.x * s;
                    ai.y += r.y * s;
                    ai.z += r.z * s;
                }
                else
                {
                    ai.x += r.x * s;
                    ai.y += r.y * s;
                    ai.z += r.z * s;
                }
    		}
    	}

    	if(childR < nObjs-1 && isOctreeNode[childR] > 0) // is octree node
    	{
    		double4 bj = {mass_x[childR],mass_y[childR],mass_z[childR],mass_w[childR]};
    		double aabb = bbox_w[childR];
    		r.x = bj.x - my.x;
    		r.y = bj.y - my.y;
    		r.z = bj.z - my.z;
	        double distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
	        distSqr += softeningSquared;

	        double invDist = 1.0f/sqrt(distSqr);

    		overlapR = (aabb * invDist) > 0.5; // default parameter for barnes-hut algorithm.
    		if(!overlapR)
    		{
    	        double invDistCube =  invDist * invDist * invDist;
                double s = graviConstant * bj.w * invDistCube; // mass * invdistcube

    			if(my.w < 0)
                {
                    s = -s;// mass * invdistcube
                    ai.x += r.x * s;
                    ai.y += r.y * s;
                    ai.z += r.z * s;
                }
                else
                {
                    ai.x += r.x * s;
                    ai.y += r.y * s;
                    ai.z += r.z * s;
                }
    		}
    	}

        bool sweepL = false;
        bool sweepR = false;
        int begin=0, end=0;
        
        if (overlapL)
        {
            if(childL > nObjs-2) // point to leaf
            {
                sweepL = true;
                begin = childL - nObjs + 1;
                end = childL - nObjs + 1;
            }
            else 
            {
                if(level[childL] > 29) 
                {
                    sweepL = true;
                    begin = (int)segmentA[childL];
                    end = (int)segmentB[childL];
                }
            }
        }
        
        if(sweepL)
        {
        	for(int i = begin; i <= end; i++)
        	{
                if(idx == sortedID[i]) {continue;}
                
				double4 sh = {posList_x[sortedID[i]],
                posList_y[sortedID[i]],
                posList_z[sortedID[i]],
                posList_w[sortedID[i]]};//childL->begin

				r.x = sh.x - my.x;
				r.y = sh.y - my.y;
				r.z = sh.z - my.z;
				double distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
				distSqr += softeningSquared;

				double invDist = 1.0f/sqrt(distSqr);
				double invDistCube =  invDist * invDist * invDist;
                double s = graviConstant * sh.w * invDistCube; // mass * invdistcube

    			if(my.w < 0)
                {
                    s = -s;// mass * invdistcube
                    ai.x += r.x * s;
                    ai.y += r.y * s;
                    ai.z += r.z * s;
                }
                else
                {
                    ai.x += r.x * s;
                    ai.y += r.y * s;
                    ai.z += r.z * s;
                }

        	}
        }

        if (overlapR)
        {
            if(childR > nObjs-2)
            {
                sweepR = true;
                begin = childR - nObjs + 1;
                end = childR - nObjs + 1;
            }
            else 
            {
                if(level[childR] > 29)
                {
                    sweepR = true;
                    begin = segmentA[childR];
                    end = segmentB[childR];
                }
            }
        }
        
        if(sweepR)
        {
        	for(int i = begin; i <= end; i++)
        	{
                if(idx == sortedID[i]) {continue;}
                    
				double4 sh = {posList_x[sortedID[i]],
                posList_y[sortedID[i]],
                posList_z[sortedID[i]],
                posList_w[sortedID[i]]};//childL->begin
                
				r.x = sh.x - my.x;
				r.y = sh.y - my.y;
				r.z = sh.z - my.z;
				double distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
				distSqr += softeningSquared;

				double invDist = 1.0f/sqrt(distSqr);
				double invDistCube =  invDist * invDist * invDist;
                double s = graviConstant * sh.w * invDistCube; // mass * invdistcube

    			if(my.w < 0)
                {
                    s = -s;// mass * invdistcube
                    ai.x += r.x * s;
                    ai.y += r.y * s;
                    ai.z += r.z * s;
                }
                else
                {
                    ai.x += r.x * s;
                    ai.y += r.y * s;
                    ai.z += r.z * s;
                }
        	}
        }

        // Query overlaps an internal node => traverse.
        bool traverseL = (overlapL && !sweepL);
        bool traverseR = (overlapR && !sweepR);

 #ifdef DEBUG
                iter++;
                if(iter > 2*nObjs - 1) 
                    mexPrintf("Current node %d\n",node);
                if(iter > 2*nObjs + 9) 
                {
                    //for(int p = 0; p < 32; p++)
                    //    mexPrintf("Stack pos %d = %d\n",p,stack[p]);
                    //mexPrintf("Current node %d\n",node);
                    mexPrintf("***\n");
                    break;
                }
#endif    
        
        if (!traverseL && !traverseR)
        {
        	node = *--stackPtr; // pop
#ifdef DEBUG
            stack_pos--;
#endif
        }
        else
        {
            node = (traverseL) ? childL : childR;
            if (traverseL && traverseR)
            {
            	*stackPtr++ = childR; // push
#ifdef DEBUG
                stack_pos++;
#endif
            }
        }
        

        
    }

    return ai;
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    double *outputArg1,*outputArg2,*outputArg3; 
    double *inputArg1; 
    double *inputArg2,*inputArg3,*inputArg4,*inputArg5; 
    double *inputArg6,*inputArg7,*inputArg8,*inputArg9; 
    double *inputArg10,*inputArg11, *inputArg12, *inputArg13; 
    size_t m1,n1; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 13) { 
	    mexErrMsgIdAndTxt( "Hierachy:invalidNumInputs",
                "12 input arguments required."); 
    } else if (nlhs > 3) {
	    mexErrMsgIdAndTxt( "Hierachy:maxlhs",
                "Too many output arguments."); 
    } 
    
    /* Check the dimensions of input arrays. input2 and 3 can be m*1. */ 
    
    m1 = mxGetM(IN_ARG1);n1 = mxGetN(IN_ARG1); 
    
    if ( n1 != 1 )
    {
        mexErrMsgIdAndTxt( "Hierachy:invalidinputArray",
                "Invalid size of inputArray.");
    }
    
    if (!mxIsUint32(IN_ARG1) || mxIsComplex(IN_ARG1) || mxIsSparse(IN_ARG1))
    {
        mexErrMsgIdAndTxt( "Hierachy:invalidinputArray",
                "Invalid type of inputArray.");
    }
    
    /* Create a matrix for the return argument */ 
    OUT_ARG1 = mxCreateDoubleMatrix( m1,1, mxREAL); 
    OUT_ARG2 = mxCreateDoubleMatrix( m1,1, mxREAL); 
    OUT_ARG3 = mxCreateDoubleMatrix( m1,1, mxREAL); 
    
    /* Assign pointers to the various parameters */ 
    outputArg1 = mxGetPr(OUT_ARG1);
    outputArg2 = mxGetPr(OUT_ARG2);
    outputArg3 = mxGetPr(OUT_ARG3);
    
    inputArg1 = mxGetPr(IN_ARG1); // sortedID
    
    inputArg2 = mxGetPr(IN_ARG2); // posList
    inputArg3 = mxGetPr(IN_ARG3);
    inputArg4 = mxGetPr(IN_ARG4);
    inputArg5 = mxGetPr(IN_ARG5);
    
    inputArg6 = mxGetPr(IN_ARG6); // child
    
    inputArg7 = mxGetPr(IN_ARG7); // tree-mass
    inputArg8 = mxGetPr(IN_ARG8); // tree-aabb
    inputArg9 = mxGetPr(IN_ARG9); // tree-isOctreeNode
    inputArg10 = mxGetPr(IN_ARG10); // tree-segment
    inputArg11 = mxGetPr(IN_ARG11); // tree-level
    
    inputArg12 = mxGetPr(IN_ARG12); // softenFactor
    inputArg13 = mxGetPr(IN_ARG13); // graviFactor
    
    /* Do the actual computations in a subroutine */
    for(int i = 0; i < m1; i++)
    {
        double4 my = {inputArg2[i], inputArg3[i], inputArg4[i], inputArg5[i]};
        
        double4 ai = computeAccel(i,my,
		(uint32_t*)inputArg1,
        inputArg2,inputArg3,inputArg4,inputArg5,
		&inputArg6[0],&inputArg6[m1-1],
        &inputArg7[0], &inputArg7[m1-1], &inputArg7[2*(m1-1)], &inputArg7[3*(m1-1)],
        &inputArg8[3*(m1-1)], inputArg9, &inputArg10[0], &inputArg10[m1-1], inputArg11,
        inputArg12[0], inputArg13[0], m1);
        
        outputArg1[i] = (double)ai.x;
        outputArg2[i] = (double)ai.y;
        outputArg3[i] = (double)ai.z;
        
#ifdef DEBUG
    //mexPrintf("Current processing particle = %d\n", i);
#endif
        
    }
    
#ifdef DEBUG
    mexPrintf("Max stack pos = %d\n", max_pos);
    max_pos = 0;
#endif
    
    return;
    
}
