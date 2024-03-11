/*=================================================================
 *
 * GENERATEHIERACHY.C	.MEX file
 *	        Generate BVH tree for object in 3d space 
 *
 * The calling syntax is:
 *
 *		[connectivity, boundingbox, segment, is_octree_node, level] ...
 *      = generateHierachy(code)
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

/* Output Arguments */

#define	OUT_ARG1	plhs[0]
#define	OUT_ARG2	plhs[1]
#define	OUT_ARG3	plhs[2]
#define	OUT_ARG4	plhs[3]
#define	OUT_ARG5	plhs[4]
#define	OUT_ARG6	plhs[5]

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

#if !defined(MIN)
#define	MIN(A, B)	((A) < (B) ? (A) : (B))
#endif

static Node* _root = NULL;

// Kernals

inline
int clz(int i, int j, int m)
{
    if(m != 0) 
        return __builtin_clz(m);
    else
        return __builtin_clz(m) + __builtin_clz(i ^ j);
}

inline
int taog(int i, int j, uint32_t* restrict sortedMortonCodes, int numObjects)
{
	if(j < 0||j >= numObjects)
		return -1;
    
    return clz(i, j, sortedMortonCodes[i] ^ sortedMortonCodes[j]);
}


void
generateHierarchy(Node* nodes,
                  uint32_t* restrict sortedMortonCodes,
                  //uint32_t*__restrict__ sortedID,
                  //typename vec4<T>::Type *__restrict__ Pos,
                  //typename vec4<T>::Type *__restrict__ boundingBox,
                  int numObjects, //bool verbose
                  int idx)
{
    //LeafNode* leafNodes = new LeafNode[numObjects];
    //InternalNode* internalNodes = new InternalNode[numObjects - 1];
	//////////////////////////////////////////////////////////////////////////////////
	// There are still problems with this kernel.
	// One inclusive scan for morton codes before building the hierachy is believed to be necessary as
	// multiple particles with the same morton code would make the direction determined by
	// (idx - 1, idx) and (idx + 1, idx) unspecified, which may yield particles to be ignored
	// when traversing the tree structure.
	//////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    // Construct leaf nodes.
    // Note: This step can be avoided by storing
    // the tree in a slightly different way.
    //////////////////////////////////////////////////////////////////////////////////

	Node* leafNodes = &nodes[numObjects - 1];
	Node* internalNodes = &nodes[0];

    //int idx = threadIdx.x + blockIdx.x * blockDim.x;

    //if(originalIdx >= N || !compactStream[originalIdx]) return;
	//int idx = compactStream[N + originalIdx];
    //int numObjects = compactStream[2 * N - 1];
        
    if(idx < numObjects)
    {
    	leafNodes[idx].begin = idx;//
    	leafNodes[idx].end = idx; //originalIdx + 1;//
    	leafNodes[idx].isLeafNode = true;
    	leafNodes[idx].isOctreeNode = false;
    }

    //////////////////////////////////////////////////////////////////////////////////
    // Construct internal nodes.
    //////////////////////////////////////////////////////////////////////////////////

    if (idx < numObjects - 1) // in parallel
    {
    	//float4 boundingBox = make_float4(0.0f, 0.0f, 0.0f, 1.0f); // remove later
    	
        //////////////////////////////////////////////////////////////////////////////////
        // Find out which range of objects the node corresponds to.
        // (This is where the magic happens!)
    	//////////////////////////////////////////////////////////////////////////////////

		int d = ( taog(idx, idx + 1, sortedMortonCodes, numObjects) - taog(idx, idx - 1, sortedMortonCodes, numObjects) ) < 0 ? -1 : 1;

		//printf("current internal node tao values : (%d, %d), direction %s\n",  tao(idx, idx - 1, sM, numObjects), tao(idx, idx + 1, sM, numObjects), d > 0 ? "+1" : "-1");
		// Compute lower & upper bound for the length of the range

		int dmin = taog(idx, idx - d, sortedMortonCodes, numObjects);
		int lmax = 2;

		while(taog(idx, idx + lmax * d, sortedMortonCodes, numObjects) > dmin)
		{
			lmax <<= 1;
		};

		// Find the other end using binary search

		int l = 0;
		for(int t = (lmax >> 1); t > 0; t >>= 1)
		{
			if(taog(idx, idx + (l + t) * d, sortedMortonCodes, numObjects) > dmin)
			{
				l += t;
			}
		}

    	// write results.
        
    	int last = idx + l * d;
        int first = idx;

        // swap if necessary
        
        if(first > last)
        {
        	int tmp = last;
        	last = first;
        	first = tmp;
        }

        //////////////////////////////////////////////////////////////////////////////////
        // Determine where to split the range.
        // Identical Morton codes => split the range in the middle.
        //////////////////////////////////////////////////////////////////////////////////

        uint32_t firstCode = sortedMortonCodes[first];
        uint32_t lastCode = sortedMortonCodes[last];

        int split = first;// initial guess

        //if (firstCode == lastCode)
        //	split =  (first + last) >> 1;

        // Calculate the number of highest bits that are the same
        // for all objects, using the count-leading-zeros intrinsic.

        //int commonPrefix = __builtin_clz(firstCode ^ lastCode);
        int commonPrefix = clz(first, last, firstCode ^ lastCode);

        // Use binary search to find where the next bit differs.
        // Specifically, we are looking for the highest object that
        // shares more than commonPrefix bits with the first one.

        int step = last - first;
        do
        {
            step = (step + 1) >> 1; // exponential decrease
            int newSplit = split + step; // proposed new position

            if (newSplit < last)
            {
                uint32_t splitCode = sortedMortonCodes[newSplit];
                int splitPrefix = clz(first, split, firstCode ^ splitCode);
                if (splitPrefix > commonPrefix)
                    split = newSplit; // accept proposal
            }
        }
        while (step > 1);

        //////////////////////////////////////////////////////////////////////////////////
        // Record parent-child relationships.
        //////////////////////////////////////////////////////////////////////////////////

        /*
        // Select childA.

        Node* childA;
        if (split == first)
            childA = &leafNodes[split];
        else
            childA = &internalNodes[split];

        // Select childB.

        Node* childB;
        if (split + 1 == last)
            childB = &leafNodes[split + 1];
        else
            childB = &internalNodes[split + 1];

        // Record parent-child relationships.

        internalNodes[idx].childA = childA;
        internalNodes[idx].childB = childB;
         */
        
        int childA;
        if (split == first)
            childA = split + numObjects - 1;
        else
            childA = split;

        // Select childB.

        int childB;
        if (split + 1 == last)
            childB = split + numObjects;
        else
            childB = split + 1;

        // Record parent-child relationships.

        internalNodes[idx].childA = childA;
        internalNodes[idx].childB = childB;
        
        // Record segementation
        
        uint32_t level = commonPrefix - 2; // last level: 29
        internalNodes[idx].begin = first;
        internalNodes[idx].end = last;
        internalNodes[idx].level = level;
        internalNodes[idx].isLeafNode = false;
        internalNodes[idx].isOctreeNode = false;


        //////////////////////////////////////////////////////////////////////////////////
        // Update bounding box and mass center for octree nodes
        //////////////////////////////////////////////////////////////////////////////////
        
        if(level%3 == 0 && level/3 > 0)
        {
            if(level >= 30){ return; }
        	//if(level == 30)// Last level
        	//{
        	//	internalNodes[idx].isLeafNode = true;
        	//	internalNodes[idx].childA = -1;
        	//	internalNodes[idx].childB = -1;
        	//}

        	internalNodes[idx].isOctreeNode = true;
            
        	//internalNodes[idx].bbox = *boundingBox;
            internalNodes[idx].bbox.x = 0;
            internalNodes[idx].bbox.y = 0;
            internalNodes[idx].bbox.z = 0;
            internalNodes[idx].bbox.w = 1;

        	uint32_t prefix = (firstCode & lastCode) >> (30 - level);
        	for(int i = 0; i < level/3; i++)
        	{
        		float current_scale = internalNodes[idx].bbox.w / (float)(2 << i);
        		internalNodes[idx].bbox.x += ((prefix >> (level - 1 - 3 * i)) & 1) * current_scale ;
        		internalNodes[idx].bbox.y += ((prefix >> (level - 2 - 3 * i)) & 1) * current_scale ;
        		internalNodes[idx].bbox.z += ((prefix >> (level - 3 - 3 * i)) & 1) * current_scale ;
        	}
        	internalNodes[idx].bbox.w /= (float)(1 << (level/3));
        	internalNodes[idx].bbox.x += internalNodes[idx].bbox.w/2 ;
        	internalNodes[idx].bbox.y += internalNodes[idx].bbox.w/2 ;
        	internalNodes[idx].bbox.z += internalNodes[idx].bbox.w/2 ;

        	// Update mass center and total mass
        	
            //internalNodes[idx].scalar.x = 0.0f;
        	//internalNodes[idx].scalar.y = 0.0f;
        	//internalNodes[idx].scalar.z = 0.0f;
        	//internalNodes[idx].scalar.w = 0.0f;

        	/*for(int i = first; i <= last; i++)
        	{
        		float4 p = {posList_x[sortedID[i]],
                posList_y[sortedID[i]],
                posList_z[sortedID[i]],
                posList_w[sortedID[i]]};//childL->begin
        		internalNodes[idx].mass.w += p.w;
        		internalNodes[idx].mass.x += p.x * p.w;
        		internalNodes[idx].mass.y += p.y * p.w;
        		internalNodes[idx].mass.z += p.z * p.w;
        	}

        	if(internalNodes[idx].mass.w) // prevent dividing 0 error
        	{
				internalNodes[idx].mass.x /= internalNodes[idx].mass.w;
				internalNodes[idx].mass.y /= internalNodes[idx].mass.w;
				internalNodes[idx].mass.z /= internalNodes[idx].mass.w;
        	}*/

        }
    }

    // Print tree
    /*if(verbose)
    {
		if (idx < numObjects - 1) // in parallel
		{
			if(nodes[idx].isOctreeNode)
			{
				printf("OctreeNode %d level %d boundingbox (%f, %f, %f, %f) range (%d, %d) mass (%f, %f, %f, %f) constructed.\n",
						idx, nodes[idx].level,
						nodes[idx].bbox.x, nodes[idx].bbox.y, nodes[idx].bbox.z, nodes[idx].bbox.w,
						nodes[idx].begin, nodes[idx].end,
						nodes[idx].mass.x, nodes[idx].mass.y, nodes[idx].mass.z, nodes[idx].mass.w);
			}

		}
		if ( idx < numObjects) // in parallel
		{
			//printf("LeafNode %d range [%d] constructed.\n", idx + numObjects - 1, nodes[idx + numObjects - 1].begin);
		}
    }*/
    // Node 0 is the root.
    //return &internalNodes[0];
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    double *outputArg1,*outputArg2,*outputArg3;
    double *outputArg4, *outputArg5, *outputArg6; 
    double *inputArg1; 
    size_t m1,n1; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 1) { 
	    mexErrMsgIdAndTxt( "Hierachy:invalidNumInputs",
                "One input arguments required."); 
    } else if (nlhs > 5) {
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
    mwSize mem = m1-1;
    OUT_ARG1 = mxCreateDoubleMatrix( mem,2, mxREAL); 
    OUT_ARG2 = mxCreateDoubleMatrix( mem,4, mxREAL); 
    OUT_ARG3 = mxCreateDoubleMatrix( mem,2, mxREAL); 
    OUT_ARG4 = mxCreateDoubleMatrix( mem,1, mxREAL); 
    OUT_ARG5 = mxCreateDoubleMatrix( mem,1, mxREAL);
    OUT_ARG6 = mxCreateDoubleMatrix( mem,4, mxREAL);
    
    /* Assign pointers to the various parameters */ 
    outputArg1 = mxGetPr(OUT_ARG1);
    outputArg2 = mxGetPr(OUT_ARG2);
    outputArg3 = mxGetPr(OUT_ARG3);
    outputArg4 = mxGetPr(OUT_ARG4);
    outputArg5 = mxGetPr(OUT_ARG5);
    outputArg6 = mxGetPr(OUT_ARG6);
    
    inputArg1 = mxGetPr(IN_ARG1); // sorted morton code
    
    /* Do the actual computations in a subroutine */
    
    if(_root)
        free(_root);
    _root = (Node *)calloc(2*m1-1,sizeof(Node));

    for(int i = 0; i < m1; i++)
    {
        generateHierarchy(_root, (uint32_t* restrict)inputArg1, m1, i);
    }

    for(int i = 0; i < mem; i++)
    {        
        outputArg1[i] = (double)_root[i].childA;
        outputArg1[i+mem] = (double)_root[i].childB;
        
        outputArg2[i] = (double)_root[i].bbox.x;
        outputArg2[i+mem] = (double)_root[i].bbox.y;
        outputArg2[i+2*mem] = (double)_root[i].bbox.z;
        outputArg2[i+3*mem] = (double)_root[i].bbox.w;

        outputArg3[i] = (double)_root[i].begin;
        outputArg3[i+mem] = (double)_root[i].end;
        
        outputArg4[i] = (double)_root[i].isOctreeNode;
        outputArg5[i] = (double)_root[i].level;
        
        outputArg6[i] = (double)_root[i].scalar.x;
        outputArg6[i+mem] = (double)_root[i].scalar.y;
        outputArg6[i+2*mem] = (double)_root[i].scalar.z;
        outputArg6[i+3*mem] = (double)_root[i].scalar.w;
        
    }
    
    return;
    
}
