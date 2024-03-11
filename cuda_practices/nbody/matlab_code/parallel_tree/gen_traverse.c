/*=================================================================
 *
 * TRAVERSE.C	.MEX file
 *	        Traverse a BVH tree. 
 *
 * The calling syntax is:
 *
 *		    [ax, ay, az] = 
 *              gen_traverse(obj_code,obj_id,x,y,z,w,[c1 c2]);
 *
 *  Author: Y. Liu
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 2020-2021 Y. L.
 *
 *=================================================================*/
#include <math.h>
#include "mex.h"
#include "defines.h"
#include <pthread.h>
/* Input Arguments */

#define	IN_ARG1	prhs[0]
#define	IN_ARG2	prhs[1]
#define	IN_ARG3	prhs[2]
#define	IN_ARG4	prhs[3]
#define	IN_ARG5	prhs[4]
#define	IN_ARG6	prhs[5]
#define	IN_ARG7	prhs[6]

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

static Node* _root = NULL;
static int _root_size = 0;

#ifdef DEBUG
static int max_pos = 0;
#endif

void generateHierarchy(Node* nodes, uint32_t* restrict sortedMortonCodes,
                  int numObjects, int idx);

void computeWeights(Node* internalNodes, uint32_t* restrict sortedID,
		double* restrict posList_x, double* restrict posList_y, 
        double* restrict posList_z, double* restrict posList_w,
        int idx);

double4 
computeAccel(Node* root, double4 my,uint32_t* restrict sortedID,
		double* restrict posList_x, double* restrict posList_y, 
        double* restrict posList_z, double* restrict posList_w,
        double softeningSquared, double graviConstant, const int nObjs);

// Multi-threading
#define NUM_THREAD (4);

typedef struct
{
    int begin,end; // processing nodes
    Node* root;
    uint32_t* obj_morton;
    uint32_t* obj_id;
    int nobjs;
    double *px, *py, *pz, *pw;
    double *soften_const, *gravi_const;
    double *ax, *ay, *az;
    
} inputPackage;

void* generateHierarchyTh(void* inputArg)
{
    inputPackage* data = (inputPackage*)inputArg;
    
    for(int i = data.begin; i < data.end; i++)
    {
        generateHierarchy(
                data.root, data.obj_morton, data.nobjs, i);
    }
    
    pthread_exit(NULL);
};

void* computeWeightsTh(void* inputArg)
{
    inputPackage* data = (inputPackage*)inputArg;
    
    for(int i = data.begin; i < data.end; i++)
    {
        computeWeights(data.root, data.obj_id,
		data.px,data.py,data.pz,data.pw, i);
    }
    
    pthread_exit(NULL);
};

void* computeAccelTh(void* inputArg)
{
    inputPackage data = *(inputPackage*)inputArg;
    
    for(int i = data.begin; i < data.end; i++)
    {
        double4 my = {data.px[i], data.py[i],
                      data.pz[i], data.pw[i]};
        
        double4 ai = computeAccel(data.root, my,
		data.obj_id,
        data.px,data.py,data.pz,data.pw,
        *(data.soften_const), *(data.gravi_const), data.nobjs);
        
        data.ax[i] = (double)ai.x;
        data.ay[i] = (double)ai.y;
        data.az[i] = (double)ai.z;
        
    }
    
    pthread_exit(NULL);
};


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    double *outputArg1,*outputArg2,*outputArg3; 
    double *inputArg1; 
    double *inputArg2,*inputArg3,*inputArg4,*inputArg5; 
    double *inputArg6,*inputArg7;
    size_t m1,n1; 
    
    /* Check for proper number of arguments */
    
    if (nrhs != 7) { 
	    mexErrMsgIdAndTxt( "Hierachy:invalidNumInputs",
                "7 input arguments required."); 
    } else if (nlhs > 3) {
	    mexErrMsgIdAndTxt( "Hierachy:maxlhs",
                "Too many output arguments."); 
    } 
    
    /* Check the dimensions of input arrays. input2 and 3 can be m*1. */ 
    
    m1 = mxGetM(IN_ARG1);n1 = mxGetN(IN_ARG1); 
    
    if ( m1 < NUM_THREAD || n1 != 1 )
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
    
    inputArg1 = mxGetPr(IN_ARG1); // morton
    inputArg2 = mxGetPr(IN_ARG2); // sortedID
    
    inputArg3 = mxGetPr(IN_ARG3); // posList
    inputArg4 = mxGetPr(IN_ARG4);
    inputArg5 = mxGetPr(IN_ARG5);
    inputArg6 = mxGetPr(IN_ARG6);
    
    inputArg7 = mxGetPr(IN_ARG7); // [softenFactor graviFactor]
    
    /* Do the actual computations in a subroutine */
    
    // Gen tree
    if(_root)
    {
        if(m1 != _root_size)
        {
            free(_root);
            _root = (Node *)calloc(2*m1-1,sizeof(Node));
            _root_size = m1;
        }
    }else
    {
        _root = (Node *)calloc(2*m1-1,sizeof(Node));
        _root_size = m1;
    }
    
    // Multi-thread function
    pthread_t threads[NUM_THREADS];
    inputPackage in[NUM_THREAD];
    
    for(int i = 0; i < NUM_THREAD; i++)
    {
        
        // Initialization thread data
        in[i].begin = i*(m1/NUM_THREAD);
        in[i].end   = (i+1)*(m1/NUM_THREAD);
        in[i].end = in[i].end > m1?m1:in[i].end;
        
        in[i].root = _root;
        in[i].obj_morton = (uint32_t*)inputArg1;
        in[i].obj_id = (uint32_t*)inputArg2;
        
        in[i].nobjs = m1;
        
        in[i].px = inputArg3;in[i].py = inputArg4;
        in[i].pz = inputArg5;in[i].pw = inputArg6;

        in[i].soften_const = &inputArg7[0];
        in[i].gravi_const = &inputArg7[1];
        
        in[i].ax = outputArg1;
        in[i].ay = outputArg2;
        in[i].az = outputArg3;
        
        // Create thread
        int ret = pthread_create(&threads[i], NULL, generateHierarchyTh, (void*)&(in[i]));
        if (ret != 0) 
        {
            printf("pthread_create error: error_code = %d\n", ret);
            exit(-1);
        }

    }
    
    // Join
    for(int i = 0; i < NUM_THREAD; i++)
    {
        // Create thread
        int ret = pthread_create(&threads[i], NULL, computeWeightsTh, (void*)&(in[i]));
        if (ret != 0) 
        {
            printf("pthread_create error: error_code = %d\n", ret);
            exit(-1);
        }
    }
    
    for(int i = 0; i < NUM_THREAD; i++)
    {
        // Create thread
        int ret = pthread_create(&threads[i], NULL, computeAccelTh, (void*)&(in[i]));
        if (ret != 0) 
        {
            printf("pthread_create error: error_code = %d\n", ret);
            exit(-1);
        }
    }
    
/*
    for(int i = 0; i < m1; i++)
    {
        generateHierarchy(_root, (uint32_t* restrict)inputArg1, m1, i);
    }
    
    // Compute weights
    for(int i = 0; i < m1; i++)
    {
        computeWeights(_root, (uint32_t* restrict)inputArg2,
		inputArg3,inputArg4,inputArg5,inputArg6, i);
    }
        
    // Traversal
    for(int i = 0; i < m1; i++)
    {
        double4 my = {inputArg3[i], inputArg4[i],
                      inputArg5[i], inputArg6[i]};
        
        double4 ai = computeAccel(_root, my,
		(uint32_t* restrict)inputArg2,
        inputArg3,inputArg4,inputArg5,inputArg6,
        inputArg7[0], inputArg7[1], m1);
        
        outputArg1[i] = (double)ai.x;
        outputArg2[i] = (double)ai.y;
        outputArg3[i] = (double)ai.z;
        
#ifdef DEBUG
    //mexPrintf("Current processing particle = %d\n", i);
#endif
        
    }
    */
#ifdef DEBUG
    mexPrintf("Max stack pos = %d\n", max_pos);
    max_pos = 0;
#endif
    
    return;
    
}

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
                  int numObjects, int idx)
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

        }
    }

    // Print tree
#ifdef DEBUG
    /*{
		if (idx < numObjects - 1) // in parallel
		{
			if(nodes[idx].isOctreeNode)
			{
				mexPrintf("OctreeNode %d level %d boundingbox (%f, %f, %f, %f) range (%d, %d) mass (%f, %f, %f, %f) constructed.\n",
						idx, nodes[idx].level,
						nodes[idx].bbox.x, nodes[idx].bbox.y, nodes[idx].bbox.z, nodes[idx].bbox.w,
						nodes[idx].begin, nodes[idx].end,
						nodes[idx].mass.x, nodes[idx].mass.y, nodes[idx].mass.z, nodes[idx].mass.w);
			}

		}
		if ( idx < numObjects) // in parallel
		{
			mexPrintf("LeafNode %d range [%d] constructed.\n", idx + numObjects - 1, nodes[idx + numObjects - 1].begin);
		}
    }*/
#endif
    // Node 0 is the root.
    //return &internalNodes[0];
}

double4 
computeAccel(Node* root, double4 my,uint32_t* restrict sortedID,
		double* posList_x, double* posList_y, 
        double* posList_z, double* posList_w,
        double softeningSquared, double graviConstant, const int nObjs)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.

    Node* stack[32];
	Node** stackPtr = &stack[0];

	*stackPtr++ = NULL; // push a stop sign

    double4 ai = {0.0f, 0.0f, 0.0f, 0.0f}, r;

    // Traverse nodes starting from the root.
    Node* node = &root[0];
    
#ifdef DEBUG
    int stack_pos = 1;
    int iter = 0;
#endif

    while (node)
    {
#ifdef DEBUG
        if(stack_pos > max_pos) {max_pos = stack_pos;}
#endif
        
        // Check each child node for overlap.

        Node* childL = &root[node->childA];
        Node* childR = &root[node->childB];

    	bool overlapL = true;
    	bool overlapR = true;

    	if(childL->isOctreeNode) // is octree node
    	{
    		float4 bj = childL->scalar;
    		float aabb = childL->bbox.w;
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

    	if(childR->isOctreeNode) // is octree node
    	{
    		float4 bj = childR->scalar;
    		float aabb = childR->bbox.w;
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

        bool sweepL = overlapL && (childL->isLeafNode || (childL->level>29));
        bool sweepR = overlapR && (childR->isLeafNode || (childR->level>29));
        
        if(sweepL)
        {
        	for(int i = childL->begin; i <= childL->end; i++)
        	{
                
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
        
        if(sweepR)
        {
        	for(int i = childR->begin; i <= childR->end; i++)
        	{
                    
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

void computeWeights(Node* internalNodes, uint32_t* restrict sortedID,
		double* posList_x, double* posList_y, 
        double* posList_z, double* posList_w,
        int idx)
{
    if(internalNodes[idx].isOctreeNode)
    {
        uint32_t first = internalNodes[idx].begin;
        uint32_t last = internalNodes[idx].end;
        
        // Update mass center and total mass
        float4 zero = {0.0f,0.0f,0.0f,0.0f};
        internalNodes[idx].scalar = zero;

        for(int i = first; i <= last; i++)
        {
            float4 p = {posList_x[sortedID[i]],
            posList_y[sortedID[i]],
            posList_z[sortedID[i]],
            posList_w[sortedID[i]]};//childL->begin
            internalNodes[idx].scalar.w += p.w;
            internalNodes[idx].scalar.x += p.x * p.w;
            internalNodes[idx].scalar.y += p.y * p.w;
            internalNodes[idx].scalar.z += p.z * p.w;
        }

        if(internalNodes[idx].scalar.w != 0) // prevent dividing 0 error
        {
            internalNodes[idx].scalar.x /= internalNodes[idx].scalar.w;
            internalNodes[idx].scalar.y /= internalNodes[idx].scalar.w;
            internalNodes[idx].scalar.z /= internalNodes[idx].scalar.w;
        }
        #ifdef DEBUG
        //mexPrintf("Current mass center (%f,%f,%f,%f),level %d\n",
        //        internalNodes[idx].scalar.x,internalNodes[idx].scalar.y,
        //        internalNodes[idx].scalar.z,internalNodes[idx].scalar.w,
        //        internalNodes[idx].level/3);
        #endif
    }
}
