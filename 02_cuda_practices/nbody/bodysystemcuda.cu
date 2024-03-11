/*
 * Copyright 1993-2015 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */
/*
 * Copyright 2019 Yuxuan Liu.  All rights reserved.
 *
 * Implemented parallel tree construction based barnes-hut algorithm.
 *
 */
#include <helper_cuda.h>
#include <math.h>
#include "timer.h"

#if defined(__APPLE__) || defined(MACOSX)
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

// CUDA standard includes
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

#include <cooperative_groups.h>

namespace cg = cooperative_groups;

#include "bodysystem.h"

#include <thrust/device_ptr.h>
#include <thrust/sort.h>

__constant__ float softeningSquared;
__constant__ double softeningSquared_fp64;

__device__ unsigned int retirementCount = 0; // for function::reduceMax
cudaError_t setRetirementCount(int retCnt) // for function::reduceMax
{
    return cudaMemcpyToSymbol(retirementCount, &retCnt, sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
}

cudaError_t setSofteningSquared(float softeningSq)
{
	return cudaMemcpyToSymbol(softeningSquared,
			&softeningSq,
			sizeof(float), 0,
			cudaMemcpyHostToDevice);
}

cudaError_t setSofteningSquared(double softeningSq)
{
	return cudaMemcpyToSymbol(softeningSquared_fp64,
			&softeningSq,
			sizeof(double), 0,
			cudaMemcpyHostToDevice);
}

template<class T>
struct SharedMemory
{
	__device__ inline operator       T *()
    		{
		extern __shared__ int __smem[];
		return (T *)__smem;
    		}

	__device__ inline operator const T *() const
    		{
		extern __shared__ int __smem[];
		return (T *)__smem;
    		}
};

template<typename T>
__device__ T rsqrt_T(T x)
{
	return rsqrt(x);
}

template<>
__device__ float rsqrt_T<float>(float x)
{
	return rsqrtf(x);
}

template<>
__device__ double rsqrt_T<double>(double x)
{
	return rsqrt(x);
}


// Macros to simplify shared memory addressing
#define SX(i) sharedPos[i+blockDim.x*threadIdx.y]
                        // This macro is only used when multithreadBodies is true (below)
#define SX_SUM(i,j) sharedPos[i+blockDim.x*j]

template <typename T>
__device__ T getSofteningSquared()
{
	return softeningSquared;
}
template <>
__device__ double getSofteningSquared<double>()
{
	return softeningSquared_fp64;
}

template <typename T>
struct Node
{
public:
	__device__
	Node(): childA(NULL), childB(NULL), begin(0), end(0) {};

//protected:
	__device__ ~Node(){};
	typename vec4<T>::Type bbox; // store min coordinates and scale for bounding box
	typename vec4<T>::Type mass;
	int begin, end;
	unsigned int level;

	Node<T>* childA;
	Node<T>* childB;
	bool isOctreeNode;
	bool isLeafNode;

};

template <typename T>
struct DeviceData
{
    T *dPos[2]; // mapped host pointers
    T *dVel;
    T *dAcc;
    T *dInfo; // contain current center and scale
    Node<T>* root;
    unsigned int *morton_code;
    unsigned int *sortedID;
    //unsigned int *compactStream;
    cudaEvent_t  event;
    unsigned int offset;
    unsigned int numBodies;
};


template <typename T>
__device__ typename vec3<T>::Type
bodyBodyInteraction(typename vec3<T>::Type ai,
		typename vec4<T>::Type bi,
		typename vec4<T>::Type bj)
{
	typename vec3<T>::Type r;

	// r_ij  [3 FLOPS]
	          r.x = bj.x - bi.x;
	          r.y = bj.y - bi.y;
	          r.z = bj.z - bi.z;

	          // distSqr = dot(r_ij, r_ij) + EPS^2  [6 FLOPS]
	          T distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
	          distSqr += getSofteningSquared<T>();

	          // invDistCube =1/distSqr^(3/2)  [4 FLOPS (2 mul, 1 sqrt, 1 inv)]
	          T invDist = rsqrt_T(distSqr);
	          T invDistCube =  invDist * invDist * invDist;

	          // s = m_j * invDistCube [1 FLOP]
	          T s = bj.w * invDistCube;

	          // a_i =  a_i + s * r_ij [6 FLOPS]
	          ai.x += r.x * s;
	          ai.y += r.y * s;
	          ai.z += r.z * s;

	          return ai;
}

template <typename T>
__device__ typename vec3<T>::Type
computeBodyAccel(typename vec4<T>::Type bodyPos,
		typename vec4<T>::Type *positions,
		int numTiles, cg::thread_block cta)
{
	typename vec4<T>::Type *sharedPos = SharedMemory<typename vec4<T>::Type>();

	typename vec3<T>::Type acc = {0.0f, 0.0f, 0.0f};

	for (int tile = 0; tile < numTiles; tile++)
	{
		sharedPos[threadIdx.x] = positions[tile * blockDim.x + threadIdx.x];

		cg::sync(cta);

		// This is the "tile_calculation" from the GPUG3 article.
#pragma unroll 128

		for (unsigned int counter = 0; counter < blockDim.x; counter++)
		{
			acc = bodyBodyInteraction<T>(acc, bodyPos, sharedPos[counter]);
		}

		cg::sync(cta);
	}

	return acc;
}

template<typename T>
__global__ void
integrateBodies(typename vec4<T>::Type *__restrict__ newPos,
		typename vec4<T>::Type *__restrict__ oldPos,
		typename vec4<T>::Type *vel,
		unsigned int deviceOffset, unsigned int deviceNumBodies,
		float deltaTime, float damping, int numTiles)
{
	// Handle to thread block group
	cg::thread_block cta = cg::this_thread_block();
	int index = blockIdx.x * blockDim.x + threadIdx.x;

	if (index >= deviceNumBodies)
	{
		return;
	}

	typename vec4<T>::Type position = oldPos[deviceOffset + index];

	typename vec3<T>::Type accel = computeBodyAccel<T>(position,
			oldPos,
			numTiles, cta);

	// acceleration = force / mass;
	// new velocity = old velocity + acceleration * deltaTime
			// note we factor out the body's mass from the equation, here and in bodyBodyInteraction
	// (because they cancel out).  Thus here force == acceleration
	typename vec4<T>::Type velocity = vel[deviceOffset + index];

	velocity.x += accel.x * deltaTime;
	velocity.y += accel.y * deltaTime;
	velocity.z += accel.z * deltaTime;

	velocity.x *= damping;
	velocity.y *= damping;
	velocity.z *= damping;

	// new position = old position + velocity * deltaTime
	position.x += velocity.x * deltaTime;
	position.y += velocity.y * deltaTime;
	position.z += velocity.z * deltaTime;

	// store new position and velocity
	newPos[deviceOffset + index] = position;
	vel[deviceOffset + index]    = velocity;
}

template<typename T, int NUM_THREAD_PER_BLOCK>
__global__ void
reduceMaxSinglePass(typename vec4<T>::Type *__restrict__ out,
					typename vec4<T>::Type *__restrict__ g_idata,
					typename vec3<T>::Type *__restrict__ g_odata, unsigned int n)
{
	cg::thread_block cta = cg::this_thread_block();
	typename vec3<T>::Type __shared__ smem1[NUM_THREAD_PER_BLOCK];
	typename vec3<T>::Type __shared__ smem2[NUM_THREAD_PER_BLOCK];
	// set thread id
	unsigned int tid = threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// boundary check
	//if (idx >= n) return;
	// convert global data pointer to the local pointer of this block
	//int* idata = g_idata + blockIdx.x * blockDim.x * 8;
	if(idx < n)
	{
		smem1[tid].x = g_idata[idx].x;
		smem1[tid].y = g_idata[idx].y;
		smem1[tid].z = g_idata[idx].z;

		smem2[tid].x = g_idata[idx].x;
		smem2[tid].y = g_idata[idx].y;
		smem2[tid].z = g_idata[idx].z;
	}
	__syncthreads();
	// unrolling 8 data blocks
	typename vec3<T>::Type tmpMax = smem1[0];
	typename vec3<T>::Type tmpMin = smem2[0];

#pragma unroll 8
	for(unsigned int i = 0; i < 8; i++)
	{
		if (idx + i * blockDim.x < n)
		{
			typename vec4<T>::Type a = g_idata[idx + i * blockDim.x];

			tmpMax.x = max(tmpMax.x, a.x);
			tmpMax.y = max(tmpMax.y, a.y);
			tmpMax.z = max(tmpMax.z, a.z);

			tmpMin.x = min(tmpMin.x, a.x);
			tmpMin.y = min(tmpMin.y, a.y);
			tmpMin.z = min(tmpMin.z, a.z);
		}
	}
	smem1[tid] = tmpMax;
	smem2[tid] = tmpMin;

	__syncthreads();

	// in-place reduction  and complete unroll
	if (blockDim.x >= 1024 && tid < 512)
	{
		smem1[tid].x = max(smem1[tid].x, smem1[tid + 512].x);
		smem1[tid].y = max(smem1[tid].y, smem1[tid + 512].y);
		smem1[tid].z = max(smem1[tid].z, smem1[tid + 512].z);

		smem2[tid].x = min(smem2[tid].x, smem2[tid + 512].x);
		smem2[tid].y = min(smem2[tid].y, smem2[tid + 512].y);
		smem2[tid].z = min(smem2[tid].z, smem2[tid + 512].z);
	}
	__syncthreads();
	if (blockDim.x >= 512 && tid < 256)
	{
		smem1[tid].x = max(smem1[tid].x, smem1[tid + 256].x);
		smem1[tid].y = max(smem1[tid].y, smem1[tid + 256].y);
		smem1[tid].z = max(smem1[tid].z, smem1[tid + 256].z);

		smem2[tid].x = min(smem2[tid].x, smem2[tid + 256].x);
		smem2[tid].y = min(smem2[tid].y, smem2[tid + 256].y);
		smem2[tid].z = min(smem2[tid].z, smem2[tid + 256].z);
	}
	__syncthreads();
	if (blockDim.x >= 256 && tid < 128)
	{
		smem1[tid].x = max(smem1[tid].x, smem1[tid + 128].x);
		smem1[tid].y = max(smem1[tid].y, smem1[tid + 128].y);
		smem1[tid].z = max(smem1[tid].z, smem1[tid + 128].z);

		smem2[tid].x = min(smem2[tid].x, smem2[tid + 128].x);
		smem2[tid].y = min(smem2[tid].y, smem2[tid + 128].y);
		smem2[tid].z = min(smem2[tid].z, smem2[tid + 128].z);
	}
	__syncthreads();
	if (blockDim.x >= 128 && tid < 64)
	{
		smem1[tid].x = max(smem1[tid].x, smem1[tid + 64].x);
		smem1[tid].y = max(smem1[tid].y, smem1[tid + 64].y);
		smem1[tid].z = max(smem1[tid].z, smem1[tid + 64].z);

		smem2[tid].x = min(smem2[tid].x, smem2[tid + 64].x);
		smem2[tid].y = min(smem2[tid].y, smem2[tid + 64].y);
		smem2[tid].z = min(smem2[tid].z, smem2[tid + 64].z);
	}
	__syncthreads();
	// unrolling warp
	if (tid < 32)
	{
		volatile typename vec3<T>::Type *vsmem1 = smem1;
		volatile typename vec3<T>::Type *vsmem2 = smem2;
		vsmem1[tid].x = max(vsmem1[tid].x, vsmem1[tid + 32].x);
		vsmem1[tid].y = max(vsmem1[tid].y, vsmem1[tid + 32].y);
		vsmem1[tid].z = max(vsmem1[tid].z, vsmem1[tid + 32].z);
		vsmem1[tid].x = max(vsmem1[tid].x, vsmem1[tid + 16].x);
		vsmem1[tid].y = max(vsmem1[tid].y, vsmem1[tid + 16].y);
		vsmem1[tid].z = max(vsmem1[tid].z, vsmem1[tid + 16].z);
		vsmem1[tid].x = max(vsmem1[tid].x, vsmem1[tid +  8].x);
		vsmem1[tid].y = max(vsmem1[tid].y, vsmem1[tid +  8].y);
		vsmem1[tid].z = max(vsmem1[tid].z, vsmem1[tid +  8].z);
		vsmem1[tid].x = max(vsmem1[tid].x, vsmem1[tid +  4].x);
		vsmem1[tid].y = max(vsmem1[tid].y, vsmem1[tid +  4].y);
		vsmem1[tid].z = max(vsmem1[tid].z, vsmem1[tid +  4].z);
		vsmem1[tid].x = max(vsmem1[tid].x, vsmem1[tid +  2].x);
		vsmem1[tid].y = max(vsmem1[tid].y, vsmem1[tid +  2].y);
		vsmem1[tid].z = max(vsmem1[tid].z, vsmem1[tid +  2].z);
		vsmem1[tid].x = max(vsmem1[tid].x, vsmem1[tid +  1].x);
		vsmem1[tid].y = max(vsmem1[tid].y, vsmem1[tid +  1].y);
		vsmem1[tid].z = max(vsmem1[tid].z, vsmem1[tid +  1].z);

		vsmem2[tid].x = min(vsmem2[tid].x, vsmem2[tid + 32].x);
		vsmem2[tid].y = min(vsmem2[tid].y, vsmem2[tid + 32].y);
		vsmem2[tid].z = min(vsmem2[tid].z, vsmem2[tid + 32].z);
		vsmem2[tid].x = min(vsmem2[tid].x, vsmem2[tid + 16].x);
		vsmem2[tid].y = min(vsmem2[tid].y, vsmem2[tid + 16].y);
		vsmem2[tid].z = min(vsmem2[tid].z, vsmem2[tid + 16].z);
		vsmem2[tid].x = min(vsmem2[tid].x, vsmem2[tid +  8].x);
		vsmem2[tid].y = min(vsmem2[tid].y, vsmem2[tid +  8].y);
		vsmem2[tid].z = min(vsmem2[tid].z, vsmem2[tid +  8].z);
		vsmem2[tid].x = min(vsmem2[tid].x, vsmem2[tid +  4].x);
		vsmem2[tid].y = min(vsmem2[tid].y, vsmem2[tid +  4].y);
		vsmem2[tid].z = min(vsmem2[tid].z, vsmem2[tid +  4].z);
		vsmem2[tid].x = min(vsmem2[tid].x, vsmem2[tid +  2].x);
		vsmem2[tid].y = min(vsmem2[tid].y, vsmem2[tid +  2].y);
		vsmem2[tid].z = min(vsmem2[tid].z, vsmem2[tid +  2].z);
		vsmem2[tid].x = min(vsmem2[tid].x, vsmem2[tid +  1].x);
		vsmem2[tid].y = min(vsmem2[tid].y, vsmem2[tid +  1].y);
		vsmem2[tid].z = min(vsmem2[tid].z, vsmem2[tid +  1].z);
	}
	// write result for this block to global mem
	if (tid == 0)
	{
		g_odata[blockIdx.x * 2    ] = smem1[0]; // max
		g_odata[blockIdx.x * 2 + 1] = smem2[0]; // min
	}

	if (gridDim.x > 1)
	{
		const unsigned int tid = threadIdx.x;
		__shared__ bool amLast;

		// wait until all outstanding memory instructions in this thread are finished
		__threadfence();

		// Thread 0 takes a ticket
		if (tid==0)
		{
			unsigned int ticket = atomicInc(&retirementCount, gridDim.x);
			// If the ticket ID is equal to the number of blocks, we are the last block!
			amLast = (ticket == gridDim.x-1);
		}

		cg::sync(cta);

		// The last block sums the results of all other blocks
		if (amLast)
		{

			int i = tid;
			volatile typename vec3<T>::Type __shared__ myMax;
			volatile typename vec3<T>::Type __shared__ myMin;

			myMax.x = g_odata[0].x;
			myMax.y = g_odata[0].y;
			myMax.z = g_odata[0].z;
			myMin.x = g_odata[1].x;
			myMin.y = g_odata[1].y;
			myMin.z = g_odata[1].z;

			while (i < gridDim.x )
			{
				myMax.x = max(myMax.x, g_odata[i * 2].x);
				myMax.y = max(myMax.y, g_odata[i * 2].y);
				myMax.z = max(myMax.z, g_odata[i * 2].z);

				myMin.x = min(myMin.x, g_odata[i * 2 + 1].x);
				myMin.y = min(myMin.y, g_odata[i * 2 + 1].y);
				myMin.z = min(myMin.z, g_odata[i * 2 + 1].z);
				i += blockDim.x;
			}

			cg::sync(cta);

			if( tid == 0)
			{
				T scale = 0;
				scale = max(scale, myMax.x - myMin.x);
				scale = max(scale, myMax.y - myMin.y);
				scale = max(scale, myMax.z - myMin.z);

				out[0].x = 0.5f * (myMax.x + myMin.x) - scale/2;
				out[0].y = 0.5f * (myMax.y + myMin.y) - scale/2;
				out[0].z = 0.5f * (myMax.z + myMin.z) - scale/2;
				out[0].w = scale;

				// reset retirement count so that next run succeeds
				retirementCount = 0;
			}
		}
	}
	else
	{
		if( tid == 0)
		{
			typename vec3<T>::Type myMax;
			typename vec3<T>::Type myMin;

			myMax.x = g_odata[0].x;
			myMax.y = g_odata[0].y;
			myMax.z = g_odata[0].z;
			myMin.x = g_odata[1].x;
			myMin.y = g_odata[1].y;
			myMin.z = g_odata[1].z;

			T scale = 0;
			scale = max(scale, myMax.x - myMin.x);
			scale = max(scale, myMax.y - myMin.y);
			scale = max(scale, myMax.z - myMin.z);

			out[0].x = 0.5f * (myMax.x + myMin.x) - scale/2;
			out[0].y = 0.5f * (myMax.y + myMin.y) - scale/2;
			out[0].z = 0.5f * (myMax.z + myMin.z) - scale/2;
			out[0].w = scale;

		}
	}
}

template<typename T>
void getBoundingBox(typename vec4<T>::Type *__restrict__ out, typename vec4<T>::Type *__restrict__ g_idata, unsigned int n)
{
	const int blockSize = 256;
	dim3 block (blockSize, 1);
	dim3 grid (max(1, (n + block.x - 1)/block.x/8), 1);//Unroll 8 -> need at least 2048 particles -> add max to make sure at least 1 block will be launched
	unsigned int retCnt = 0;
	setRetirementCount(retCnt);

	typename vec3<T>::Type * g_odata = NULL;
	checkCudaErrors(cudaMalloc((void**)&g_odata, 2 * grid.x * sizeof(typename vec3<T>::Type)));
	reduceMaxSinglePass<T, blockSize><<<grid, block>>>(out, g_idata, g_odata, n);
	cudaFree(g_odata);
}

template<typename T>
__global__ void
assignMortonCode(unsigned int *__restrict__ morton_code,
		typename vec4<T>::Type *__restrict__ oldPos,
		typename vec4<T>::Type *__restrict__ info,
		unsigned int deviceNumBodies)
{

	unsigned int idx = threadIdx.x + blockIdx.x * blockDim.x;
	typename vec4<T>::Type center; // store center and scale
	center.x = info[0].x + info[0].w/2;
	center.y = info[0].y + info[0].w/2;
	center.z = info[0].z + info[0].w/2;
	center.w = info[0].w;

	if(idx < deviceNumBodies)
	{
		bool region[3];
		morton_code[idx] = 0;

		for(unsigned int i = 0; i < 10; i++)
		{

			region[0] = oldPos[idx].x > center.x;
			region[1] = oldPos[idx].y > center.y;
			region[2] = oldPos[idx].z > center.z;

			unsigned int m_code = (region[0] << 2)
									+ (region[1] << 1)
									+ (region[2]);

			morton_code[idx] += (m_code << (3*(9 - i)));

			center.w /= 2.0;
			center.x += center.w * 0.5 * (region[0] > 0 ? 1 : -1);
			center.y += center.w * 0.5 * (region[1] > 0 ? 1 : -1);
			center.z += center.w * 0.5 * (region[2] > 0 ? 1 : -1);
		}
		if(morton_code[idx] == 0)
		{
			//printf("Morton code can not be 0!\n");
		}
		if(idx < 32)
		{
			//printf("Assign morton code for point %d (%f, %f, %f), center (%f, %f, %f), code: %o\n",
			//		idx, oldPos[idx].x, oldPos[idx].y, oldPos[idx].z, center.x, center.y, center.z, morton_code[idx]);
		}
	}
}

__device__ __forceinline__
int taog(int i, int j, unsigned int* sortedMortonCodes, int numObjects)
{
	if(j < 0||j >= numObjects)
	{
		return -1;
	}
	return __clz(sortedMortonCodes[i] ^ sortedMortonCodes[j]);
}

template<typename T>
__global__ void
generateHierarchy(Node<T>* nodes,
						  unsigned int*__restrict__ sortedMortonCodes,
						  unsigned int*__restrict__ sortedID,
						  //unsigned int*__restrict__ compactStream,
						  typename vec4<T>::Type *__restrict__ Pos,
						  typename vec4<T>::Type *__restrict__ boundingBox,
						  const int	numObjects, bool verbose)
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

	Node<T>* leafNodes = &nodes[numObjects - 1];
	Node<T>* internalNodes = &nodes[0];

    int idx = threadIdx.x + blockIdx.x * blockDim.x;

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

		while ( taog(idx, idx + lmax * d, sortedMortonCodes, numObjects) > dmin )
		{
			lmax <<= 1;
		};

		// Find the other end using binary search

		int l = 0;

		for(int t = (lmax >> 1); t > 0; t >>= 1)
		{
			if (taog(idx, idx + (l + t) * d, sortedMortonCodes, numObjects) > dmin)
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

        unsigned int firstCode = sortedMortonCodes[first];
        unsigned int lastCode = sortedMortonCodes[last];

        int split = first;// initial guess

        if (firstCode == lastCode)
        	split =  (first + last) >> 1;

        // Calculate the number of highest bits that are the same
        // for all objects, using the count-leading-zeros intrinsic.

        int commonPrefix = __clz(firstCode ^ lastCode);

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
                unsigned int splitCode = sortedMortonCodes[newSplit];
                int splitPrefix = __clz(firstCode ^ splitCode);
                if (splitPrefix > commonPrefix)
                    split = newSplit; // accept proposal
            }
        }
        while (step > 1);

        //////////////////////////////////////////////////////////////////////////////////
        // Record parent-child relationships.
        //////////////////////////////////////////////////////////////////////////////////

        // Select childA.

        Node<T>* childA;
        if (split == first)
            childA = &leafNodes[split];
        else
            childA = &internalNodes[split];

        // Select childB.

        Node<T>* childB;
        if (split + 1 == last)
            childB = &leafNodes[split + 1];
        else
            childB = &internalNodes[split + 1];

        // Record parent-child relationships.

        internalNodes[idx].childA = childA;
        internalNodes[idx].childB = childB;

        // Record segementation
        unsigned int level = commonPrefix - 2;
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
        	if(level == 30)// Last level
        	{
        		internalNodes[idx].isLeafNode = true;
        		internalNodes[idx].childA = NULL;
        		internalNodes[idx].childB = NULL;
        	}

        	internalNodes[idx].isOctreeNode = true;
        	internalNodes[idx].bbox = *boundingBox;

        	unsigned int prefix = (firstCode & lastCode) >> (30 - level);
        	for(int i = 0; i < level/3; i++)
        	{
        		T current_scale = internalNodes[idx].bbox.w / (T)(2 << i);
        		internalNodes[idx].bbox.x += ((prefix >> (level - 1 - 3 * i)) & 1) * current_scale ;
        		internalNodes[idx].bbox.y += ((prefix >> (level - 2 - 3 * i)) & 1) * current_scale ;
        		internalNodes[idx].bbox.z += ((prefix >> (level - 3 - 3 * i)) & 1) * current_scale ;
        	}
        	internalNodes[idx].bbox.w /= (T)(1 << (level/3));
        	internalNodes[idx].bbox.x += internalNodes[idx].bbox.w/2 ;
        	internalNodes[idx].bbox.y += internalNodes[idx].bbox.w/2 ;
        	internalNodes[idx].bbox.z += internalNodes[idx].bbox.w/2 ;

        	// Update mass center and total mass
        	internalNodes[idx].mass.x = 0.0f;
        	internalNodes[idx].mass.y = 0.0f;
        	internalNodes[idx].mass.z = 0.0f;
        	internalNodes[idx].mass.w = 0.0f;

        	for(int i = first; i <= last; i++)
        	{
        		typename vec4<T>::Type p = Pos[sortedID[i]];
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
        	}

        }
    }

    // Print tree
    if(verbose)
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
    }
    // Node 0 is the root.
    //return &internalNodes[0];
}

template<typename T>
__device__ typename vec4<T>::Type
computeAccel(typename vec4<T>::Type my,
		typename vec4<T>::Type * posList,
		unsigned int*__restrict__ sortedID,
			Node<T>*        root)
{
    // Allocate traversal stack from thread-local memory,
    // and push NULL to indicate that there are no postponed nodes.
	Node<T>* stack[32];
	Node<T>** stackPtr = stack;

	*stackPtr++ = nullptr; // push

    typename vec4<T>::Type ai = {0.0f, 0.0f, 0.0f, 0.0f}, r;

    // Traverse nodes starting from the root.
    Node<T>* node = root;

    while (node != nullptr)
    {
        // Check each child node for overlap.

        Node<T>* childL = node->childA;
        Node<T>* childR = node->childB;

    	bool overlapL = true;
    	bool overlapR = true;

    	if(childL->isOctreeNode) // is octree node
    	{
    		typename vec4<T>::Type mass = childL->mass;
    		T aabb = childL->bbox.w;
    		r.x = mass.x - my.x;
    		r.y = mass.y - my.y;
    		r.z = mass.z - my.z;
	        T distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
	        distSqr += getSofteningSquared<T>();

	        T invDist = rsqrt_T(distSqr);

    		overlapL = __any_sync(__activemask(), (aabb * invDist) > 0.5); // default parameter for barnes-hut algorithm.
    		if(!overlapL)
    		{
    			T invDistCube = invDist * invDist * invDist;
    			T dirt = mass.w;
    			if(my.w < 0 && dirt > 0) dirt = 0.01f*-dirt;
    			T s = dirt * invDistCube; // mass * invdistcube
    			ai.x += r.x * s;
    			ai.y += r.y * s;
    			ai.z += r.z * s;
    		}
    	}

    	if(childR->isOctreeNode) // is octree node
    	{
    		typename vec4<T>::Type mass = childR->mass;
    		T aabb = childR->bbox.w;
    		r.x = mass.x - my.x;
    		r.y = mass.y - my.y;
    		r.z = mass.z - my.z;
	        T distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
	        distSqr += getSofteningSquared<T>();

	        T invDist = rsqrt_T(distSqr);

    		overlapR = __any_sync(__activemask(), (aabb * invDist) > 0.5); // default parameter for barnes-hut algorithm.
    		if(!overlapR)
    		{
    	        T invDistCube =  invDist * invDist * invDist;
    			T dirt = mass.w;
    			if(my.w < 0 && dirt > 0) dirt = 0.01f*-dirt;
    			T s = dirt * invDistCube; // mass * invdistcube
    			ai.x += r.x * s;
    			ai.y += r.y * s;
    			ai.z += r.z * s;
    		}
    	}

        if (overlapL && childL->isLeafNode)
        {
        	for(int i = childL->begin; i <= childL->end; i++)
        	{
				typename vec4<T>::Type sh = posList[sortedID[i]];//childL->begin

				r.x = sh.x - my.x;
				r.y = sh.y - my.y;
				r.z = sh.z - my.z;
				T distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
				distSqr += getSofteningSquared<T>();

				T invDist = rsqrt_T(distSqr);
				T invDistCube =  invDist * invDist * invDist;
    			T dirt = sh.w;
    			if(my.w < 0 && dirt > 0) dirt = 0.01f*-dirt;
    			T s = dirt * invDistCube; // mass * invdistcube
				ai.x += r.x * s;
				ai.y += r.y * s;
				ai.z += r.z * s;

        	}
        }

        if (overlapR && childR->isLeafNode)
        {
        	for(int i = childR->begin; i <= childR->end; i++)
        	{
				typename vec4<T>::Type sh = posList[sortedID[i]];//childL->begin

				r.x = sh.x - my.x;
				r.y = sh.y - my.y;
				r.z = sh.z - my.z;
				T distSqr = r.x * r.x + r.y * r.y + r.z * r.z;
				distSqr += getSofteningSquared<T>();

				T invDist = rsqrt_T(distSqr);
				T invDistCube =  invDist * invDist * invDist;
    			T dirt = sh.w;
    			if(my.w < 0 && dirt > 0) dirt = 0.01f*-dirt;
    			T s = dirt * invDistCube; // mass * invdistcube
				ai.x += r.x * s;
				ai.y += r.y * s;
				ai.z += r.z * s;
        	}
        }

        //if(ind == 0)  printf("%d BREAKPOINT 5\n", ind);
        // Query overlaps an internal node => traverse.
        bool traverseL = (overlapL && !childL->isLeafNode);
        bool traverseR = (overlapR && !childR->isLeafNode);

        if (!traverseL && !traverseR)
        {
        	node = *--stackPtr; // pop
        }
        else
        {
            node = (traverseL) ? childL : childR;
            if (traverseL && traverseR)
            {
            	*stackPtr++ = childR; // push
            }
        }
    }

    return ai;
}

template<typename T>
__global__ void
integrateBodies(Node<T>* root,
				unsigned int*__restrict__ sortedID,
				typename vec4<T>::Type *__restrict__ newPos,
        		//typename vec4<T>::Type *__restrict__ oldPos,
        		typename vec4<T>::Type *__restrict__ vel,
        		typename vec4<T>::Type *__restrict__ acc,
        		unsigned int deviceOffset, unsigned int deviceNumBodies,
        		float deltaTime, float damping)
{
	// Handle to thread block group
	// utilize shared memory to store local tree for entire block, may include halos, position list etc
	// compute essential tree for each warp or block
	// to improve global load efficiency
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= deviceNumBodies)
	{
		return;
	}

	unsigned int sortedIdx = sortedID[idx];

	typename vec4<T>::Type position = newPos[sortedIdx]; //deviceOffset +
	typename vec4<T>::Type velocity = vel[sortedIdx]; //deviceOffset +

	typename vec4<T>::Type accel = computeAccel<T>(position, newPos, sortedID, root);

	velocity.x += 0.5 * accel.x * deltaTime;
	velocity.y += 0.5 * accel.y * deltaTime;
	velocity.z += 0.5 * accel.z * deltaTime;

	// store new position and velocity
	//newPos[sortedIndex] = position; //deviceOffset +
	vel[sortedIdx]    = velocity;
	acc[sortedIdx]    = accel;
}

template<typename T>
__global__ void
LeapFrog(typename vec4<T>::Type *__restrict__ newPos,
		 typename vec4<T>::Type *__restrict__ oldPos,
		 typename vec4<T>::Type *__restrict__ vel,
		 typename vec4<T>::Type *__restrict__ acc,
		 unsigned int deviceOffset, unsigned int deviceNumBodies,
		 float deltaTime)
{
	// Handle to thread block group
	int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx >= deviceNumBodies)
	{
		return;
	}

	typename vec4<T>::Type position = oldPos[idx]; //deviceOffset +
	typename vec4<T>::Type velocity = vel[idx]; //deviceOffset +
	typename vec4<T>::Type oldAcc = acc[idx]; //deviceOffset +

	// Leapfrog
	velocity.x += 0.5 * oldAcc.x * deltaTime;
	velocity.y += 0.5 * oldAcc.y * deltaTime;
	velocity.z += 0.5 * oldAcc.z * deltaTime;

	position.x += velocity.x * deltaTime;
	position.y += velocity.y * deltaTime;
	position.z += velocity.z * deltaTime;

	vel[idx]    = velocity;
	newPos[idx] = position; //deviceOffset +
}

template <typename T>
void integrateNbodySystem(DeviceData<T> *deviceData,
		cudaGraphicsResource **pgres,
		unsigned int currentRead,
		float deltaTime,
		float damping,
		unsigned int numBodies,
		unsigned int numDevices,
		int blockSize,
		bool bUsePBO)
{
	if (bUsePBO)
	{
		checkCudaErrors(cudaGraphicsResourceSetMapFlags(pgres[currentRead], cudaGraphicsMapFlagsReadOnly));
		checkCudaErrors(cudaGraphicsResourceSetMapFlags(pgres[1-currentRead], cudaGraphicsMapFlagsWriteDiscard));
		checkCudaErrors(cudaGraphicsMapResources(2, pgres, 0));
		size_t bytes;
		checkCudaErrors(cudaGraphicsResourceGetMappedPointer((void **)&(deviceData[0].dPos[currentRead]), &bytes, pgres[currentRead]));
		checkCudaErrors(cudaGraphicsResourceGetMappedPointer((void **)&(deviceData[0].dPos[1-currentRead]), &bytes, pgres[1-currentRead]));
	}

	blockSize = 512;
	int numBlocks = (deviceData[0].numBodies + blockSize-1) / blockSize;

	LeapFrog<T><<<numBlocks, blockSize>>>(
			(typename vec4<T>::Type *)deviceData[0].dPos[1-currentRead],
			(typename vec4<T>::Type *)deviceData[0].dPos[currentRead],
			(typename vec4<T>::Type *)deviceData[0].dVel,
			(typename vec4<T>::Type *)deviceData[0].dAcc,
			deviceData[0].offset, deviceData[0].numBodies,
			deltaTime);

	getBoundingBox<T>((typename vec4<T>::Type *)deviceData[0].dInfo, (typename vec4<T>::Type *)deviceData[0].dPos[1 - currentRead], deviceData[0].numBodies);

	assignMortonCode<T><<<numBlocks, blockSize>>>
				   (deviceData[0].morton_code,
					(typename vec4<T>::Type *)deviceData[0].dPos[1 - currentRead],
					(typename vec4<T>::Type *)deviceData[0].dInfo,
					deviceData[0].numBodies);

	thrust::device_ptr<unsigned int> morton_ptr = thrust::device_pointer_cast(deviceData[0].morton_code);
	thrust::device_ptr<unsigned int> sort_ptr = thrust::device_pointer_cast(deviceData[0].sortedID);

	thrust::sequence(sort_ptr, sort_ptr + deviceData[0].numBodies);
	thrust::sort_by_key(morton_ptr, morton_ptr + deviceData[0].numBodies, sort_ptr);

	generateHierarchy<T><<<numBlocks, blockSize>>>
				   (deviceData[0].root,
					deviceData[0].morton_code,
					deviceData[0].sortedID,
					//deviceData[0].compactStream,
					(typename vec4<T>::Type *)deviceData[0].dPos[1 - currentRead],
					(typename vec4<T>::Type *)deviceData[0].dInfo,
					deviceData[0].numBodies, false);

	integrateBodies<T><<< numBlocks, blockSize>>>
				   (deviceData[0].root,
					deviceData[0].sortedID,
					(typename vec4<T>::Type *)deviceData[0].dPos[1-currentRead],
					(typename vec4<T>::Type *)deviceData[0].dVel,
					(typename vec4<T>::Type *)deviceData[0].dAcc,
					deviceData[0].offset, deviceData[0].numBodies,
					deltaTime, damping);

	for (unsigned int dev = 0; dev != numDevices; dev++)
	{
		if (numDevices > 1)
		{
			cudaSetDevice(dev);
		}

		if (numDevices > 1)
		{
			checkCudaErrors(cudaEventRecord(deviceData[dev].event));
			// MJH: Hack on older driver versions to force kernel launches to flush!
			cudaStreamQuery(0);
		}

		// check if kernel invocation generated an error
		getLastCudaError("Kernel execution failed");
	}

	if (numDevices > 1)
	{
		for (unsigned int dev = 0; dev < numDevices; dev++)
		{
			checkCudaErrors(cudaEventSynchronize(deviceData[dev].event));
		}
	}

	if (bUsePBO)
	{
		checkCudaErrors(cudaGraphicsUnmapResources(2, pgres, 0));
	}
}


// Explicit specializations needed to generate code
template void integrateNbodySystem<float>(DeviceData<float> *deviceData,
		cudaGraphicsResource **pgres,
		unsigned int currentRead,
		float deltaTime,
		float damping,
		unsigned int numBodies,
		unsigned int numDevices,
		int blockSize,
		bool bUsePBO);

template void integrateNbodySystem<double>(DeviceData<double> *deviceData,
		cudaGraphicsResource **pgres,
		unsigned int currentRead,
		float deltaTime,
		float damping,
		unsigned int numBodies,
		unsigned int numDevices,
		int blockSize,
		bool bUsePBO);
