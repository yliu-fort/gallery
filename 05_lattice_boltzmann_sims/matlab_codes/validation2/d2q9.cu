/**
* @file streamingd2q9.cu
*
* CUDA code to perform 2D LBM simulation on a GPU.
*
* Copyright Yuxuan Liu 2019
*/

/** Work out which piece of the global array this thread should operate on */
__device__ size_t calculateGlobalIndex() {
    // Which block are we?
    size_t const globalBlockIndex = blockIdx.x + blockIdx.y * gridDim.x;
    // Which thread are we within the block?
    size_t const localThreadIdx = threadIdx.x + blockDim.x * threadIdx.y;
    // How big is each block?
    size_t const threadsPerBlock = blockDim.x*blockDim.y;
    // Which thread are we overall?
    return localThreadIdx + globalBlockIndex*threadsPerBlock;

}

/** Kernel constants */
__constant__ int tileX;
__constant__ int tileY;
__constant__ int nElem;
__constant__ int nDir;
__constant__ double wi[9];
__constant__ int dirX[9];
__constant__ int dirY[9];
__constant__ int bi[9];

__constant__ double cs[4];// csSqrInv, csSqr, csInv, cs
__constant__ double tau[3];// tau, 1/tau, 1-1/tau

// Periodic fetching
__device__ __forceinline__
int2 fetch2D(int ind) {

    int2 ij;
    ij.x = (ind%nElem)%tileX;
    ij.y = (ind%nElem)/tileX;
    return ij;

};

// Periodic fetching
__device__ __forceinline__
int fetchInd(int2 ij) {

    return (ij.x + tileX*ij.y + nElem)%nElem;

};

/** Streaming kernal */
__device__ void streaming(int * out, const int2& uv)
{
    for(int i = 0; i < nDir; i++)
    {
        int2 dir;
        dir.x = uv.x - dirX[i];
        dir.y = uv.y - dirY[i];
        out[i] = fetchInd(dir);
    }
}

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 */
__global__ void processStreaming(
                      double * out,
                      const double * in,
                      const double * r,
                      const double * u,
                      const double * v,
                      const bool * mask ) {
    // Work out which thread we are
    size_t const globalThreadIdx = calculateGlobalIndex();
    int2 ij = fetch2D(globalThreadIdx);

    // If we're off the end, return now
    if (globalThreadIdx >= nElem) {
        return;
    }

    // Get our X and Y coords
    int neighbours[9];
    streaming(neighbours, ij);

    // Processing boundary nodes
    if(mask[globalThreadIdx] == true) {
        // no-slip (Bounce_back)
        for(int i = 0; i < nDir; i++)
        {
            double cidotw = dirX[i]*u[globalThreadIdx]
                          + dirY[i]*v[globalThreadIdx];

         // rhow ignored for incompressible flow
            out[globalThreadIdx + bi[i]*nElem] =
                    in[neighbours[i] + i*nElem] 
                            - 2.0*wi[i]*r[globalThreadIdx]*(3.0*cidotw);
        }

        return;
    }

    // Processing internal nodes
    for(int i = 0; i < nDir; i++)
    {
        out[globalThreadIdx + i*nElem] = in[neighbours[i] + i*nElem];
    }

}

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 */
__device__ double compute_equilibrium(
                          const double& r,
                          const double& u,
                          const double& v,
                          const int& dim) {

        // calculate dot product
        double cidotu = dirX[dim]*u + dirY[dim]*v;

        // calculate equilibrium
        return wi[dim]*r*(1.0 + 3.0*cidotu+4.5*cidotu*cidotu-1.5*(u*u+v*v));

        // BGK collider
        //out = tau[2]*out + tau[1]*feq;

}

__constant__ double constant_Fx[9]; // Fi = wi*c_ia*Fa/cs^2
__constant__ double constant_Fy[9];

__global__ void processCollision(
                      double * out,
                      const double * in,
                      const double * r,
                      const double * u,
                      const double * v,
                      const bool * mask ) {
    // Work out which thread we are
    size_t const globalThreadIdx = calculateGlobalIndex();
    int2 ij = fetch2D(globalThreadIdx);

    // If we're off the end, return now
    if (globalThreadIdx >= nElem) {
        return;
    }

    // no collision for boundary nodes
    if(mask[globalThreadIdx] == true) {
        return;
    }

    // Read inputs
    double rho_eq = r[globalThreadIdx];
    double u_eq = u[globalThreadIdx];
    double v_eq = v[globalThreadIdx];

    // Force: Shan&Chen scheme
    for(int i = 0; i < nDir; i++)
    {
        u_eq += (tau[0]-0.5)*constant_Fx[i]*dirX[i]/rho_eq;
        v_eq += (tau[0]-0.5)*constant_Fy[i]*dirY[i]/rho_eq;
    }

    // Compute collision
    for(int i = 0; i < nDir; i++)
    {
        double f = in[globalThreadIdx + i*nElem];

        double feq = compute_equilibrium(rho_eq, u_eq, v_eq, i);

        // BGK Collider
        f = tau[2]*f + tau[1]*feq;

        out[globalThreadIdx + i*nElem] = f;
    }

}

/** Main entry point.
 * Works out where the current thread should read/write to global memory
 * and calls doIterations to do the actual work.
 */
__global__ void processCompute(
                      double * r,
                      double * u,
                      double * v,
                      const double * fin,
                      const bool * mask ) {
    // Work out which thread we are
    size_t const globalThreadIdx = calculateGlobalIndex();
    int2 ij = fetch2D(globalThreadIdx);

    // If we're off the end, return now
    if (globalThreadIdx >= nElem) {
        return;
    }

    // no collision for boundary nodes
    if(mask[globalThreadIdx] == true) {
        return;
    }

    // compute r, u, v
    double _rho=0, _u=0, _v=0;
    for(int i = 0; i < nDir; i++)
    {
        double f = fin[globalThreadIdx + i*nElem];
        _rho += f;
        _u += dirX[i]*f;
        _v += dirY[i]*f;

        // Force contribution
        _u += 0.5*constant_Fx[i]*dirX[i];
        _v += 0.5*constant_Fy[i]*dirY[i];
    }

    _u /= _rho;
    _v /= _rho;

    // Output
    r[globalThreadIdx] = _rho;
    u[globalThreadIdx] = _u;
    v[globalThreadIdx] = _v;

}

__global__ void processPeriodic(
              double * out,
              const double * in,
              const double * a,
              const double * rp,
              const double * r,
              const double * u,
              const double * v,
              const int nBElem) {
    // Work out which thread we are
    size_t const globalThreadIdx = calculateGlobalIndex();
    int2 ij = fetch2D(globalThreadIdx);

    // If we're off the end, return now
    if (globalThreadIdx >= nBElem) {
        return;
    }

    // Read inputs
    double alpha = a[globalThreadIdx];
    double rhop_eq = rp[globalThreadIdx];
    double rho_eq = r[globalThreadIdx];
    double u_eq = u[globalThreadIdx];
    double v_eq = v[globalThreadIdx];
    double u_eq_s = alpha*u_eq;
    double v_eq_s = alpha*v_eq;
    
    // Processing boundary nodes
        // periodic with no pressure drop
        for(int i = 0; i < nDir; i++)
        {
            
            double f = in[globalThreadIdx + i*nBElem];
            double feq = compute_equilibrium(rho_eq, u_eq, v_eq, i);
            double feq_s = compute_equilibrium(rhop_eq, u_eq_s, v_eq_s, i);

            // rhow ignored for incompressible flow
            out[globalThreadIdx + i*nBElem] = f - feq + feq_s;
        }

}