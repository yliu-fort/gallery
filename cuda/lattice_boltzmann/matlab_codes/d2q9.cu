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

// Periodic fetching
__device__ __forceinline__
int fetchStencil(int2 ij, int dirx, int diry) {

    return (ij.x+dirx + tileX*(ij.y+diry) + nElem)%nElem;

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
        /*double req = r[globalThreadIdx];
        double ueq = u[globalThreadIdx];
        double veq = v[globalThreadIdx];
        // no-slip (Bounce_back)
        for(int i = 0; i < nDir; i++)
        {
            double cidotw = dirX[i]*ueq + dirY[i]*veq;

         // rhow ignored for incompressible flow
            out[globalThreadIdx + bi[i]*nElem] =
                in[neighbours[i] + i*nElem] - 2.0*wi[i]*req*(3.0*cidotw);
        }*/

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

}

__device__ __forceinline__ double collision(
                          const double& f,
                          const double& feq,
                          const double& a,
                          const double& b) {

        return (f - a*b*(f - feq));

}

__device__ __forceinline__ double compute_g(
                          const double * f,
                          const double * feq,
                          const double& a,const double& b) {

        double result = 0.0;
        for(int i = 0;i < nDir; i++)
        {
            double c = collision(f[i],feq[i],a,b);
            if(c < 0.0) {c = 0.00000000001;}
            result += c*log(c/wi[i]) - f[i]*log(f[i]/wi[i]);
        }
        return result;
}

__device__ __forceinline__ double compute_gradg(
                          const double * f,
                          const double * feq,
                          const double& a,const double& b)
{

        double result = 0.0;
        for(int i = 0;i < nDir; i++)
        {
            double c = collision(f[i],feq[i],a,b);
            if(c < 0.0) {c = 0.00000000001;}
            result += -b*(f[i] - feq[i])*(log(c/wi[i]) + 1.0);
        }
        return result;
}

__device__ __forceinline__ void swap(double& a, double& b)
{

            double tmp = a;
            a = b;
            b = tmp;
}

__constant__ double stableDeviation;
__constant__ double alphaMin;
__constant__ bool enableEntropyConstraint;
__device__ double constrain_entropy(
                          const double * f,
                          const double * feq,const double& b) {
        // calculate deviation
        double amin=alphaMin, amax=2.0;
        double maxDeviation = 0.0;
        for(int i = 0;i < nDir; i++)
        {
            double deviation = abs(f[i]-feq[i])/feq[i];
            if(deviation > maxDeviation)
                maxDeviation = deviation;
        }

        // if deviation is too large
        //double stableDeviation = 0.2;
        if(maxDeviation < stableDeviation) return amax;

        // compute G value
        double Gmin = compute_g(f,feq,amin,b);
        double Gmax = compute_g(f,feq,amax,b);
        double gradGmin = compute_gradg(f,feq,amin,b);
        double gradGmax = compute_gradg(f,feq,amax,b);
        if(Gmin*Gmax > 0) return amax;
        if(Gmin > 0) swap(amin, amax);

        double a = 0.5*(amin + amax);
        double da = abs(amax - amin);
        double a_o = a;
        //double da_o = da;
        double G = compute_g(f,feq,a,b);
        double gradG = compute_gradg(f,feq,a,b);

        int maxIter = 20;
        double tolerance = 0.0001;
        for(int it = 0; it < maxIter; it++)
        {
            if( ( ((a-amax)*gradG-G)*((a-amin)*gradG-G) >= 0 )
            ||  ( abs(a_o*gradG-G-1.0) > 1.0 ) )
            {
                // bisection
                //da_o = da;
                da = 0.5*(amax - amin);
                a = amin-amax;
                if(amin == a) return a;
            }else
            {
                //da_o = da;
                da = G/gradG;
                a_o = a;
                a -= da;
                if(a_o == a) return a;
            }
            if(abs(da) < tolerance) return a;

            G = compute_g(f,feq,a,b);
            gradG = compute_gradg(f,feq,a,b);
            if(G < 0.0) {amin = a;}
            else {amax = a;}
        }

        return amax;

}

__constant__ double constant_Fx[9]; // Fi = wi*c_ia*Fa/cs^2
__constant__ double constant_Fy[9];

__global__ void processCollision(
                      double * out,
                      const double * in,
                      const double * r,
                      const double * u,
                      const double * v,
                      const bool * mask )
{
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
    double f[9];
    double feq[9];
    for(int i = 0; i < nDir; i++)
    {
        f[i] = in[globalThreadIdx + i*nElem];
        feq[i] = compute_equilibrium(rho_eq, u_eq, v_eq, i);
    }

    // Entropic LBM implementation
    double alpha = 2.0;
    if(enableEntropyConstraint)
    	alpha = constrain_entropy(f, feq, 0.5*tau[1]);

    for(int i = 0; i < nDir; i++)
    {
        // BGK Collider
        f[i] = collision(f[i],feq[i],alpha,0.5*tau[1]);

        out[globalThreadIdx + i*nElem] = f[i];
    }
    //out2[globalThreadIdx] = alpha;

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

__global__ void processNoslip(
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
    if(mask[globalThreadIdx] == false) return;

    // Get our X and Y coords
    int neighbours[9];
    streaming(neighbours, ij);

    // Processing boundary nodes
    double req = r[globalThreadIdx];
    double ueq = u[globalThreadIdx];
    double veq = v[globalThreadIdx];

    // no-slip (Bounce_back)
    for(int i = 0; i < nDir; i++)
    {
        double cidotw = dirX[i]*ueq + dirY[i]*veq;

      //rhow ignored for incompressible flow
        out[globalThreadIdx + bi[i]*nElem] =
                in[neighbours[i] + i*nElem] - 2.0*wi[i]*req*(3.0*cidotw);
    }

}

__device__ double one_sided_diff(double a1, double a2, double a3, double dir)
{
    return dir*(0.5*(-3.0*a1+4.0*a2-a3));
}

__device__ double central_diff(double ap, double am)
{
    return 0.5*(ap-am);
}


__global__ void processNRoutlet(
              double * out,
              const double * r,
              const double * u,
              const double * v,
              const bool * mask) {
    // Work out which thread we are
    size_t const globalThreadIdx = calculateGlobalIndex();
    int2 ij = fetch2D(globalThreadIdx);

    // If we're off the end, return now
    if (globalThreadIdx >= nElem) return;
    if(mask[globalThreadIdx] == false) return;

    // Read inputs
    double r1 = r[globalThreadIdx];
    double u1 = u[globalThreadIdx];
    double v1 = v[globalThreadIdx];

    // Read stencil
    int x1 = fetchStencil(ij,-1, 0);
    int x2 = fetchStencil(ij,-2, 0);
    int y1 = fetchStencil(ij, 0, 1);
    int y2 = fetchStencil(ij, 0,-1);

    double drdx = one_sided_diff(r1,r[x1],r[x2],-1);
    double dudx = one_sided_diff(u1,u[x1],u[x2],-1);
    double dvdx = one_sided_diff(v1,v[x1],v[x2],-1);

    double drdy = central_diff(r[y1],r[y2]);
    double dudy = central_diff(u[y1],u[y2]);
    double dvdy = central_diff(v[y1],v[y2]);

    double L1 = 0;
    double L2 = u1*dvdx;
    double L3 = (u1 + cs[3])*(cs[1]*drdx + cs[3]*r1*dudx);

    // dmdt = -Pxinv*Lx'
    double P1 = 0.5*cs[0];
    double P2 = 0.5/r1*cs[2];
    double dm1dt =  P1*L1 + P1*L3;
    double dm2dt = -P2*L1 + P2*L3;
    double dm3dt = L2;

    double Ydm1dy = v1*drdy         + r1*dvdy;
    double Ydm2dy =         v1*dudy;
    double Ydm3dy = cs[1]/r1*drdy   + v1*dvdy;

    // compute new m
    double gamma = 0.0;
    r1 -= (dm1dt + gamma*Ydm1dy);
    u1 -= (dm2dt + gamma*Ydm2dy);
    v1 -= (dm3dt + gamma*Ydm3dy);

    // Processing boundary nodes
    // periodic with no pressure drop
    for(int i = 0; i < nDir; i++)
    {
        double feq = compute_equilibrium(r1,u1,v1, i);

        // replace f by feq
        out[globalThreadIdx + i*nElem] = feq;
    }
    
}

__global__ void processRoutlet(
              double * out,
              const double * in,
              const double * r,
              const double * u,
              const double * v,
              const bool * mask) {
    // Work out which thread we are
    size_t const globalThreadIdx = calculateGlobalIndex();
    int2 ij = fetch2D(globalThreadIdx);

    // If we're off the end, return now
    if (globalThreadIdx >= nElem) return;
    if(mask[globalThreadIdx] == false) return;

    // Read stencil
    int x1 = fetchStencil(ij,-1, 0);
    double r1 = r[x1];
    double u1 = u[x1];
    double v1 = v[x1];

    // Processing boundary nodes
    // periodic with no pressure drop
    for(int i = 0; i < nDir; i++)
    {
        //double feq = compute_equilibrium(r1,u1,v1, i);

        // replace f by feq
        out[globalThreadIdx + i*nElem] = in[x1 + i*nElem];
    }
    
}