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
                          const double alphaOld,
                          const double * f,
                          const double * feq,const double& b) {
        // calculate deviation
        double amin=alphaMin, amax=alphaOld;
        double maxDeviation = 0.0;
        for(int i = 0;i < nDir; i++)
        {
            double deviation = abs(f[i]-feq[i])/feq[i];
            if(deviation > maxDeviation)
                maxDeviation = deviation;
        }

        // if deviation is too large
        //double stableDeviation = 0.2;
        if(maxDeviation < stableDeviation) return min(2.0,1.01*amax);

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

//__constant__ double constant_Fx[9]; // Fi = wi*c_ia*Fa/cs^2
//__constant__ double constant_Fy[9];

__global__ void processCollision(
                      double * out,
                      double * alpha_in,
                      const double * in,
                      const double * r,
                      const double * u,
                      const double * v,
                      const double * fx,
                      const double * fy,
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
    double _fx = fx[globalThreadIdx];
    double _fy = fy[globalThreadIdx];
    //u_eq += tau[0]*fx[globalThreadIdx]/rho_eq;
    //v_eq += tau[0]*fy[globalThreadIdx]/rho_eq;

    u_eq += 0.5*_fx/rho_eq;
    v_eq += 0.5*_fy/rho_eq;

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
    	alpha = constrain_entropy(alpha_in[globalThreadIdx], f, feq, 0.5*tau[1]);

    for(int i = 0; i < nDir; i++)
    {
        // Forcing
        double cidotu = dirX[i]*u_eq + dirY[i]*v_eq;
        double Si = (1 - 0.5*tau[1])*wi[i]*((dirX[i] - u_eq)*cs[0] + 
                cidotu*dirX[i]*cs[0]*cs[0])*_fx + 
                    (1 - 0.5*tau[1])*wi[i]*((dirY[i] - v_eq)*cs[0] + 
                cidotu*dirY[i]*cs[0]*cs[0])*_fy;

        // BGK Collider
        f[i] = collision(f[i],feq[i],alpha,0.5*tau[1]) + Si;

        out[globalThreadIdx + i*nElem] = f[i];
    }

    // Record alpha value
    alpha_in[globalThreadIdx] = alpha;

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
        if(f<0) continue; // add for stability
        _rho += f;
        _u += dirX[i]*f;
        _v += dirY[i]*f;
    }

    _u /= _rho;
    _v /= _rho;

    // Output
    r[globalThreadIdx] = _rho;
    u[globalThreadIdx] = _u;
    v[globalThreadIdx] = _v;

}

// IB node forcing
__device__ __forceinline__ double 
biLinear(double q11, double q12, double q21, double q22, int x1, int x2, int y1, int y2, float x, float y) 
{
    double x2x1, y2y1, x2x, y2y, yy1, xx1;
    x2x1 = x2 - x1; // = 1
    y2y1 = y2 - y1; // = 1
    x2x = x2 - x;
    y2y = y2 - y;
    yy1 = y - y1;
    xx1 = x - x1;
    return 1.0 / (x2x1 * y2y1) * (
        q11 * x2x * y2y +
        q21 * xx1 * y2y +
        q12 * x2x * yy1 +
        q22 * xx1 * yy1
    );
}

 __global__ void processInterpIbNode(
                      double * r_star,
                      double * u_star,
                      double * v_star,
                      double * r,
                      double * u,
                      double * v,
                      const double * fx,
                      const double * fy,
                      const double * px,
                      const double * py,
                      const int np) {
    // Work out which thread we are
    size_t const globalThreadIdx = calculateGlobalIndex(); // N-th L-point

    // If we're off the end, return now
    if (globalThreadIdx >= np) {
        // Todo: edge cases
        return;
    }

    float2 p_ij;
    p_ij.x = px[globalThreadIdx];
    p_ij.y = py[globalThreadIdx];
    int2 ij0, ij1, ij2, ij3;
    ij0.x = floor(p_ij.x);   ij0.y = floor(p_ij.y);
    ij1.x = floor(p_ij.x)+1; ij1.y = floor(p_ij.y);
    ij2.x = floor(p_ij.x);   ij2.y = floor(p_ij.y)+1;
    ij3.x = floor(p_ij.x)+1; ij3.y = floor(p_ij.y)+1;
    size_t const ind0 = fetchInd(ij0);
    size_t const ind1 = fetchInd(ij1);
    size_t const ind2 = fetchInd(ij2);
    size_t const ind3 = fetchInd(ij3);

    // Interp from Eulerian Lattice
    double r0[4];
    r0[0] = r[ind0];r0[1] = r[ind1];
    r0[2] = r[ind2];r0[3] = r[ind3];
    double u0[4];
    u0[0] = u[ind0];u0[1] = u[ind1];
    u0[2] = u[ind2];u0[3] = u[ind3];
    double v0[4];
    v0[0] = v[ind0];v0[1] = v[ind1];
    v0[2] = v[ind2];v0[3] = v[ind3];

    double _r = biLinear(r0[0],r0[2],r0[1],r0[3], ij0.x, ij0.x+1,ij0.y,ij0.y+1,p_ij.x,p_ij.y);
    double _u = biLinear(u0[0] + 0.5*fx[ind0]/r0[0],
                         u0[2] + 0.5*fx[ind2]/r0[2],
                         u0[1] + 0.5*fx[ind1]/r0[1],
                         u0[3] + 0.5*fx[ind3]/r0[3], ij0.x, ij0.x+1,ij0.y,ij0.y+1,p_ij.x,p_ij.y);
    double _v = biLinear(v0[0] + 0.5*fy[ind0]/r0[0],
                         v0[2] + 0.5*fy[ind2]/r0[2],
                         v0[1] + 0.5*fy[ind1]/r0[1],
                         v0[3] + 0.5*fy[ind3]/r0[3], ij0.x, ij0.x+1,ij0.y,ij0.y+1,p_ij.x,p_ij.y);

    // Output
    r_star[globalThreadIdx] = _r;
    u_star[globalThreadIdx] = _u;
    v_star[globalThreadIdx] = _v;

}

__global__ void processComputeIbNode(
                      double * fx,
                      double * fy,
                      const double* u_star,
                      const double* v_star,
                      const double* px,
                      const double* py,
                      const double dr,
                      const int np) {
    // Work out which thread we are
    size_t const globalThreadIdx = calculateGlobalIndex(); // Nth L-point
            
                // If we're off the end, return now
    if (globalThreadIdx < np) 
    {

        double _px = px[globalThreadIdx];
        double _py = py[globalThreadIdx];

        double _ustar = u_star[globalThreadIdx];
        double _vstar = v_star[globalThreadIdx];

        // Loop all neighbours
        for(int nb = 0; nb < 16; nb++)
        {
            // Compute E-neighbour coords
            int2 ij;
            ij.x = floor(_px) + nb%4-1;
            ij.y = floor(_py) + nb/4-1;
            int ind = fetchInd(ij);

            // Compute distance to E-neighbour
            double d = sqrt((ij.x - _px)*(ij.x - _px) 
                          + (ij.y - _py)*(ij.y - _py));
            if(d < 2.0) 
            {

                // compute L integration
                double dh = dr*0.25*(1.0+cos(3.141592654*d/2.0)); // assume dr = 1.0

                // Assume ud = 0...
                // Output, atomic add
                // Caution: add flag -arch=sm_61 when compile code
                double xcorr = -_ustar*dh;
                double ycorr = -_vstar*dh;

                //fx[ind] += xcorr;
                //fy[ind] += ycorr;

                atomicAdd(&fx[ind], xcorr);
                atomicAdd(&fy[ind], ycorr);
            }
        }
    }

}

__global__ void processComputeIbNodeS(
                      double * fx,
                      double * fy,
                      const double u_star,
                      const double v_star,
                      const double px,
                      const double py,
                      const int pid,
                      const double dr) {
    // Work out which thread we are
    size_t const globalThreadIdx = calculateGlobalIndex(); // 0 - 15
    int2 ij;
    ij.x = floor(px) + globalThreadIdx%4-1;
    ij.y = floor(py) + globalThreadIdx/4-1;
    size_t const ind = fetchInd(ij);

    // If we're off the end, return now
    if (globalThreadIdx >= 16) {
        // Todo: edge cases
        return;
    }

    double d = sqrt((ij.x - px)*(ij.x - px) 
                  + (ij.y - py)*(ij.y - py));
    if(d > 2.0) return; // out of range

    // compute L integration
    double dh = dr*0.25*(1+cos(3.1415926*d/2.0)); // assume dr = 1.0

    // Assume ud = 0...
    // Output
    double xcorr = -u_star*dh;
    double ycorr = -v_star*dh;

    fx[ind] += xcorr;
    fy[ind] += ycorr;

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
    if(~mask[globalThreadIdx]) return;

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
    //double r1 = r[x1];
    //double u1 = u[x1];
    //double v1 = v[x1];

    // Processing boundary nodes
    // periodic with no pressure drop
    for(int i = 0; i < nDir; i++)
    {
        //double feq = compute_equilibrium(r1,u1,v1, i);

        // replace f by feq
        out[globalThreadIdx + i*nElem] = in[x1 + i*nElem];
    }
    
}