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

#include <stdio.h>
#include <stdlib.h>

#include <cuda_runtime.h>
#include <cufft.h>          // CUDA FFT Libraries
#include <helper_cuda.h>    // Helper functions for CUDA Error handling

// OpenGL Graphics includes
#define HELPERGL_EXTERN_GL_FUNC_IMPLEMENTATION
#include <helper_gl.h>


// FluidsGL CUDA kernel definitions
#include "fluidsGL_kernels.cuh"

// Texture reference for reading velocity field
texture<float2, 2> texref;
static cudaArray *array = NULL;

texture<float2, 2> tex_n_1_ref;
static cudaArray *array_n_1 = NULL;

// Particle data
extern GLuint vbo;
extern GLuint vbo_color;                 // OpenGL vertex buffer object
extern struct cudaGraphicsResource *cuda_vbo_resource;
extern struct cudaGraphicsResource *cuda_vbo_color_resource; // handles OpenGL-CUDA exchange

// Texture pitch
extern size_t tPitch; // Should be no difference since we allocate same type & amount of memory
extern cufftHandle planr2c;
extern cufftHandle planc2r;
cData *vxfield = NULL;
cData *vyfield = NULL;

void setupTexture(int x, int y)
{
    // Wrap mode appears to be the new default
    texref.normalized = 1; // Use normalized addressing to validate wrap mode
    texref.filterMode = cudaFilterModeLinear;
    texref.addressMode[0] = cudaAddressModeWrap;
    texref.addressMode[1] = cudaAddressModeWrap;
    texref.addressMode[2] = cudaAddressModeWrap;
    cudaChannelFormatDesc desc = cudaCreateChannelDesc<float2>();

    cudaMallocArray(&array, &desc, y, x);
    getLastCudaError("cudaMalloc failed");

    tex_n_1_ref.normalized = 1; // Use normalized addressing to validate wrap mode
    tex_n_1_ref.filterMode = cudaFilterModeLinear;
    tex_n_1_ref.addressMode[0] = cudaAddressModeWrap;
    tex_n_1_ref.addressMode[1] = cudaAddressModeWrap;
    tex_n_1_ref.addressMode[2] = cudaAddressModeWrap;
    cudaChannelFormatDesc desc_n_1 = cudaCreateChannelDesc<float2>();

    cudaMallocArray(&array_n_1, &desc_n_1, y, x);
    getLastCudaError("cudaMalloc failed");
}

void bindTexture(void)
{
    cudaBindTextureToArray(texref, array);
    getLastCudaError("cudaBindTexture failed");

    cudaBindTextureToArray(tex_n_1_ref, array_n_1);
    getLastCudaError("cudaBindTexture failed");
}

void unbindTexture(void)
{
    cudaUnbindTexture(texref);
    cudaUnbindTexture(tex_n_1_ref);
}

void updateTexture(cData *data, size_t wib, size_t h, size_t pitch)
{
    cudaMemcpy2DToArray(array, 0, 0, data, pitch, wib, h, cudaMemcpyDeviceToDevice);
    getLastCudaError("cudaMemcpy failed");
}

void updateTexture_n_1(cData *data, size_t wib, size_t h, size_t pitch)
{
    cudaMemcpy2DToArray(array_n_1, 0, 0, data, pitch, wib, h, cudaMemcpyDeviceToDevice);
    getLastCudaError("cudaMemcpy failed");
}

void deleteTexture(void)
{
    cudaFreeArray(array);
    cudaFreeArray(array_n_1);
}

// Note that these kernels are designed to work with arbitrary
// domain sizes, not just domains that are multiples of the tile
// size. Therefore, we have extra code that checks to make sure
// a given thread location falls within the domain boundaries in
// both X and Y. Also, the domain is covered by looping over
// multiple elements in the Y direction, while there is a one-to-one
// mapping between threads in X and the tile size in X.
// Nolan Goodnight 9/22/06

// This method adds constant force vectors to the velocity field
// stored in 'v' according to v(x,t+1) = v(x,t) + dt * f.
__global__ void
addForces_k(cData *v, int dx, int dy, int spx, int spy, float fx, float fy, float dt, int r, size_t pitch)
{

    int tx = threadIdx.x;
    int ty = threadIdx.y;
    cData *fj = (cData *)((char *)v + (ty + spy) * pitch) + tx + spx;

    cData vterm = *fj;
    tx -= r;
    ty -= r;
    float s = 1.f / (1.f + tx*tx*tx*tx + ty*ty*ty*ty);
    vterm.x += s * fx;
    vterm.y += s * fy;

    // Bound by CFL
    float CFL = 0.75f;
    cData phiB = make_float2(CFL/(float)dx/dt, CFL/(float)dy/dt);
    vterm.x = max(min(vterm.x, phiB.x), -phiB.x);
    vterm.y = max(min(vterm.y, phiB.y), -phiB.y);

    *fj = vterm;
}

// This method performs the velocity advection step, where we
// trace velocity vectors back in time to update each grid cell.
// That is, v(x,t+1) = v(p(x,-dt),t). Here we perform bilinear
// interpolation in the velocity space.
__global__ void
advectVelocity_k(cData *v, float *vx, float *vy,
                 int dx, int pdx, int dy, float dt, int lb, size_t pitch)
{

    int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
    int gtidy = blockIdx.y * (lb * blockDim.y) + threadIdx.y * lb;
    int p;

    cData vterm, uprime, vprime, ploc, tv;
    int2 upwind;
    //float vxterm, vyterm;

    // gtidx is the domain location in x for this thread
    if (gtidx < dx)
    {
        for (p = 0; p < lb; p++)
        {
            // fi is the domain location in y for this thread
            int fi = gtidy + p;

            if (fi < dy)
            {

              // Read cell velocity
              ploc.x = (gtidx+0.5f)/(float)dx;// Why add 0.5f??
              ploc.y = (fi+0.5f)/(float)dy;
              vterm = tex2D(texref, ploc.x, ploc.y);
              upwind.x = vterm.x > 0 ? -1:1;
              upwind.y = vterm.y > 0 ? -1:1;

              // Forward difference, upwind
              ploc.x = (gtidx+0.5f+upwind.x)/(float)dx;
              ploc.y = (fi+0.5f)/(float)dy;
              uprime = tex2D(texref, ploc.x, ploc.y);
              uprime.x = (uprime.x - vterm.x) * (float)dx * upwind.x; // dudx
              uprime.y = (uprime.y - vterm.y) * (float)dx * upwind.x; // dvdx

              ploc.x = (gtidx+0.5f)/(float)dx;// Why add 0.5f??
              ploc.y = (fi+0.5f+upwind.y)/(float)dy;
              vprime = tex2D(texref, ploc.x, ploc.y);
              vprime.x = (vprime.x - vterm.x) * (float)dy * upwind.y; // dudy
              vprime.y = (vprime.y - vterm.y) * (float)dy * upwind.y; // dvdy

              tv.x = vterm.x - dt * (vterm.x * uprime.x + vterm.y * vprime.x);
              tv.y = vterm.y - dt * (vterm.x * uprime.y + vterm.y * vprime.y);

              //ploc.x -= dt * vterm.x;// Forward advection
              //ploc.y -= dt * vterm.y;
              //vterm = tex2D(texref, ploc.x, ploc.y);

              // Update to v
              cData *mv = (cData *)((char *)v + fi * pitch) + gtidx;
              *mv = tv;
            }
        }
    }
}

__global__ void
advectVelocityRev_k(cData *v, float *vx, float *vy,
                 int dx, int pdx, int dy, float dt, int lb, size_t pitch)
{

    int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
    int gtidy = blockIdx.y * (lb * blockDim.y) + threadIdx.y * lb;
    int p;

    cData vterm, tvterm, uprime, vprime, ploc, tv;
    int2 upwind;
    float vxterm, vyterm;

    // gtidx is the domain location in x for this thread
    if (gtidx < dx)
    {
        for (p = 0; p < lb; p++)
        {
            // fi is the domain location in y for this thread
            int fi = gtidy + p;

            if (fi < dy)
            {
                int fj = fi * pdx + gtidx;

                // Read cell velocity
                ploc.x = (gtidx+0.5f)/(float)dx;
                ploc.y = (fi+0.5f)/(float)dy;
                vterm = tex2D(tex_n_1_ref, ploc.x, ploc.y);
                //cData *mv = (cData *)((char *)v + fi * pitch) + gtidx;
                //vterm = v[fj];

                // Upwind direction
                tvterm = tex2D(texref, ploc.x, ploc.y);
                upwind.x = vterm.x > 0 ? -1:1;
                upwind.y = vterm.y > 0 ? -1:1;

                // Backward difference, upwind in reverse time direction
                ploc.x = (gtidx+0.5f-upwind.x)/(float)dx;
                ploc.y = (fi+0.5f)/(float)dy;
                uprime = tex2D(tex_n_1_ref, ploc.x, ploc.y);
                uprime.x = (vterm.x - uprime.x) * (float)dx * upwind.x; //dudx
                uprime.y = (vterm.y - uprime.y) * (float)dx * upwind.x;

                ploc.x = (gtidx+0.5f)/(float)dx;
                ploc.y = (fi+0.5f-upwind.y)/(float)dy;
                vprime = tex2D(tex_n_1_ref, ploc.x, ploc.y);
                vprime.x = (vterm.x - vprime.x) * (float)dy * upwind.y;
                vprime.y = (vterm.y - vprime.y) * (float)dy * upwind.y;

                tv.x = vterm.x + dt * (vterm.x * uprime.x + vterm.y * vprime.x);
                tv.y = vterm.y + dt * (vterm.x * uprime.y + vterm.y * vprime.y);

/*
                //ploc.x = (gtidx+0.5f)/(float)dx;// Why add 0.5f??
                //ploc.y = (fi+0.5f)/(float)dy;
                //vterm = tex2D(tex_n_1_ref, ploc.x, ploc.y);
                //ploc.x += dt * vterm.x;// Reverse advection
                //ploc.y += dt * vterm.y;
                //vterm = tex2D(tex_n_1_ref, ploc.x, ploc.y);

                // MacCormack
                tv.x = 0.5f * (tvterm.x + vterm.x)
                 - 0.5f * dt * (tvterm.x * uprime.x + tvterm.y * vprime.x);
                tv.y = 0.5f * (tvterm.y + vterm.y)
                 - 0.5f * dt * (tvterm.x * uprime.y + tvterm.y * vprime.y);

                vx[fj] = tv.x;
                vy[fj] = tv.y;
*/

                // Update phi_n - phi_n_hat to vx & vy
                vxterm = tv.x - tvterm.x; // phi_n - phi_n_hat
                vyterm = tv.y - tvterm.y;

                vterm.x = tvterm.x - 0.5f* vxterm; // u_wavebar
                vterm.y = tvterm.y - 0.5f* vyterm;
                // Write u_wavebar into v
                cData *mv = (cData *)((char *)v + fi * pitch) + gtidx;
                *mv = vterm;

                vx[fj] = vterm.x;
                vy[fj] = vterm.y;

            }
        }
    }
}

// Limited BFECC Scheme -> second order, unconditionally stable
__global__ void
advectVelocityBFECC_k(cData *v, float *vx, float *vy,
                 int dx, int pdx, int dy, float dt, int lb, size_t pitch)
{

    int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
    int gtidy = blockIdx.y * (lb * blockDim.y) + threadIdx.y * lb;
    int p;

    cData center, npos, uprime, vprime;
    float vxterm, vyterm;
    int2 upwind;
    // gtidx is the domain location in x for this thread
    if (gtidx < dx)
    {
        for (p = 0; p < lb; p++)
        {
            // fi is the domain location in y for this thread
            int fi = gtidy + p;

            if (fi < dy)
            {
                int fj = fi * pdx + gtidx;

                vxterm = vx[fj];
                vyterm = vy[fj];
                // Trace back along the initial characteristic - we'll use
                // values near this semi-Lagrangian "particle" to clamp our
                // final advected value.
                center.x = gtidx+0.5f;// Why add 0.5f??
                center.y = fi+0.5f;
                cData cellVelocity = tex2D(texref, center.x/(float)dx, center.y/(float)dy);

                npos.x = center.x - dt * cellVelocity.x * (float)dx;
                npos.y = center.y - dt * cellVelocity.y * (float)dy;
                // Find the cell corner closest to the "particle" and compute the
                // texture coordinate corresponding to that location.
                npos.x = floor(npos.x + 0.5f);
                npos.y = floor(npos.y + 0.5f);

                // Get the values of nodes that contribute to the interpolated value.
                // Texel centers will be a half-texel away from the cell corner.
                cData ht = make_float2(0.5f , 0.5f);
                cData nodeValues[4];
                nodeValues[0] = tex2D(texref, (npos.x-ht.x)/(float)dx, (npos.y-ht.y)/(float)dy);
                nodeValues[1] = tex2D(texref, (npos.x-ht.x)/(float)dx, (npos.y+ht.y)/(float)dy);
                nodeValues[2] = tex2D(texref, (npos.x+ht.x)/(float)dx, (npos.y-ht.y)/(float)dy);
                nodeValues[3] = tex2D(texref, (npos.x+ht.x)/(float)dx, (npos.y+ht.y)/(float)dy);

                // Determine a valid range for the result.
                cData phiMin, phiMax;
                phiMin.x = min(min(min(
                  nodeValues[0].x,  nodeValues[1].x), nodeValues[2].x), nodeValues[3].x);
                phiMin.y = min(min(min(
                  nodeValues[0].y,  nodeValues[1].y), nodeValues[2].y), nodeValues[3].y);
                phiMax.x = max(max(max(
                  nodeValues[0].x,  nodeValues[1].x), nodeValues[2].x), nodeValues[3].x);
                phiMax.y = max(max(max(
                  nodeValues[0].y,  nodeValues[1].y), nodeValues[2].y), nodeValues[3].y);
                // Perform final advection, combining values from intermediate
                // advection steps.
                //cData r_n_1 = tex2D(tex_n_1_ref, center.x, center.y);
                //cData *r_n_1 = (cData *)((char *)v + fi * pitch) + gtidx;
                //vxterm = r_n_1.x - 0.5f * vxterm;
                //vyterm = r_n_1.y - 0.5f * vyterm;

                // Read cell velocity
                //cellVelocity = tex2D(tex_n_1_ref, center.x, center.y);
                upwind.x = vxterm > 0 ? -1:1;
                upwind.y = vyterm > 0 ? -1:1;
                // Forward difference, upwind
                center.x = (gtidx+0.5f+upwind.x)/(float)dx;
                center.y = (fi+0.5f)/(float)dy;
                uprime = tex2D(tex_n_1_ref, center.x, center.y);
                uprime.x = (uprime.x - vxterm) * (float)dx * upwind.x; // dudx
                uprime.y = (uprime.y - vyterm) * (float)dx * upwind.x; // dvdx

                center.x = (gtidx+0.5f)/(float)dx;
                center.y = (fi+0.5f+upwind.y)/(float)dy;
                vprime = tex2D(tex_n_1_ref, center.x, center.y);
                vprime.x = (vprime.x - vxterm) * (float)dy * upwind.y; // dudy
                vprime.y = (vprime.y - vyterm) * (float)dy * upwind.y; // dvdy

                // Forward advection...
                cellVelocity.x = vxterm - dt * (vxterm * uprime.x + vyterm * vprime.x);
                cellVelocity.y = vyterm - dt * (vxterm * uprime.y + vyterm * vprime.y);
                vxterm = cellVelocity.x;
                vyterm = cellVelocity.y;

                // Clamp result to the desired range.
                vxterm = max(min(vxterm, phiMax.x), phiMin.x);
                vyterm = max(min(vyterm, phiMax.y), phiMin.y);

                // Bound by CFL
                //float CFL = 1.0f;
                //cData phiB = make_float2(CFL*(float)dx/dt, CFL*(float)dy/dt);
                //vxterm = max(min(vxterm, phiB.x), -phiB.x);
                //vyterm = max(min(vyterm, phiB.y), -phiB.y);

                vx[fj] = vxterm;
                vy[fj] = vyterm;
            }
        }
    }
}

// This method performs velocity diffusion and forces mass conservation
// in the frequency domain. The inputs 'vx' and 'vy' are complex-valued
// arrays holding the Fourier coefficients of the velocity field in
// X and Y. Diffusion in this space takes a simple form described as:
// v(k,t) = v(k,t) / (1 + visc * dt * k^2), where visc is the viscosity,
// and k is the wavenumber. The projection step forces the Fourier
// velocity vectors to be orthogonal to the vectors for each
// wavenumber: v(k,t) = v(k,t) - ((k dot v(k,t) * k) / k^2.
__global__ void
diffuseProject_k(cData *vx, cData *vy, int dx, int dy, float dt,
                 float visc, int lb)
{

    int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
    int gtidy = blockIdx.y * (lb * blockDim.y) + threadIdx.y * lb;
    int p;

    cData xterm, yterm;

    // gtidx is the domain location in x for this thread
    if (gtidx < dx)
    {
        for (p = 0; p < lb; p++)
        {
            // fi is the domain location in y for this thread
            int fi = gtidy + p;

            if (fi < dy)
            {
                int fj = fi * dx + gtidx;
                xterm = vx[fj];
                yterm = vy[fj];

                // Compute the index of the wavenumber based on the
                // data order produced by a standard NN FFT.
                int iix = gtidx;
                int iiy = (fi>dy/2)?(fi-(dy)):fi;
                //int iix = (gtidx>dx/2)?(gtidx - dx):gtidx;
                //int iiy = (fi>dy/2)?(fi-(dy)):fi;

                // Velocity diffusion
                float kk = (float)(iix * iix + iiy * iiy); // k^2
                float diff = 1.f / (1.f + visc * dt * kk);
                xterm.x *= diff;
                xterm.y *= diff;
                yterm.x *= diff;
                yterm.y *= diff;

                // Velocity projection
                if (kk > 0.f)
                {
                    float rkk = 1.f / kk;
                    // Real portion of velocity projection
                    float rkp = (iix * xterm.x + iiy * yterm.x);
                    // Imaginary portion of velocity projection
                    float ikp = (iix * xterm.y + iiy * yterm.y);
                    xterm.x -= rkk * rkp * iix;
                    xterm.y -= rkk * ikp * iix;
                    yterm.x -= rkk * rkp * iiy;
                    yterm.y -= rkk * ikp * iiy;
                }

                vx[fj] = xterm;
                vy[fj] = yterm;
            }
        }
    }
}

// This method updates the velocity field 'v' using the two complex
// arrays from the previous step: 'vx' and 'vy'. Here we scale the
// real components by 1/(dx*dy) to account for an unnormalized FFT.
__global__ void
updateVelocity_k(cData *v, float *vx, float *vy,
                 int dx, int pdx, int dy, int lb, size_t pitch)
{

    int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
    int gtidy = blockIdx.y * (lb * blockDim.y) + threadIdx.y * lb;
    int p;

    float vxterm, vyterm;
    cData nvterm;

    // gtidx is the domain location in x for this thread
    if (gtidx < dx)
    {
        for (p = 0; p < lb; p++)
        {
            // fi is the domain location in y for this thread
            int fi = gtidy + p;

            if (fi < dy)
            {
                int fjr = fi * pdx + gtidx;
                vxterm = vx[fjr];
                vyterm = vy[fjr];

                // Normalize the result of the inverse FFT
                float scale = 1.f / (dx * dy);
                nvterm.x = vxterm * scale;
                nvterm.y = vyterm * scale;

                vx[fjr] *= scale;
                vy[fjr] *= scale;
                cData *fj = (cData *)((char *)v + fi * pitch) + gtidx;
                *fj = nvterm;
            }
        } // If this thread is inside the domain in Y
    } // If this thread is inside the domain in X
}

// This method updates the particles by moving particle positions
// according to the velocity field and time step. That is, for each
// particle: p(t+1) = p(t) + dt * v(p(t)).
__global__ void
advectParticles_k(cData *part, cData *v, int dx, int dy,
                  float dt, int lb, size_t pitch)
{

    int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
    int gtidy = blockIdx.y * (lb * blockDim.y) + threadIdx.y * lb;
    int p;

    // gtidx is the domain location in x for this thread
    cData pterm, vterm;

    if (gtidx < dx)
    {
        for (p = 0; p < lb; p++)
        {
            // fi is the domain location in y for this thread
            int fi = gtidy + p;

            if (fi < dy)
            {
                int fj = fi * dx + gtidx;
                pterm = part[fj];

                int xvi = ((int)(pterm.x * dx)
              );
                int yvi = ((int)(pterm.y * dy));
                vterm = *((cData *)((char *)v + yvi * pitch) + xvi);

                pterm.x += dt * vterm.x;
                pterm.x = pterm.x - (int)pterm.x;
                pterm.x += 1.f;
                pterm.x = pterm.x - (int)pterm.x;
                pterm.y += dt * vterm.y;
                pterm.y = pterm.y - (int)pterm.y;
                pterm.y += 1.f;
                pterm.y = pterm.y - (int)pterm.y;

                part[fj] = pterm;
            }
        } // If this thread is inside the domain in Y
    } // If this thread is inside the domain in X
}

__global__ void
displayVelocityMag_k(cData *part, float4 *pcolor, cData *v, int dx, int dy,
                  float dt, int lb, size_t pitch)
{

    int gtidx = blockIdx.x * blockDim.x + threadIdx.x;
    int gtidy = blockIdx.y * (lb * blockDim.y) + threadIdx.y * lb;
    int p;

    // gtidx is the domain location in x for this thread
    cData pterm, vterm;
    float4 pcterm;

    if (gtidx < dx)
    {
        for (p = 0; p < lb; p++)
        {
            // fi is the domain location in y for this thread
            int fi = gtidy + p;

            if (fi < dy)
            {
                int fj = fi * dx + gtidx;
                pterm = part[fj];
                pcterm = pcolor[fj];

                // Get velocity
                int xvi = ((int)(pterm.x * dx));
                int yvi = ((int)(pterm.y * dy));
                vterm = *((cData *)((char *)v + yvi * pitch) + xvi);

                // Update color
                float CFL = 1.00f;
                cData phiB = make_float2(vterm.x/(CFL/(float)dx/dt),
                                         vterm.y/(CFL/(float)dy/dt));
                float umag = sqrtf(phiB.x*phiB.x + phiB.y*phiB.y) + 1E-7f;
                pcterm.x = 0.5f*sinf(40*umag        )+0.5f;
                pcterm.y = 0.5f*sinf(40*umag+2.0943f)+0.5f;
                pcterm.z = 0.5f*sinf(40*umag+4.1887f)+0.5f;
                pcterm.w = 1.0f;

                //float norm = 1.f/sqrtf(pcterm.x*pcterm.x + pcterm.y*pcterm.y + pcterm.z*pcterm.z);
                //pcterm.x *= norm;
                //pcterm.y *= norm;
                //pcterm.z *= norm;

                pcolor[fj] = pcterm;

            }
        } // If this thread is inside the domain in Y
    } // If this thread is inside the domain in X
}
// These are the external function calls necessary for launching fluid simulation
extern "C"
void addForces(cData *v, int dx, int dy, int spx, int spy, float fx, float fy, float dt, int r)
{

    dim3 tids(2*r+1, 2*r+1);

    addForces_k<<<1, tids>>>(v, dx, dy, spx, spy, fx, fy, dt, r, tPitch);
    getLastCudaError("addForces_k failed.");
}

extern "C"
void advectVelocity(cData *v, float *vx, float *vy, int dx, int pdx, int dy, float dt)
{
    dim3 grid((dx/TILEX)+(!(dx%TILEX)?0:1), (dy/TILEY)+(!(dy%TILEY)?0:1));
    dim3 tids(TIDSX, TIDSY);

    // Update texture n
    updateTexture(v, DIM*sizeof(cData), DIM, tPitch);
    // Forward advection
    advectVelocity_k<<<grid, tids>>>(v, vx, vy, dx, pdx, dy, dt, TILEY/TIDSY, tPitch);
    getLastCudaError("advectVelocity_k failed.");
    // Update texture n_1_hat
    updateTexture_n_1(v, DIM*sizeof(cData), DIM, tPitch);

    // Forward advection phi_n_1_hat -> phi_n_hat
    advectVelocityRev_k<<<grid, tids>>>(v, vx, vy, dx, pdx, dy, dt, TILEY/TIDSY, tPitch);
    getLastCudaError("advectVelocity_k failed.");
    // Update texture u_wavebar
    updateTexture_n_1(v, DIM*sizeof(cData), DIM, tPitch);

    advectVelocityBFECC_k<<<grid, tids>>>(v, vx, vy, dx, pdx, dy, dt, TILEY/TIDSY, tPitch);

}

extern "C"
void diffuseProject(cData *vx, cData *vy, int dx, int dy, float dt, float visc)
{
    // Forward FFT
    checkCudaErrors(cufftExecR2C(planr2c, (cufftReal *)vx, (cufftComplex *)vx));
    checkCudaErrors(cufftExecR2C(planr2c, (cufftReal *)vy, (cufftComplex *)vy));

    uint3 grid = make_uint3((dx/TILEX)+(!(dx%TILEX)?0:1),
                            (dy/TILEY)+(!(dy%TILEY)?0:1), 1);
    uint3 tids = make_uint3(TIDSX, TIDSY, 1);

    diffuseProject_k<<<grid, tids>>>(vx, vy, dx, dy, dt, visc, TILEY/TIDSY);
    getLastCudaError("diffuseProject_k failed.");

    // Inverse FFT
    checkCudaErrors(cufftExecC2R(planc2r, (cufftComplex *)vx, (cufftReal *)vx));
    checkCudaErrors(cufftExecC2R(planc2r, (cufftComplex *)vy, (cufftReal *)vy));
}

extern "C"
void updateVelocity(cData *v, float *vx, float *vy, int dx, int pdx, int dy)
{
    dim3 grid((dx/TILEX)+(!(dx%TILEX)?0:1), (dy/TILEY)+(!(dy%TILEY)?0:1));
    dim3 tids(TIDSX, TIDSY);

    updateVelocity_k<<<grid, tids>>>(v, vx, vy, dx, pdx, dy, TILEY/TIDSY, tPitch);
    getLastCudaError("updateVelocity_k failed.");
}

extern "C"
void advectParticles(GLuint vbo, cData *v, int dx, int dy, float dt)
{
    dim3 grid((dx/TILEX)+(!(dx%TILEX)?0:1), (dy/TILEY)+(!(dy%TILEY)?0:1));
    dim3 tids(TIDSX, TIDSY);

    cData *p;
    cudaGraphicsMapResources(1, &cuda_vbo_resource, 0);
    getLastCudaError("cudaGraphicsMapResources failed");

    size_t num_bytes;
    cudaGraphicsResourceGetMappedPointer((void **)&p, &num_bytes,
                                         cuda_vbo_resource);
    getLastCudaError("cudaGraphicsResourceGetMappedPointer failed");

    advectParticles_k<<<grid, tids>>>(p, v, dx, dy, dt, TILEY/TIDSY, tPitch);
    getLastCudaError("advectParticles_k failed.");

    cudaGraphicsUnmapResources(1, &cuda_vbo_resource, 0);
    getLastCudaError("cudaGraphicsUnmapResources failed");
}

extern "C"
void displayVelocityMag(GLuint vbo, GLuint vbo_color, cData *v, int dx, int dy, float dt)
{
    dim3 grid((dx/TILEX)+(!(dx%TILEX)?0:1), (dy/TILEY)+(!(dy%TILEY)?0:1));
    dim3 tids(TIDSX, TIDSY);

    cData *p;
    cudaGraphicsMapResources(1, &cuda_vbo_resource, 0);
    getLastCudaError("cudaGraphicsMapResources failed");

    float4 *pc;
    cudaGraphicsMapResources(1, &cuda_vbo_color_resource, 0);
    getLastCudaError("cudaGraphicsMapResources failed");

    size_t num_bytes;
    cudaGraphicsResourceGetMappedPointer((void **)&pc, &num_bytes,
                                         cuda_vbo_color_resource);

    getLastCudaError("cudaGraphicsResourceGetMappedPointer failed");

    cudaGraphicsResourceGetMappedPointer((void **)&p, &num_bytes,
                                         cuda_vbo_resource);

    getLastCudaError("cudaGraphicsResourceGetMappedPointer failed");

    displayVelocityMag_k<<<grid, tids>>>(p, pc, v, dx, dy, dt, TILEY/TIDSY, tPitch);
    getLastCudaError("advectParticles_k failed.");

    cudaGraphicsUnmapResources(1, &cuda_vbo_resource, 0);
    cudaGraphicsUnmapResources(1, &cuda_vbo_color_resource, 0);
    getLastCudaError("cudaGraphicsUnmapResources failed");
}
