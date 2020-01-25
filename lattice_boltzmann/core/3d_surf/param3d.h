#ifndef PARAM3D_H
#define PARAM3D_H

#include <iostream>
// Need another class RuntimeParam etc to expose some runtime configuration..
template<typename T>
class Lattice3DParam
{
public:
    T wi[28]; // 3*3=9,3*(3+1)=12, 3*9=27, 4*9=36 for uniform block padding
    T cx[28];
    T cy[28];
    T cz[28];
    uint nx, ny, nz, nelem;
    T csSqrInv;
    T csSqr;
    T Reg;
    T u_max;
    T tau;

    Lattice3DParam(uint w, uint h, uint d)
    {
        T w0 = 8.0/27.0; //zero weight
        T ws = 2.0/27.0; //adjacent weight
        T wd = 1.0/54.0; //diagonal weight
        T wc = 1.0/216.0;//tri-diagonal weight
        T wi_in[] = {ws, ws, ws, ws,wd, wd, wd, wd,    wd,wd,wd,wd,wc,wc,wc,wc,          wd,wd,wd,wd,wc,wc, wc, wc,    ws,ws,w0, 0};
        T cx_in[] = { 1, -1,  0,  0, 1, -1,  1, -1,     1, -1,  0,  0, 1, -1,  1, -1,    1, -1,  0,  0, 1, -1,  1,  -1, 0,  0,0, 0};
        T cy_in[] = { 0,  0,  1, -1, 1, -1, -1,  1,     0,  0,  1, -1, 1, -1, -1,  1,    0,  0,  1, -1, 1, -1, -1,   1, 0,  0,0, 0};
        T cz_in[] = { 0,  0,  0,  0, 0,  0,  0,  0,     1,  1,  1,  1, 1,  1,  1,  1,   -1, -1, -1, -1,-1, -1, -1,  -1, 1, -1,0, 0};
        std::copy(wi_in,wi_in+28,wi);
        std::copy(cx_in,cx_in+28,cx);
        std::copy(cy_in,cy_in+28,cy);
        std::copy(cz_in,cz_in+28,cz);

        csSqrInv = 3.0;
        csSqr = 1.0/3.0;

        u_max = 0.075;
        tau = 0.55;
        Reg = u_max/(tau - 0.5)*csSqrInv;

        nx = w;
        ny = h;
        nz = d;
        nelem = nx*ny*nz;

        // for acoustic simulations, tau -> 0.5 improves accuarcy (also introducing dispersive pattern)
        // for fluid simulations, tau > 0.6 improves stability
        std::cout << "tau = " << tau << std::endl;
        std::cout << "Reg = " << Reg << std::endl;
        std::cout << "Re = " << Reg * (ny-2) << std::endl;
    }
};

// Need another class RuntimeParam etc to expose some runtime configuration..
template<typename T>
class Lattice3DParamRuntime
{
public:
    uint ntask_i, ntask_j, ntask_k;
    const uint tile;
    uint timestep;

    Lattice3DParamRuntime(uint t, uint nbstep = 0):ntask_i(0),ntask_j(0),ntask_k(0),timestep(nbstep),tile(t){}
    bool internalIterationCounter(const uint& i, const uint& j, const uint& k, uint& nbstep)
    {
        bool callsForASwap = false;
        ntask_i++;
        if(ntask_i >= i) { ntask_i = 0; ntask_j++;}
        if(ntask_j >= j) { ntask_j = 0; ntask_k++;}
        if(ntask_k >= k) { ntask_k = 0; nbstep++;timestep = nbstep;callsForASwap = true;}
        return callsForASwap;
    }
    void reset(uint nbstep = 0)
    {
        ntask_i = 0;
        ntask_j = 0;
        ntask_k = 0;
        timestep = nbstep;
    }
    bool isInInternalCycle()
    {
        if(ntask_i != 0 || ntask_j != 0 || ntask_k != 0) return true;
        return false;
    }
};

#endif
