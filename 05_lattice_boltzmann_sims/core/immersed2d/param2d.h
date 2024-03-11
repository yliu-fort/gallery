#ifndef PARAM2D_H
#define PARAM2D_H

#include <iostream>

// Need another class RuntimeParam etc to expose some runtime configuration..
template<typename T>
class Lattice2DParam
{
public:
    T wi[12]; // 3*3=9,3*(3+1)=12, 3*9=27, 4*9=36 for uniform block padding
    T cx[12];
    T cy[12];
    uint nx, ny, nz, nelem;
    T csSqrInv;
    T csSqr;
    T Reg;
    T u_max;
    T tau;

    Lattice2DParam(uint w, uint h)
    {
        T w0 = 4.0/9.0; // zero weight
        T ws = 1.0/9.0; // adjacent weight
        T wd = 1.0/36.0; // diagonal weight
        T wi_in[] = {w0,ws,ws,   0,ws,ws,wd,     0,wd,wd,wd,      0};
        T cx_in[] = {0.0,1.0,0.0,0,-1.0, 0.0,1.0,0,-1.0,-1.0, 1.0,0};
        T cy_in[] = {0.0,0.0,1.0,0, 0.0,-1.0,1.0,0, 1.0,-1.0,-1.0,0};
        std::copy(wi_in,wi_in+12,wi);
        std::copy(cx_in,cx_in+12,cx);
        std::copy(cy_in,cy_in+12,cy);

        csSqrInv = 3.0;
        csSqr = 1.0/3.0;

        u_max = 0.075;
        tau = 0.505;
        Reg = u_max/(tau - 0.5)*csSqrInv;

        nx = w;
        ny = h;
        nz = 1;
        nelem = nx*ny*nz;

        // for acoustic simulations, tau -> 0.5 improves accuarcy (also introducing dispersive pattern)
        // for fluid simulations, tau > 0.6 improves stability
        std::cout << "Param2D:: tau = " << tau << std::endl;
        std::cout << "Param2D:: Reg = " << Reg << std::endl;
        std::cout << "Param2D:: Re = " << Reg * (ny-2) << std::endl;
    }
};

// Need another class RuntimeParam etc to expose some runtime configuration..
template<typename T>
class Lattice2DParamRuntime
{
public:
    uint ntask_i, ntask_j, ntask_k;
    const uint tile;
    uint timestep;

    Lattice2DParamRuntime(uint t, uint nbstep = 0):ntask_i(0),ntask_j(0),ntask_k(0),timestep(nbstep),tile(t){}
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
