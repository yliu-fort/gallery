#ifndef PARAM2D_H
#define PARAM2D_H

#include <iostream>

// Need another class RuntimeParam etc to expose some runtime configuration..
template<typename T>
class Lattice2DParam
{
public:
    T wi[12]; // 3*4 for uniform block padding
    T cx[12];
    T cy[12];
    T csSqrInv;
    T csSqr;
    T Reg;
    T u_max;
    T tau;
    int nx, ny;

    Lattice2DParam(int w, int h)
    {
        T w0 = 4.0/9.0; // zero weight
        T ws = 1.0/9.0; // adjacent weight
        T wd = 1.0/36.0; // diagonal weight
        T wi_in[] = {w0,ws,ws,0,ws,ws,wd,0,wd,wd,wd,0};
        T cx_in[] = {0.0,1.0,0.0,0,-1.0, 0.0,1.0,0,-1.0,-1.0, 1.0,0};
        T cy_in[] = {0.0,0.0,1.0,0, 0.0,-1.0,1.0,0, 1.0,-1.0,-1.0,0};
        std::copy(wi_in,wi_in+12,wi);
        std::copy(cx_in,cx_in+12,cx);
        std::copy(cy_in,cy_in+12,cy);
        csSqrInv = 3.0;
        csSqr = 1.0/3.0;

        u_max = 0.075;
        tau = 0.51;
        Reg = u_max/(tau - 0.5)*csSqrInv;

        nx = w;
        ny = h;

        // for acoustic simulations, tau -> 0.5 improves accuarcy (also introducing dispersive pattern)
        // for fluid simulations, tau > 0.6 improves stability
        std::cout << "tau = " << tau << std::endl;
        std::cout << "Reg = " << Reg << std::endl;
        std::cout << "Re = " << Reg * (ny-2) << std::endl;
    }
};

#endif
