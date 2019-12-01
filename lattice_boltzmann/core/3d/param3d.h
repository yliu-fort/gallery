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
    T csSqrInv;
    T csSqr;
    T Reg;
    T u_max;
    T tau;
    int nx, ny, nz, nelem;

    Lattice3DParam(int w, int h, int d)
    {
        T w0 = 8.0/27.0; //zero weight
        T ws = 2.0/27.0; //adjacent weight
        T wd = 1.0/54.0; //diagonal weight
        T wc = 1.0/216.0;//tri-diagonal weight
        // Arrays of the lattice weights and direction components
        T wi_in[] = {w0, ws, ws, ws, ws, ws, ws, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wd, wc, wc, wc, wc, wc, wc, wc, wc, 0};
        T cx_in[] = {0, 1, -1,  0,  0,  0,  0,  1, -1,  1, -1, 1, -1,  1, -1, 0,  0,  0,  0, 1, -1,  1, -1, -1,  1, -1,  1, 0};
        T cy_in[] = {0, 0,  0,  1, -1,  0,  0,  1, -1, -1,  1, 0,  0,  0,  0, 1, -1,  1, -1, 1, -1, -1,  1,  1, -1, -1,  1, 0};
        T cz_in[] = {0, 0,  0,  0,  0,  1, -1,  0,  0,  0,  0, 1, -1, -1,  1, 1, -1, -1,  1, 1, -1,  1, -1,  1, -1,  1, -1, 0};
        std::copy(wi_in,wi_in+28,wi);
        std::copy(cx_in,cx_in+28,cx);
        std::copy(cy_in,cy_in+28,cy);
        std::copy(cz_in,cz_in+28,cz);
        csSqrInv = 3.0;
        csSqr = 1.0/3.0;

        u_max = 0.075;
        tau = 0.505;
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

#endif
