#include <iostream>
#include <random>
#include <CL/sycl.hpp>

#include "param.h"

namespace algo
{

template<typename T>
void kinetic_simulation_configuration3D(std::vector<T>& buf, const Param& param)
{
    // Magic numbers
    //float pi = 3.141592654f;
    //float L = 10.0;
    //float eta = L/91.0f; // 991: unstable when t > 5
    //const uint nk = 64; // explicitly defined in kernels
    //float E0 = 3.2f;
    //float lambda = 0.5;

    //float LI = L/3.0f;
    //float TL = 0.2f*LI/sqrtf(E0);

    //srand(param.randomSeed);
    std::mt19937 rng(param.randomSeed);
    std::uniform_real_distribution<float> gen(0, 1);

    std::vector<double> theta,psi,alpha,k,E,dk,a;

    for(int i = 0; i < param.nk; i++){ theta.push_back(2.0*param.pi*gen(rng));/* skip alternate values */rng(); }
    for(int i = 0; i < param.nk; i++){ psi.push_back(2.0*param.pi*gen(rng));/* skip alternate values */rng(); }

    double k1 = 2.0*param.pi/param.L;
    //float k_Nk = 2*pi/eta;

    //float alpha = (L/eta)^(1/(nk-1));
    for(int i = 0; i < param.nk; i++){ alpha.push_back(pow(param.L/param.eta,1.0/(param.nk-1.0))); }
    //float k = k1*alpha.^((1:nk).'-1);
    for(int i = 0; i < param.nk; i++){ k.push_back(k1*pow(alpha[i],(double)i)); }
    //float E = E0*L*(k*L).^(-5/3);
    for(int i = 0; i < param.nk; i++){ E.push_back(param.E0*param.L*pow((k[i]*param.L),-param.energy_cascade)); }
    //float dk = gradient(k);
    for(int i = 0; i < param.nk; i++){
        if(i == 0){ dk.push_back(k[1]-k[0]);continue; }
        if(i == param.nk-1){ dk.push_back(k[param.nk-1]-k[param.nk-2]);break; }
        dk.push_back(0.5*(k[i+1]-k[i-1]));
    }

    //float a = sqrt(E.*dk);
    for(int i = 0; i < param.nk; i++){ a.push_back(sqrt(E[i]*dk[i])); }
    //float b = a;

    //float omega = lambda*sqrt(k.^3.*E).';
    // Convert from double to float in order to squeeze into the buffer
    std::vector<double> a_xt, a_yt, a_zt, k_xt, k_yt, k_zt, omega;
    for(int i = 0; i < param.nk; i++)
    {
        omega.push_back( param.lambda*sqrt(E[i]*pow(k[i],3.0)) );
    }

    // Vectors
    //An = a.*[ cos(theta).*cos(psi) -sin(theta)  cos(theta).*sin(psi)];
    //Bn = b.*[-cos(theta).*cos(psi)  sin(theta) -cos(theta).*sin(psi)];
    //Kn = k.*[ sin(theta).*cos(psi)  cos(theta)  sin(theta).*sin(psi)];

    //std::vector<float> b_xt, b_yt, b_zt;
    for(int i = 0; i < param.nk; i++){
        a_xt.push_back( a[i]*cos(theta[i])*cos(psi[i]) );
        a_yt.push_back(-a[i]*sin(theta[i]) );
        a_zt.push_back( a[i]*cos(theta[i])*sin(psi[i]) );

        k_xt.push_back( k[i]*sin(theta[i])*cos(psi[i]) );
        k_yt.push_back( k[i]*cos(theta[i]) );
        k_zt.push_back( k[i]*sin(theta[i])*sin(psi[i]) );
    }

    buf.resize(7*param.nk);
    //assert(param.nk % 4 == 0, "Number of modes must be a multiply of 4 for alignment purpose!\n");
    for(int i = 0; i < param.nk; i++)
    {
        buf[i+0*param.nk] = a_xt[i];
        buf[i+1*param.nk] = a_yt[i];
        buf[i+2*param.nk] = a_zt[i];

        buf[i+3*param.nk] = k_xt[i];
        buf[i+4*param.nk] = k_yt[i];
        buf[i+5*param.nk] = k_zt[i];

        buf[i+6*param.nk] = omega[i];

    }

    // push parameters
    param.push_params(buf);
}

template<typename T, int N,
         cl::sycl::access::target DATA_STORAGE_TYPE = cl::sycl::access::target::global_buffer,
         cl::sycl::access::target CONST_STORAGE_TYPE = cl::sycl::access::target::constant_buffer
         >
class KineticTurb3D
{
public:
    using Vec = cl::sycl::vec<T,N>;
    using vec4 = cl::sycl::vec<T,4>;

    using acc_t = cl::sycl::accessor
    <Vec,1,
    cl::sycl::access::mode::read_write,
    DATA_STORAGE_TYPE>;

    using cst_t = cl::sycl::accessor
    <vec4,1,
    cl::sycl::access::mode::read,
    CONST_STORAGE_TYPE>;

public:
    KineticTurb3D(acc_t p, cst_t c, bool have_event = true) :
        ptr(p)
      , param(c)
      , haveEventFcn(have_event)
    {}

    void operator()( const Vec &x , Vec &dxdt , const T t ) const
    {
        T tau_inv = (T)param[7*chunk].y();
        T gtau = (T)param[7*chunk].z();
        T stride = (T)param[7*chunk].w();

        vec4 uMat = vec4(0);
        vec4 vMat = vec4(0);
        vec4 wMat = vec4(0);

        // A large L-vortex appeared near origin?
        // Add some constant in position to avoid that...
        T px = x.s0() + stride; // note: need a random machine...
        T py = x.s1() + stride;
        T pz = x.s2() + stride;

        // Unrolled loops
        for(int i = 0; i < chunk; i++)
        {
            vec4 KS,dKS;

            KS = px*param[i+3*chunk] + py*param[i+4*chunk] + pz*param[i+5*chunk] + param[i+6*chunk]*(t + stride);
            dKS = cl::sycl::cos( KS ) - cl::sycl::sin( KS );

            uMat += param[i+0*chunk]*dKS;
            vMat += param[i+1*chunk]*dKS;
            wMat += param[i+2*chunk]*dKS;

        }

        T vx = x.s4();
        T vy = x.s5();
        T vz = x.s6();

        // note: should have better approach to read parameters
        dxdt = Vec(vx, vy, vz, 0.0, tau_inv*(dot(uMat,vec4(1)) - vx) , tau_inv*(dot(vMat,vec4(1)) - vy  + gtau), tau_inv*(dot(wMat,vec4(1)) - vz), 0.0);
        //dxdt = Vec(vec4(x.hi()), tau_inv*(vec4(dot(uMat,vec4(1)),dot(vMat,vec4(1)),dot(wMat,vec4(1)),0) - x.hi()));

    }

    bool event ( const Vec &x , const Vec &dxdt , const T t ) const
    {
        bool trigger = (T)x.s1() < (T)0.0 && (T)dxdt.s1() < (T)0; // passing through x = 1.0 from bottom only

        // when a event is triggered, program will refine the timestep until convergence
        return trigger;
    }

    void output ( const Vec &x , const Vec &dxdt , const T t ) const
    {
        //out << "debug ostream" << sycl::endl;
    }

    // ODEoptions
    void setODEoption(T& rtol, T& atol, T& max_dt, bool& have_event_fcn, bool& have_output_fcn) const
    {
        rtol = param[7*chunk+1].x(); // less than 1e-8 would be very laggy
        atol = param[7*chunk+1].y();
        max_dt = cl::sycl::min((T)0.1, param[7*chunk+1].z());
        have_event_fcn = haveEventFcn;
        have_output_fcn = haveOutputFcn;
    }


private:

    bool haveEventFcn;
    bool haveOutputFcn = false;

    // Kinetic simu
    int chunk = 16;

    acc_t ptr;
    cst_t param;
};

template<typename T, int N,
         cl::sycl::access::target DATA_STORAGE_TYPE = cl::sycl::access::target::global_buffer,
         cl::sycl::access::target CONST_STORAGE_TYPE = cl::sycl::access::target::constant_buffer
         >
class KineticTurb3DPassive
{
public:
    using Vec = cl::sycl::vec<T,N>;
    using vec4 = cl::sycl::vec<T,4>;

    using acc_t = cl::sycl::accessor
    <Vec,1,
    cl::sycl::access::mode::read_write,
    DATA_STORAGE_TYPE>;

    using cst_t = cl::sycl::accessor
    <vec4,1,
    cl::sycl::access::mode::read,
    CONST_STORAGE_TYPE>;

public:
    KineticTurb3DPassive(acc_t p, cst_t c) :
        ptr(p)
      , param(c)
    {}

    void operator()( const Vec &x , Vec &dxdt , const T t ) const
    {
        //T tau_inv = (T)param[7*chunk].y();
        //T gtau = (T)param[7*chunk].z();
        T stride = (T)param[7*chunk].w();

        vec4 uMat = vec4(0);
        vec4 vMat = vec4(0);
        vec4 wMat = vec4(0);

        // A large L-vortex appeared near origin?
        // Add some constant in position to avoid that...
        T px = x.s0() + stride; // note: need a random machine...
        T py = x.s1() + stride;
        T pz = x.s2() + stride;

        // Unrolled loops
        for(int i = 0; i < chunk; i++)
        {
            vec4 KS,dKS;

            KS = px*param[i+3*chunk] + py*param[i+4*chunk] + pz*param[i+5*chunk] + param[i+6*chunk]*(t + stride);
            dKS = cl::sycl::cos( KS ) - cl::sycl::sin( KS );

            uMat += param[i+0*chunk]*dKS;
            vMat += param[i+1*chunk]*dKS;
            wMat += param[i+2*chunk]*dKS;

        }

        // note: should have better approach to read parameters
        dxdt = Vec(dot(uMat,vec4(1)),
                   dot(vMat,vec4(1)),
                   dot(wMat,vec4(1)), 0.0,0.0,0.0,0.0,0.0);
        //dxdt = Vec(vec4(x.hi()), tau_inv*(vec4(dot(uMat,vec4(1)),dot(vMat,vec4(1)),dot(wMat,vec4(1)),0) - x.hi()));

    }

    bool event ( const Vec &x , const Vec &dxdt , const T t ) const
    {
        return false;
    }

    void output ( const Vec &x , const Vec &dxdt , const T t ) const
    {}

    // ODEoptions
    void setODEoption(T& rtol, T& atol, T& max_dt, bool& have_event_fcn, bool& have_output_fcn) const
    {
        rtol = param[7*chunk+1].x(); // less than 1e-8 would be very laggy
        atol = param[7*chunk+1].y();
        max_dt = cl::sycl::min((T)0.1, param[7*chunk+1].z());
        have_event_fcn = haveEventFcn;
        have_output_fcn = haveOutputFcn;
    }


private:

    bool haveEventFcn = false;
    bool haveOutputFcn = false;

    // Kinetic simu
    int chunk = 16;

    acc_t ptr;
    cst_t param;
};
}
