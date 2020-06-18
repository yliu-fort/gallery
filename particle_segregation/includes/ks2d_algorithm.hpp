#include <iostream>
#include <random>
#include <CL/sycl.hpp>

#include "param.h"

namespace algo
{

template<typename T>
void kinetic_simulation_configuration2D(std::vector<T>& buf, const Param& param)
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

    // Do not remove extra rng()s' (drop alternative value per call help reproduce the behaviour of matlab rand())
    for(int i = 0; i < param.nk; i++){ theta.push_back(2.0*param.pi*gen(rng));/* skip alternate values */rng(); }
    for(int i = 0; i < param.nk; i++){ psi.push_back(2.0*param.pi*gen(rng));/* skip alternate values */rng(); }

    double k1 = 2.0*param.pi/param.L;
    //float k_Nk = 2*pi/eta;

    //float alpha = (L/eta)^(1/(nk-1));
    for(int i = 0; i < param.nk; i++){ alpha.push_back(pow(param.L/param.eta,1.0/(param.nk-1.0))); }
    //float k = k1*alpha.^((1:nk).'-1);
    for(int i = 0; i < param.nk; i++){ k.push_back(k1*pow(alpha[i],(double)i));}

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
    std::vector<double> a_xt, a_yt, k_xt, k_yt, omega;
    for(int i = 0; i < param.nk; i++)
    {
        omega.push_back( param.lambda*sqrt(E[i]*pow(k[i],3.0)) );
    }

    // Vectors
    //An = (a.*[ cos(theta) -sin(theta)]).';
    //Bn = (b.*[-cos(theta)  sin(theta)]).';
    //Kn = (k.*[ sin(theta)  cos(theta)]).';
    //std::vector<float> b_xt, b_yt, b_zt;
    for(int i = 0; i < param.nk; i++){
        a_xt.push_back( a[i]*cos(theta[i]) );
        a_yt.push_back(-a[i]*sin(theta[i]) );

        k_xt.push_back( k[i]*sin(theta[i]) );
        k_yt.push_back( k[i]*cos(theta[i]) );
    }

    buf.resize(5*param.nk);
    //assert(param.nk % 4 == 0, "Number of modes must be a multiply of 4 for alignment purpose!\n");
    for(int i = 0; i < param.nk; i++)
    {
        buf[i+0*param.nk] = a_xt[i];
        buf[i+1*param.nk] = a_yt[i];
        buf[i+2*param.nk] = k_xt[i];
        buf[i+3*param.nk] = k_yt[i];
        buf[i+4*param.nk] = omega[i];

    }

    // push parameters
    param.push_params(buf);

}

template<typename T, int N,
         cl::sycl::access::target DATA_STORAGE_TYPE = cl::sycl::access::target::global_buffer,
         cl::sycl::access::target CONST_STORAGE_TYPE = cl::sycl::access::target::constant_buffer
         >
class KineticTurb2D
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

    KineticTurb2D(acc_t p, cst_t c, bool have_event = true, bool have_output = false) :
        ptr(p)
      , param(c)
      , haveEventFcn(have_event)
      , haveOutputFcn(have_output)
    {}

    void operator()( const Vec &x , Vec &dxdt , const T t ) const
    {
        T tau_inv = (T)param[5*chunk].y();
        T gtau = (T)param[5*chunk].z();
        T stride = (T)param[5*chunk].w();

        Vec uMat = Vec(0);
        Vec vMat = Vec(0);

        // A large L-vortex appeared near origin?
        // Add some constant in position to avoid that...
        T px = x.x() + stride;
        T py = x.y() + stride;

        // Unrolled loops
        for(int i = 0; i < chunk; i++)
        {
            Vec KS,dKS;

            KS = px*param[i+2*chunk] + py*param[i+3*chunk] + param[i+4*chunk]*(t + stride);
            dKS = cl::sycl::cos( KS ) - cl::sycl::sin( KS );

            uMat += param[i+0*chunk]*dKS;
            vMat += param[i+1*chunk]*dKS;

        }

        T vx = x.z();
        T vy = x.w();

        // note: should have better approach to read parameters
        dxdt = Vec( vx, vy, tau_inv*(dot(uMat,Vec(1)) - vx) , tau_inv*(dot(vMat,Vec(1)) - vy + gtau) );

    }

    bool event ( const Vec &x , const Vec &dxdt , const T t ) const
    {
        bool trigger = (T)x.y() < (T)0.0 && (T)dxdt.y() < (T)0; // passing through x = 1.0 from bottom only

        // when a event is triggered, program will refine the timestep until convergence
        return trigger;
    }

    void output ( const Vec &x , const Vec &dxdt , const T t ) const
    {
        //out << "debug ostream" << sycl::endl;
    }

    void setODEoption(T& rtol, T& atol, T& max_dt, bool& have_event_fcn, bool& have_output_fcn) const
    {
        rtol = param[5*chunk+1].x(); // less than 1e-8 would be very laggy
        atol = param[5*chunk+1].y();
        max_dt = cl::sycl::min((T)0.1, param[5*chunk+1].z());
        have_event_fcn = haveEventFcn;
        have_output_fcn = haveOutputFcn;
    };

private:

    bool haveEventFcn;
    bool haveOutputFcn;

    // Kinetic simu
    const int chunk = 16;

    acc_t ptr;
    cst_t param;
};

template<typename T, int N,
         cl::sycl::access::target DATA_STORAGE_TYPE = cl::sycl::access::target::global_buffer,
         cl::sycl::access::target CONST_STORAGE_TYPE = cl::sycl::access::target::constant_buffer
         >
class KineticTurb2DPassive
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

    KineticTurb2DPassive(acc_t p, cst_t c) :
        ptr(p)
      , param(c)
    {}

    void operator()( const Vec &x , Vec &dxdt , const T t ) const
    {
        //T tau_inv = (T)param[5*chunk].y();
        //T gtau = (T)param[5*chunk].z();
        T stride = (T)param[5*chunk].w();

        Vec uMat = Vec(0);
        Vec vMat = Vec(0);

        // A large L-vortex appeared near origin?
        // Add some constant in position to avoid that...
        T px = x.x() + stride; // note: need a random machine...
        T py = x.y() + stride;

        // Unrolled loops
        for(int i = 0; i < chunk; i++)
        {
            Vec KS,dKS;

            KS = px*param[i+2*chunk] + py*param[i+3*chunk] + param[i+4*chunk]*(t + stride);
            dKS = cl::sycl::cos( KS ) - cl::sycl::sin( KS );

            uMat += param[i+0*chunk]*dKS;
            vMat += param[i+1*chunk]*dKS;

        }

        // note: should have better approach to read parameters
        dxdt = Vec( dot(uMat,Vec(1)), dot(vMat,Vec(1)), 0, 0 );

    }

    bool event ( const Vec &x , const Vec &dxdt , const T t ) const
    {
        return false;
    }

    void output ( const Vec &x , const Vec &dxdt , const T t ) const
    {
    }

    void setODEoption(T& rtol, T& atol, T& max_dt, bool& have_event_fcn, bool& have_output_fcn) const
    {
        rtol = param[5*chunk+1].x(); // less than 1e-8 would be very laggy
        atol = param[5*chunk+1].y();
        max_dt = cl::sycl::min((T)0.1, param[5*chunk+1].z());
        have_event_fcn = haveEventFcn;
        have_output_fcn = haveOutputFcn;
    };

private:

    bool haveEventFcn = false;
    bool haveOutputFcn = false;

    // Kinetic simu
    const int chunk = 16;

    acc_t ptr;
    cst_t param;

};

}

//// Current SYCL does not support virtual function in device code
/*template<typename T, int N,
         cl::sycl::access::target DATA_STORAGE_TYPE = cl::sycl::access::target::global_buffer,
         cl::sycl::access::target CONST_STORAGE_TYPE = cl::sycl::access::target::constant_buffer
         >
class KineticTurb2DBase
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


    KineticTurb2DBase()
    {}

    virtual void operator()( const Vec &x , Vec &dxdt , const T t ) const = 0;

    virtual bool event ( const Vec &x , const Vec &dxdt , const T t ) const { return false };

    virtual void output ( const Vec &x , const Vec &dxdt , const T t ) const {};

    virtual void setODEoption(T& rtol, T& atol, T& max_dt, bool& have_event_fcn, bool& have_output_fcn) const = 0;

};*/
