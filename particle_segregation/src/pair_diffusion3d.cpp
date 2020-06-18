#include <iostream>
#include <string>

#include "syclode45.hpp"
#define ABSOLUTE_ERR_CONTROL
#include "sycl_integrator_rk45.hpp"
#include "sycl_integrator_adm4.hpp"
#include "sycl_integrator_rkf45.hpp"
#include "param.h"
#include "ks3d_algorithm.hpp"

#define HAVE_EVENT_FUNC false

using namespace cl::sycl;

// state_type = double
constexpr access::mode sycl_read = access::mode::read;
constexpr access::mode sycl_write = access::mode::write;
constexpr access::mode sycl_rw = access::mode::read_write;

/* This is the class used to name the kernel for the runtime.
 * This must be done when the kernel is expressed as a lambda. */
template <typename T, int N>
class _lambda_ode45;


template <typename T, int N>
bool ode45<T, N>::initialize(const char* config_file)
{

    // read parameters
    array_size = imported_param.dim[0];

    // set/get params
    T L = (T)imported_param.L;
    T fdist = (T)imported_param.b*(T)imported_param.L; // falling distance
    T zdist = (T)imported_param.b*(T)imported_param.eta; // falling distance
    T gtau = (T)imported_param.get_drift_velocity();
    imported_param.initDist = imported_param.eta;

    // RAII scalar buffer (time)
    buf_s = buffer<T, 1>(range<1>(array_size));

    // init constant buffer
    std::vector<T> B;
    algo::kinetic_simulation_configuration3D(B, imported_param);
    param = buffer<T, 1>(B.begin(), B.end()).template reinterpret<cl::sycl::vec<T,4>>(range<1>(B.size()/4));

    // Debug write param to disk
    //imported_param.dump(B, 1, "param.dat.1");
    //std::cout << "param has been dumped\n" << std::endl;

    // gen vel
    auto uf = [&](const T& x, const T& y, const T& z, const T t, const int dim)
    {

        // refer parameters
        const auto& buf = B;
        const auto& param = imported_param;

        // shifted pos
        T px = x + param.xoffset; // note: need a random machine...
        T py = y + param.xoffset;
        T pz = z + param.xoffset;

        T vel[3]{0};

        // compute ks parameters
        for(unsigned int i = 0; i < param.nk; i++)
        {
            T KS,dKS;

            KS =      px*buf[i+3*param.nk]
                    + py*buf[i+4*param.nk]
                    + pz*buf[i+5*param.nk]
                    + buf[i+6*param.nk]*(t + param.xoffset);
            dKS = cos( KS ) - sin( KS );

            vel[0] += buf[i+0*param.nk]*dKS;
            vel[1] += buf[i+1*param.nk]*dKS;
            vel[2] += buf[i+2*param.nk]*dKS;

        }

        // add external force
        //vel[1] += gtau;

        return vel[dim];
    };

    // init state vector [pos vel]
    int i = 0;
    std::vector<vec> A(array_size);

    // fill buffer
    switch (imported_param.initApproach) {
    case algo::ZERO_ACC:
    {
        std::generate(A.begin(), A.end(), [&](){
            T px = (++i)*imported_param.initDist;
            T py = 0.0;
            T pz = zdist;
            return vec(px,py,pz,0,uf(px,py,pz,0,0),uf(px,py,pz,0,1) + gtau,uf(px,py,pz,0,2),0);
        });
        break;
    }
    case algo::ZERO_VEL:
    {
        std::generate(A.begin(), A.end(), [&](){
            T px = (++i)*imported_param.initDist;
            T py = 0.0;
            T pz = zdist;
            return vec(px,py,pz,0,0,0,0,0);
        });
        break;
    }
    case algo::FLUID_VEL:
    {
        std::generate(A.begin(), A.end(), [&](){
            T px = (++i)*imported_param.initDist;
            T py = 0.0;
            T pz = zdist;
            return vec(px,py,pz,0,uf(px,py,pz,0,0),uf(px,py,pz,0,1),uf(px,py,pz,0,2),0);
        });
        break;
    }
    default:
        break;
    }

    // transfer buffer from host to device
    buf = buffer<vec, 1>(A.begin(), A.end());

    return true;
}

template <typename T>
class printkernel;

template <typename T, int N>
void ode45<T, N>::calc(float tf)
{

    // get max dt
    tf *= imported_param.get_kolmogorov_time();
    T init_dt = imported_param.get_max_allowed_dt();
    bool isPassive = imported_param.is_passive();

    // submit job
    submit([&](handler& cgh) {
        auto val_t = buf.template get_access<sycl_rw>(cgh);
        auto time_t = buf_s.template get_access<sycl_rw>(cgh);
        auto param_t = param.template get_access<sycl_read, access::target::constant_buffer>(cgh);

        algo::KineticTurb3DPassive<T,N> func_passive(val_t, param_t);
        algo::KineticTurb3D<T,N> func_inertial(val_t, param_t, false);

        auto kern = [=](id<1> idx) {

            auto x = val_t[idx];
            auto t = time_t[idx];

            if(isPassive)
                algo::integrator_rk45(func_passive, t, (T)tf, x, init_dt);
            else
                algo::integrator_rk45(func_inertial, t, (T)tf, x, init_dt);


            val_t[idx] = x;
            time_t[idx] = t;

        };
        cgh.parallel_for<class _lambda_ode45<T, N>>(range<1>(array_size), kern);
    });

}

template <typename T, int N>
void ode45<T, N>::finalize(const char* filename)
{
    // query host pointer will hit the performance heavily. only do it when debugging
    wait_and_throw();
    auto h_data = buf.template get_access<sycl_read>();
    auto t_data = buf_s.template get_access<sycl_read>();

    // write to host mem
    int stride = 6+1;
    std::vector<T> out(stride*array_size);
    for(int i = 0; i < array_size; i++)
    {
        out[stride*i+0] = h_data[i].s0();
        out[stride*i+1] = h_data[i].s1();
        out[stride*i+2] = h_data[i].s2();
        out[stride*i+3] = h_data[i].s4();
        out[stride*i+4] = h_data[i].s5();
        out[stride*i+5] = h_data[i].s6();
        out[stride*i+6] = t_data[i];
    }

    // output to disk
    imported_param.dump(out, stride, filename);
}

// explicit initiate class templates
template class ode45<float,8>;
template class ode45<double,8>;
