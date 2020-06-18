#include <iostream>
#include "syclode45.hpp"
#include "sycl_integrator_rk45.hpp"
#include "sycl_integrator_adm4.hpp"
#include "param.h"
#include "ks2d_algorithm.hpp"

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
    if(!imported_param.config_param(config_file)) { return false; }
    imported_param.print_param();
    array_size = imported_param.dim[0];

    T L = (T)imported_param.L;
    T fdist = (T)imported_param.b*(T)imported_param.L; // falling distance
    T grav_acc = (T)imported_param.g;
    T tau_inv = imported_param.tau_inv;
    T E0 = imported_param.E0;

    // custom setting
    //array_size = 1<<20;
    //fdist *= 0.25; // 9.81 by default
    //grav_acc *= 0.5;
    //tau_inv *= 2;

    int i = 0;
    std::vector<vec> A(array_size);
    std::generate(A.begin(), A.end(), [&](){ return vec(-0.5*L + (++i)*L/(T)array_size,fdist,0,0); });
    buf = buffer<vec, 1>(A.begin(), A.end());//.template reinterpret<vec>(range<1>(array_size));

    buf_s = buffer<T, 1>(range<1>(array_size));

    std::vector<T> B;
    algo::kinetic_simulation_configuration2D(B, imported_param);
    param = buffer<T, 1>(B.begin(), B.end()).template reinterpret<cl::sycl::vec<T,4>>(range<1>(B.size()/4));

    return true;
}

template <typename T>
class printkernel;

template <typename T, int N>
void ode45<T, N>::calc(float tf)
{

    T init_dt = 0.01*imported_param.get_kolmogorov_time();

    submit([&](handler& cgh) {
        auto val_t = buf.template get_access<sycl_rw>(cgh);
        auto time_t = buf_s.template get_access<sycl_rw>(cgh);
        auto param_t = param.template get_access<sycl_read, access::target::constant_buffer>(cgh);

        //sycl::stream out(1024, 256, cgh);

        algo::KineticTurb2D<T,N> odeHandle(val_t, param_t);

        auto kern = [=](id<1> idx) {

            auto x = val_t[idx];
            auto t = time_t[idx];

            algo::integrator_rk45(odeHandle, t, (T)tf, x, init_dt);
            //algo::integrator_adm4(odeHandle, t, tf, x, init_dt);

            val_t[idx] = x;
            time_t[idx] = t;
            //out << x << ", " << t << sycl::endl;

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

    int stride = 4+1;
    std::vector<T> out(stride*array_size);
    for(int i = 0; i < array_size; i++)
    {
        out[stride*i+0] = h_data[i].x();
        out[stride*i+1] = h_data[i].y();
        out[stride*i+2] = h_data[i].z();
        out[stride*i+3] = h_data[i].w();
        out[stride*i+4] = t_data[i];
    }
    imported_param.dump(out, stride, filename);
}

// explicit initiate class templates
template class ode45<float,4>;
template class ode45<double,4>;
