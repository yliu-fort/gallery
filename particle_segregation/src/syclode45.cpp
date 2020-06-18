#include "syclode45.hpp"

#include <iostream>

#include "sycl_integrator_rk45.hpp"
#include "sycl_ode_system.hpp"

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
    array_size = 16;
    std::vector<T> A(N*array_size);
    srand(2393425);
    auto gen = [](){ return (2.0*rand()/(double(RAND_MAX))-1.0); };
    std::generate(A.begin(), A.end(), [=](){ return gen(); });
    buf = buffer<T, 1>(A.begin(), A.end()).template reinterpret<vec>(range<1>(array_size));

    return true;

}

template <typename T, int N>
void ode45<T, N>::calc(float tf)
{
    // pass 1: host run
    //T x = 0;
    //algo::integrator_rk45(rhs, 1.0, 10.0, x);
    //std::cout << x << std::endl;

    // pass 2: kernel run
    submit([&](handler& cgh) {
        auto accessorA = buf.template get_access<sycl_rw>(cgh);

        algo::ODEsystem<T,N> odeHandle(accessorA);

        auto kern = [=](id<1> wiID) {
            auto x = accessorA[wiID];
            T t = 1;

            algo::integrator_rk45(odeHandle, t, (T)tf, x);

            accessorA[wiID] = x;
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
    for(int i = 0; i < array_size; i++)
    {
        T s[4];
        s[0] = h_data[i].x(); s[1] = h_data[i].y();
        s[2] = h_data[i].z(); s[3] = h_data[i].w();
        printf("[%d]: %06f, %06f, %06f, %06f\n",i, s[0], s[1], s[2], s[3]);
    }
}
// explicit initiate class templates
template class ode45<float,4>;
template class ode45<double,4>;
