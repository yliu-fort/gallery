#include <iostream>
#include <string>

#define CL_GL_INTEROP
#include "syclode45.hpp"
#include "sycl_integrator_rk45.hpp"
#include "sycl_integrator_adm4.hpp"
#include "ks3d_algorithm.hpp"

using namespace cl::sycl;

// state_type = double
constexpr access::mode sycl_read = access::mode::read;
constexpr access::mode sycl_write = access::mode::write;
constexpr access::mode sycl_rw = access::mode::read_write;

// flags
#define HAVE_EVENT_FUN false

/* This is the class used to name the kernel for the runtime.
 * This must be done when the kernel is expressed as a lambda. */
template <typename T, int N>
class _lambda_ode45;


template <typename T, int N>
void ode45<T, N>::setup()
{
    // Override default setting
    imported_param.L /= 10;
    imported_param.eta /= 1000;
    imported_param.tau_inv = std::numeric_limits<T>::infinity();
    imported_param.g = 0;
    imported_param.xoffset = 13;
    imported_param.lambda = 0.5;
    imported_param.tol = 1e-3;
    imported_param.absTol = 1e-2;
    imported_param.set_t_factor(1.0);

    buf_s = buffer<T, 1>(range<1>(array_size));

    std::vector<T> B;
    algo::kinetic_simulation_configuration3D(B, imported_param);
    param = buffer<T, 1>(B.begin(), B.end()).template reinterpret<cl::sycl::vec<T,4>>(range<1>(B.size()/4));

}

template <typename T, int N>
bool ode45<T, N>::initialize(const char* num_item)
{
    // read parameters
    array_size = std::stoul(std::string(num_item));

    setup();

    std::vector<vec> A(array_size);
    auto gen = [=](){ return 2*rand()/(double)RAND_MAX-1; };
    std::generate(A.begin(), A.end(), [&](){ return vec(gen(),gen(),gen(),0,0,0,0,0); });
    buf = buffer<vec, 1>(A.begin(), A.end());//.template reinterpret<vec>(range<1>(array_size));

    return true;
}


#ifdef CL_GL_INTEROP
template <typename T, int N>
bool ode45<T, N>::initialize(T* gl_ptr, unsigned int gl_buf, const char* num_item)
{
    // read parameters
    array_size = std::stoul(std::string(num_item));

    setup();

    T RADIUS = 0.1*imported_param.eta;
    srand(3876953);
    auto gen = [=](){ return RADIUS*(2*rand()/(double)RAND_MAX-1); };
    //std::generate(gl_ptr, gl_ptr+array_size*N, [&](){ return gen(); });
    std::generate(reinterpret_cast<vec*>(gl_ptr), reinterpret_cast<vec*>(gl_ptr)+array_size,
                  [&](){ return vec(gen(),gen(),gen(),0,0,0,0,0); });


    SYCLbase::create_from_gl_buffer(mem_obj, CL_MEM_READ_WRITE, static_cast<cl_GLuint>(gl_buf)); // err handled by base class
    buf = cl::sycl::buffer<vec, 1>(mem_obj, getContext());

    clFinish(getQueue().get());

    return true;
}
#endif

template <typename T>
class printkernel;

template <typename T, int N>
void ode45<T, N>::calc(float tf)
{
#ifdef CL_GL_INTEROP
    clEnqueueAcquireGLObjects(getQueue().get(), 1, &mem_obj, 0,0,0);
#endif

    bool isPassive = imported_param.is_passive();

    submit([&](handler& cgh) {
        auto val_t = buf.template get_access<sycl_rw>(cgh);
        auto time_t = buf_s.template get_access<sycl_rw>(cgh);
        auto param_t = param.template get_access<sycl_read, access::target::constant_buffer>(cgh);

        //cl::sycl::stream out(1024, 256, cgh);

        algo::KineticTurb3DPassive<T,N> func_passive(val_t, param_t);
        algo::KineticTurb3D<T,N> func_inertial(val_t, param_t, HAVE_EVENT_FUN);

        auto kern = [=](id<1> idx) {

            auto x = val_t[idx];
            auto t = time_t[idx];

            if(isPassive)
                algo::integrator_rk45(func_passive, t, time_t[0]+(T)0.002, x, (T)0.002);
            else
                algo::integrator_rk45(func_inertial, t, time_t[0]+(T)0.002, x, (T)0.002);

            val_t[idx] = x;
            time_t[idx] = t;
            //out << x << ", " << t << cl::sycl::endl;

        };
        cgh.parallel_for<class _lambda_ode45<T, N>>(range<1>(array_size), kern);
    });

}

template <typename T, int N>
void ode45<T, N>::finalize(const char* filename)
{

    // query host pointer will hit the performance heavily. only do it when debugging
    wait_and_throw();

#ifdef CL_GL_INTEROP
    clEnqueueReleaseGLObjects(getQueue().get(), 1, &mem_obj, 0,0,0);
#endif
}

template <typename T, int N>
void ode45<T, N>::finalize(T* ptr, int stride)
{
    // query host pointer will hit the performance heavily. only do it when debugging

    auto h_data = buf.template get_access<sycl_read>();
    auto t_data = buf_s.template get_access<sycl_read>();

    int size = stride==0?N:stride;
    for(int i = 0; i < array_size; i++)
    {
        ptr[size*i+0] = h_data[i].s0();
        ptr[size*i+1] = h_data[i].s1();
        ptr[size*i+2] = h_data[i].s2();

        ptr[size*i+4] = h_data[i].s4();
        ptr[size*i+5] = h_data[i].s5();
        ptr[size*i+6] = h_data[i].s6();

    }
}

// explicit initiate class templates
template class ode45<float,8>;
template class ode45<double,8>;
