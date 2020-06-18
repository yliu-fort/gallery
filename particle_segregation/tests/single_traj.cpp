#include <iostream>
#include "syclode45.hpp"
#include "param.h"
#include "ks2d_algorithm.hpp"
#include "sycl_integrator_rk45.hpp"
#include "sycl_integrator_adm4.hpp"
#include <cmath>
#include <string>
#include <memory>

using real = double;
constexpr int NEQN = 4;
#define IDX 3680


int main()
{
    ode45<real, NEQN> m_app;

    if(m_app.initialize("validate_input.txt"))
    {
        m_app.calc();
        m_app.finalize();
    }

    return 0;
}


using namespace cl::sycl;


template<typename T, int N>
class InlineKS2D : public algo::KineticTurb2D<T,N,access::target::host_buffer,access::target::host_buffer>
{
    using Vec = vec<T,N>;
    using vec4 = vec<T,4>;

    using acc_t = accessor
    <Vec,1,
    access::mode::read_write,
    access::target::host_buffer>;

    using cst_t = accessor
    <vec4,1,
    access::mode::read,
    access::target::host_buffer>;

    using base_class = algo::KineticTurb2D<T,N,access::target::host_buffer,access::target::host_buffer>;

public:
    InlineKS2D(acc_t p, cst_t c): base_class(p, c) {}

    void output ( const Vec &x , const Vec &dxdt , const T t ) const
    {
        printf("%e\t%e\t%e\t%e\t%e\n",double(x.x()),double(x.y()),double(x.z()),double(x.w()),t);
    }

    void setODEoption ( T& rtol, T& atol, T& max_dt, bool& have_event_fcn, bool& have_output_fcn ) const
    {
        base_class::setODEoption(rtol, atol, max_dt, have_event_fcn, have_output_fcn);

        have_event_fcn = false;
        have_output_fcn = true;
    }

};


// state_type = double
constexpr access::mode sycl_read = access::mode::read;
constexpr access::mode sycl_write = access::mode::write;
constexpr access::mode sycl_rw = access::mode::read_write;

/* This is the class used to name the kernel for the runtime.
 * This must be done when the kernel is expressed as a lambda. */
template <typename T, int N>
class _lambda_ode45;


template <typename T, int N>
void ode45<T, N>::calc(float tf)
{

    T init_dt = imported_param.get_max_allowed_dt();

    {
        auto val_t = buf.template get_access<sycl_rw>();
        auto time_t = buf_s.template get_access<sycl_rw>();
        auto param_t = param.template get_access<sycl_read>();


        InlineKS2D<T,N> odeHandle(val_t, param_t);


        auto x = val_t[IDX];
        auto t = time_t[IDX];
        t = 0; // no RAII when run on host side

        algo::integrator_rk45(odeHandle, t, (T)(100.0), x, init_dt);

        val_t[IDX] = x;
        time_t[IDX] = t;

    }

}

template <typename T, int N>
void ode45<T, N>::finalize(const char* filename)
{
    // query host pointer will hit the performance heavily. only do it when debugging
    wait_and_throw();
}

template <typename T, int N>
bool ode45<T, N>::initialize(const char* config_file)
{
    // read parameters
    if(!imported_param.config_param(config_file)) { return false; }
    array_size = imported_param.dim[0];

    T L = (T)imported_param.L;
    T fdist = T(imported_param.b*imported_param.L); // falling distance
    T distx = T(imported_param.eta*0.1);

    int i = 0;
    std::vector<vec> A(array_size);
    std::generate(A.begin(), A.end(), [&](){ return vec(-0.5*L + (++i)*distx,fdist,0,0); });
    buf = buffer<vec, 1>(A.begin(), A.end());//.template reinterpret<vec>(range<1>(array_size));

    buf_s = buffer<T, 1>(range<1>(array_size));


    std::vector<T> B;
    algo::kinetic_simulation_configuration2D(B, imported_param);
    param = buffer<T, 1>(B.begin(), B.end()).template reinterpret<cl::sycl::vec<T,4>>(range<1>(B.size()/4));

    return true;
}

// explicit initiate class templates
template class ode45<float,4>;
template class ode45<double,4>;
