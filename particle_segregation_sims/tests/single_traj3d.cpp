#include <iostream>
#include "syclode45.hpp"
#include "param.h"
#include "ks3d_algorithm.hpp"
#include "sycl_integrator_rk45.hpp"
#include "sycl_integrator_adm4.hpp"
#include <cmath>
#include <string>
#include <memory>

using real = double;
constexpr int NEQN = 8;
#define IDX (0)


int main()
{
    ode45<real, NEQN> m_app;

    if(m_app.initialize("input.txt"))
    {
        m_app.calc(100);
        m_app.finalize();
    }

    return 0;
}


using namespace cl::sycl;


template<typename T, int N>
class InlineKS3D : public algo::KineticTurb3D<T,N,access::target::host_buffer,access::target::host_buffer>
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

    using base_class = algo::KineticTurb3D<T,N,access::target::host_buffer,access::target::host_buffer>;

public:
    InlineKS3D(acc_t p, cst_t c): base_class(p, c) {}

    void output ( const Vec &x , const Vec &dxdt , const T t ) const
    {
        printf("%e\t%e\t%e\t%e\n",double(x.s0()),double(x.s1()),double(x.s2()),t);
    }

    void setODEoption ( T& rtol, T& atol, T& max_dt, bool& have_event_fcn, bool& have_output_fcn ) const
    {
        base_class::setODEoption(rtol, atol, max_dt, have_event_fcn, have_output_fcn);

        have_event_fcn = false;
        have_output_fcn = true;
    }

};

template<typename T, int N>
class InlineKS3DPassive : public algo::KineticTurb3DPassive<T,N,access::target::host_buffer,access::target::host_buffer>
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

    using base_class = algo::KineticTurb3DPassive<T,N,access::target::host_buffer,access::target::host_buffer>;

public:
    InlineKS3DPassive(acc_t p, cst_t c): base_class(p, c) {}

    void output ( const Vec &x , const Vec &dxdt , const T t ) const
    {
        printf("%e\t%e\t%e\t%e\n",double(x.s0()),double(x.s1()),double(x.s2()),t);
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
    tf *= imported_param.get_kolmogorov_time();
    bool isPassive = imported_param.is_passive();

    {
        auto val_t = buf.template get_access<sycl_rw>();
        auto time_t = buf_s.template get_access<sycl_rw>();
        auto param_t = param.template get_access<sycl_read>();


        InlineKS3D<T,N> func_inertial(val_t, param_t);
        InlineKS3DPassive<T,N> func_passive(val_t, param_t);

        auto x = val_t[IDX];
        auto t = time_t[IDX];
        t = 0; // no RAII when run on host side

        if(isPassive)
            algo::integrator_rk45(func_passive, t, (T)tf, x, init_dt);
        else
            algo::integrator_rk45(func_inertial, t, (T)tf, x, init_dt);

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
    imported_param.print_param();
    array_size = imported_param.dim[0];

    // set/get params
    T L = (T)imported_param.L;
    T fdist = (T)imported_param.b*(T)imported_param.L; // falling distance
    T zdist = (T)imported_param.b*(T)imported_param.eta; // falling distance
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

        T gtau = buf[7*param.nk+2];
        auto stride = buf[7*param.nk+3];

        // shifted pos
        T px = x + stride; // note: need a random machine...
        T py = y + stride;
        T pz = z + stride;

        T vel[3]{0};

        // compute ks parameters
        for(unsigned int i = 0; i < param.nk; i++)
        {
            T KS,dKS;

            KS =      px*buf[i+3*param.nk]
                    + py*buf[i+4*param.nk]
                    + pz*buf[i+5*param.nk]
                    + buf[i+6*param.nk]*(t + stride);
            dKS = cos( KS ) - sin( KS );

            vel[0] += buf[i+0*param.nk]*dKS;
            vel[1] += buf[i+1*param.nk]*dKS;
            vel[2] += buf[i+2*param.nk]*dKS;

        }

        // add external force
        vel[1] += gtau;

        return vel[dim];
    };

    // init state vector [pos vel]
    int i = 0;
    std::vector<vec> A(array_size);

    // fill buffer
    std::generate(A.begin(), A.end(), [&](){
        T px = (++i)*imported_param.initDist;
        T py = 0.0;
        T pz = zdist;
        return vec(px,py,pz,0,uf(px,py,pz,0,0),uf(px,py,pz,0,1),uf(px,py,pz,0,2),0);
    });

    // transfer buffer from host to device
    buf = buffer<vec, 1>(A.begin(), A.end());

    return true;
}

// explicit initiate class templates
template class ode45<float,8>;
template class ode45<double,8>;
