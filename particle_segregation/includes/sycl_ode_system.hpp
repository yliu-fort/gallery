#pragma once

#include <CL/sycl.hpp>

namespace algo
{


template<typename T, int N>
class ODEsystem
{
    using acc_t = cl::sycl::accessor
    <cl::sycl::vec<T,N>,1,
    cl::sycl::access::mode::read_write,
    cl::sycl::access::target::global_buffer>;

public:
    ODEsystem(acc_t p) : ptr(p) {}

    void operator()( const cl::sycl::vec<T,N> &x , cl::sycl::vec<T,N> &dxdt , const T t ) const
    {
        dxdt = (T)3.0/((T)2.0*t*t) + x/(2.0*t);
    }

    bool event ( const cl::sycl::vec<T,N> &x , const cl::sycl::vec<T,N> &dxdt , const T t ) const
    {
        bool trigger = (T)x.y() > (T)1.0 && (T)dxdt.x() > (T)0; // passing through x = 1.0 from bottom only

        // when a event is triggered, program will refine the timestep until convergence
        return trigger;
    }

    void output ( const cl::sycl::vec<T,N> &x , const cl::sycl::vec<T,N> &dxdt , const T t ) const
    {
        //out << "debug ostream" << sycl::endl;
    }

    // ODEoptions
    void setODEoption(T& rtol, T& atol, T& max_dt, bool& have_event_fcn, bool& have_output_fcn) const
    {
        rtol = relTol;
        atol = absTol;
        max_dt = maxDeltaT;
        have_event_fcn = haveEventFcn;
    }


private:
    T relTol = 1e-6;
    T absTol = 1e-6;
    T maxDeltaT = 0.01;
    bool haveEventFcn = true;

    acc_t ptr;
};

}


