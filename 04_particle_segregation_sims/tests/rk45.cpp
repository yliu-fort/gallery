#include <iostream>
#include <algorithm>
#include <math.h>

//#define SIMPLE_ODE_ONLY
#include "sycl_integrator_rk45.hpp"
#include "sycl_integrator_adm4.hpp"
#include "sycl_integrator_rkf45.hpp"
#include "sycl_ode_system.hpp"


using num_t = double;
constexpr int num_ode = 4;
using vec = cl::sycl::vec<num_t,num_ode>;

void rhs( const vec &x , vec &dxdt , const num_t t )
{
    dxdt = (num_t)3.0/((num_t)2.0*t*t) + x/((num_t)2.0*t);
}

class odefun
{
public:
    void operator()( const vec &x , vec &dxdt , const num_t t )
    {
        //dxdt = (num_t)3.0/((num_t)2.0*t*t) + x/((num_t)2.0*t);
        dxdt = vec(x.y(),(1-x.x()*x.x())*x.y()-x.x(),0,0);
    }
    bool event ( const vec &x , const vec &dxdt , const num_t t ) const
    {
        return false;
    }

    void output ( const vec &x , const vec &dxdt , const num_t t )
    {
        //std::cout << t << "\t" << x.x() << ", " << x.y() << std::endl;
        if(t > count)
        {
            count++;
            printf("%f\t%f\t%f\n",t,double(x.x()),double(x.y()));
        }
    }

    // ODEoptions
    void setODEoption(num_t& rtol, num_t& atol, num_t& max_dt, bool& have_event_fcn, bool& have_output_fcn) const
    {
        rtol = 1e-10;
        atol = 1e-10;
        have_output_fcn = true;
        return;
    }
private:
    int count = 0;

};

int main()
{
    std::cout.precision(17);
    vec y(2,0,0,0);
    num_t t0 = 0.0;

    algo::integrator_rk45(odefun(), t0, (num_t)20.0, y);
    std::cout << "RK45::t = " << t0 << " y final = " << y.x() << ", " << y.y() << "\n";

    y = vec(2,0,0,0);
    t0 = 0.0;

    algo::integrator_adm4(odefun(), t0, (num_t)20.0, y);
    std::cout << "ADM4::t = " << t0 << " y final = " << y.x() << ", " << y.y() << "\n";

    y = vec(2,0,0,0);
    t0 = 0.0;

    algo::integrator_rkf45(odefun(), t0, (num_t)20.0, y);
    std::cout << "RKF45::t = " << t0 << " y final = " << y.x() << ", " << y.y() << "\n";

}
