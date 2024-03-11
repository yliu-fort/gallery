#pragma once
// Adaptive Runge-Kutta 4/5th Fehlberg method
// notice that event system is not fully implemented
#include <CL/sycl.hpp>

namespace algo
{

#define ABS(x) ( cl::sycl::sign(x)*(x) ) // workaround for unsupported operations on NV card
// T must be scalar
// T must be scalar or vector with support of math opertion
// following operation has been used
// element-wise +, -, *, /, =
// sign, abs, min, max, pow, norm(reduction in vec)
// scalar-scalar, scalar-vector, vector-vector operations
template<typename Func, typename T, int N>
void integrator_rkf45(Func&& odeFun,
                     T& t0,
                     const T tf,
                     cl::sycl::vec<T,N>& y0,
                     const T init_t = 0
        )
{
    T t = t0;
    T tfinal = tf;
    T tspan = tf-t0;
    auto y = y0;
    decltype(y) f1;
    odeFun(y,f1,t);

    // ODE config
    T rtol = 1e-6;
    T atol = 1e-6;
    T tmax = 0.1;
    bool haveEventFcn = false;
    bool haveOutoutFcn = false;
    bool output_sol = true;


#ifndef SIMPLE_ODE_ONLY
    odeFun.setODEoption(rtol, atol, tmax, haveEventFcn, haveOutoutFcn);
#endif


    bool NNreset_f7 = false;
    decltype(y) threshold(atol / rtol);
    T tdir = cl::sycl::sign(tfinal - t0); // sign

    // Initialize method parameters.
    T exponent = static_cast<T>(1)/static_cast<T>(4);
    //std::array<T> A = {1/4, 3/8, 12/13, 1, 1/2}; // Still used by restarting criteria
    // B = [
    //     1/4         3/32    1932/2197     439/216          -8/27        25/216
    //     0           9/32   -7200/2197       -8               2            0
    //     0           0       7296/2197    3680/513       -3544/2565    1408/2565
    //     0           0       0            -845/4104       1859/4104    2197/4104
    //     0           0       0                0            -11/40        -1/5
    //     0           0       0                0               0            0
    //     0           0       0                0               0            0
    //     ];
    // E = [ 1/360, 0, -128/4275, -2197/75240, 1/50, 2/55];

    // Same values as above extracted as scalars (1 and 0 values ommitted)
    T a2=static_cast<T>(1)/static_cast<T>(4);
    T a3=static_cast<T>(3)/static_cast<T>(8);
    T a4=static_cast<T>(12)/static_cast<T>(13);
    //T a5=static_cast<T>(1)/static_cast<T>(1);
    T a6=static_cast<T>(1)/static_cast<T>(2);

    T b11=static_cast<T>(1     )/static_cast<T>(4    );
    T b21=static_cast<T>(3     )/static_cast<T>(32   );
    T b31=static_cast<T>(1932  )/static_cast<T>(2197 );
    T b41=static_cast<T>(439   )/static_cast<T>(216  );
    T b51=static_cast<T>(-8    )/static_cast<T>(27   );
    T b61=static_cast<T>(25    )/static_cast<T>(216  );
    T b22=static_cast<T>(9     )/static_cast<T>(32   );
    T b32=static_cast<T>(-7200 )/static_cast<T>(2197 );
    T b42=static_cast<T>(-8    )/static_cast<T>(1    );
    T b52=static_cast<T>(2     )/static_cast<T>(1    );
    T b33=static_cast<T>(7296  )/static_cast<T>(2197 );
    T b43=static_cast<T>(3680  )/static_cast<T>(513  );
    T b53=static_cast<T>(-3544 )/static_cast<T>(2565 );
    T b63=static_cast<T>(1408  )/static_cast<T>(2565 );
    T b44=static_cast<T>(-845  )/static_cast<T>(4104 );
    T b54=static_cast<T>(1859  )/static_cast<T>(4104 );
    T b64=static_cast<T>(2197  )/static_cast<T>(4104 );
    T b55=static_cast<T>(-11   )/static_cast<T>(40   );
    T b65=static_cast<T>(-1    )/static_cast<T>(5    );

    T e1=static_cast<T>(1     )/static_cast<T>(360   );
    T e3=static_cast<T>(-128  )/static_cast<T>(4275  );
    T e4=static_cast<T>(-2197 )/static_cast<T>(75240 );
    T e5=static_cast<T>(1     )/static_cast<T>(50    );
    T e6=static_cast<T>(2     )/static_cast<T>(55    );

    //T hmin = 16*std::nextafter(0,1.f); see matlab ode45 definition
    T hmin = 1e-7;
    T hmax = cl::sycl::min((T)1e15, tmax);
    T htry = cl::sycl::max((T)hmin, init_t);


    T absh = cl::sycl::min(hmax, cl::sycl::max(hmin, htry));

    decltype(y) y2, y3, y4, y5, y6, y7;
    decltype(y) f2=f1, f3=f1, f4=f1, f5=f1, f6=f1, f7=f1;
    T t2 = t, t3 = t, t4 = t, t5 = t, t6 = t, t7 = t;
    T err=0;
    T tnew = t, tout_new = t;
    decltype(y) ynew=y0, yout_new=y0;

    bool done = false;
    while (~done)
    {
        // By default, hmin is a small number such that t+hmin is only slightly
        // different than t.  It might be 0 if t is 0.
        //hmin = 16*eps(t);
        absh = cl::sycl::min(hmax, cl::sycl::max(hmin, absh));    // couldn't limit absh until new hmin
        T h = tdir * absh;

        // Stretch the step if within 10% of tfinal-t.
        if (1.1*absh >= ABS(tfinal-t)) //cl::sycl::abs???
        {
            h = tfinal-t;
            absh = ABS(h); //cl::sycl::abs???
            done = true;
        }

        // LOOP FOR ADVANCING ONE STEP.
        bool nofailed = true;// no failed attempts
        while (true)
        {

            y2 = y + h * (b11*f1 );
            t2 = t + h * a2;
            odeFun(y2,f2,t2);

            y3 = y + h * (b21*f1 + b22*f2 );
            t3 = t + h * a3;
            odeFun(y3,f3,t3);

            y4 = y + h * (b31*f1 + b32*f2 + b33*f3 );
            t4 = t + h * a4;
            odeFun(y4,f4,t4);

            y5 = y + h * (b41*f1 + b42*f2 + b43*f3 + b44*f4 );
            t5 = t + h;
            odeFun(y5,f5,t5);

            y6 = y + h * (b51*f1 + b52*f2 + b53*f3 + b54*f4 + b55*f5 );
            t6 = t + h * a6;
            odeFun(y6,f6,t6);

            tnew = t + h;
            if (done)
                tnew = tfinal;   // Hit end point exactly.
            h = tnew - t;      // Purify h.

            ynew = y + h* ( b61*f1 + b63*f3 + b64*f4 + b65*f5 );
            odeFun(ynew,f7,tnew);

            //nfevals = nfevals + 6;

            // Estimate the error.
            bool NNrejectStep = false;
            auto fE = f1*e1 + f3*e3 + f4*e4 + f5*e5 + f6*e6;


#ifdef ABSOLUTE_ERR_CONTROL
            auto verr = absh * ABS(fE);
#else
            auto verr = absh * ABS(fE) / cl::sycl::max(cl::sycl::max(ABS(y),ABS(ynew)),threshold);
#endif
            T va = cl::sycl::length(cl::sycl::vec<T,(N+1)/2>(verr.lo()));//8-4
            T vb = cl::sycl::length(cl::sycl::vec<T,(N+1)/2>(verr.hi()));//8-4
            err = cl::sycl::max(va,vb);

            // Accept the solution only if the weighted error is no more than the
            // tolerance rtol.  Estimate an h that will yield an error of rtol on
            // the next step or the next try at taking this step, as the case may be,
            // and use 0.8 of this value to avoid failures.
            if (err > rtol)                       // Failed step
            {

                //nfailed = nfailed + 1;
                if (absh <= hmin)
                {
                    //printf('MATLAB:ode45:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin ));
                    return;
                }


                if (nofailed)
                {
                    nofailed = false;
                    if (NNrejectStep)
                        absh = cl::sycl::max(hmin, (T)0.5*absh);
                    else
                        absh = cl::sycl::max(hmin, absh * cl::sycl::max((T)0.1, cl::sycl::pow((T)0.84*(rtol/err),exponent)));
                }
                else
                {
                    absh = cl::sycl::max(hmin, (T)0.5 * absh);
                }


                h = tdir * absh;
                done = false;
            }
            else                                // Successful step
            {
                NNreset_f7 = false;
                break;
            }
        } // main loop end
        //nsteps = nsteps + 1;

        bool stop = false;
        if (haveEventFcn)
        {
#ifndef SIMPLE_ODE_ONLY
            stop = odeFun.event(ynew,f7,tnew);

            if(stop)
            {

                tnew = t;
                ynew = y;

                if(absh > hmin)
                {
                    absh /= 2.0;
                    NNreset_f7 = true; // reset f1
                    nofailed = false;
                    stop = false;
                }else {
                    done = true;
                }
            }

#endif
            // decide how to refine timestep..
        }

        //if (output_ty )
        {
#ifndef SIMPLE_ODE_ONLY
            if(haveOutoutFcn && nofailed)
                odeFun.output(ynew,f7,tnew);
#endif

            //while (next <= ntspan)
            {
                if( tdir * (tnew - tfinal) < 0)
                {
                    if (haveEventFcn && stop )    // output tstop,ystop
                    {

                        //nout_new = nout_new + 1;
                        tout_new = tnew;
                        yout_new = ynew;
                    }
                }

                //nout_new = nout_new + 1;
                //tout_new = [tout_new, tspan(next)];
                if (tfinal == tnew)
                {
                    tout_new = tfinal;
                    yout_new = ynew;
                }
                //next = next + 1;
            }

            if(stop)
            {
                done = true;
            }

        }

        if (done)
            break;

        // If there were no failures compute a new h.
        if (nofailed)
        {
            // Note that absh may shrink by 0.8, and that err may be 0.
            T temp = 1.25*cl::sycl::pow(err/rtol, exponent);
            if (temp > 0.2)
            {
                absh = absh / temp;
            }
            else
            {
                absh = 5.0*absh;
            }
            //std::cout << pow(err/rtol, exponent) << " ::absh adjust=" << absh << std::endl;
        }

        // Advance the integration one step.
        t = tnew;
        y = ynew;

        //std::cout << "advance---\n";
        if (NNreset_f7)
        {
            // Used f7 for unperturbed solution to interpolate.
            // Now reset f7 to move along constraint.
            odeFun(ynew,f7,tnew);
            //nfevals = nfevals + 1;

        }
        f1 = f7;  // Already have f(tnew,ynew)
    }

    t0 = tout_new;
    y0 = yout_new;
}

#undef ABS
}
