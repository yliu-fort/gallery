#pragma once
// Adams-Bashforth-Moulton 4th-order multistep prediction correction method
// notice that event system is not fully implemented
#include <CL/sycl.hpp>

#include "sycl_integrator_rk45.hpp"

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
void integrator_adm4(Func&& odeFun,
                     T& t0,
                     const T tf,
                     cl::sycl::vec<T,N>& y0,
                     const T init_t = 0.001
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
    T b13 = static_cast<T>(  55)/static_cast<T>(24);
    T b12 = static_cast<T>( -59)/static_cast<T>(24);
    T b11 = static_cast<T>(  37)/static_cast<T>(24);
    T b10 = static_cast<T>(  -9)/static_cast<T>(24);

    T b24 = static_cast<T>( 251)/static_cast<T>(720);
    T b23 = static_cast<T>( 646)/static_cast<T>(720);
    T b22 = static_cast<T>(-264)/static_cast<T>(720);
    T b21 = static_cast<T>( 106)/static_cast<T>(720);
    T b20 = static_cast<T>( -19)/static_cast<T>(720);


    // initial timestep
    T hmin = 1e-6;
    T absh = init_t;
    T h = tdir * absh;

    decltype(y) y1=y;
    T t1 = t;

    decltype(y) y2=y;
    integrator_rk45(odeFun, t, t + h, y, absh);
    T t2 = t;

    decltype(y) y3=y;
    integrator_rk45(odeFun, t, t + h, y, absh);
    T t3 = t;

    //decltype(y) y4=y;
    integrator_rk45(odeFun, t, t + h, y, absh);
    //T t4 = t;



    decltype(y) f2=f1, f3=f1, f4=f1, f5=f1, f6=f1, f7=f1;
    //T t2 = t, t3 = t, t4 = t, t5 = t, t6 = t, t7 = t;
    T err=0;
    T tnew = t, tout_new = t;
    decltype(y) ynew=y0, yout_new=y0;

    bool done = false;

    while (~done)
    {
        // No expand to the exact tfinal
        // This will cause truncation error, that is, tout != tfinal, i.e.
        // 9.9999999999987 != 10
        //if (1.1*absh >= ABS(tfinal-t)) //cl::sycl::abs???
        //{
        //    h = tfinal-t;
        //    absh = ABS(h); //cl::sycl::abs???
        //    done = true;
        //}

        // LOOP FOR ADVANCING ONE STEP.
        bool nofailed = true;// no failed attempts

        // Predictor
        odeFun(y1,f1,t1);
        odeFun(y2,f2,t2);
        odeFun(y3,f3,t3);
        odeFun(y,f4,t);

        decltype(y) yp = y + h*(b13*f4 + b12*f3 + b11*f2 + b10*f1);

        tnew = t + h;

        //while (true)
        for(int corr = 0; corr < 3; corr++)
        {

            //if (done)
            //    tnew = tfinal;   // Hit end point exactly.
            //h = tnew - t;      // Purify h.

            // corrector
            odeFun(yp,f5,tnew);
            //ynew = y + h* ( b61*f1 + b63*f3 + b64*f4 + b65*f5 + b66*f6 );
            ynew = y + h*(b24*f5 + b23*f4 + b22*f3 + b21*f2 + b20*f1);
            odeFun(ynew, f7, tnew);

            //nfevals = nfevals + 6;

            // Estimate the error.
            bool NNrejectStep = false;
            //auto fE = f1*e1 + f3*e3 + f4*e4 + f5*e5 + f6*e6 + f7*e7;
            auto fE = ynew - yp;

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
                done = false;
                yp = ynew;
            }
            else                                // Successful step
            {
                //y = ynew;
                //NNreset_f7 = false;
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


                integrator_rk45(odeFun, tnew, tfinal, ynew);

                done = true;
            }

#endif
            // decide how to refine timestep..
        }

        //if (output_ty )
        {
#ifndef SIMPLE_ODE_ONLY
            if(haveOutoutFcn)
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
                if ((tfinal - hmin) < tnew)
                {
                    tout_new = tnew;
                    yout_new = ynew;
                    done = true;
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

        // Advance the integration one step.
        // todo: use warp indexing to reduce actual copy & register usage
        t1 = t2;
        t2 = t3;
        t3 = t;

        y1 = y2;
        y2 = y3;
        y3 = y;

        t = tnew;
        y = ynew;
    }

    t0 = tout_new;
    y0 = yout_new;
}

#undef ABS
}
