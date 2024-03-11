#include "syclode45.hpp"
#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <memory>
#include "filesystem.hpp"

constexpr int NEQN = 8;
constexpr int ITER = 128;
static int REALIZATION = 1;

// A simple wrapper class
// todo: need a command line argument parser
class App
{
public:
    App() {}
    void initialize();
    int exec();
    //void finalize();

    void set_config_file(const char* in) { if(in) { config_file = in; std::cout << "Set config file from cmd:" << config_file << std::endl;} }
    void set_batch_option(bool flag) { batch_run = flag; }
    void set_precision(bool flag) { fp64 = flag; }
    void set_tendlog(float arg) { l0 = log10(arg); }
    void set_tbeginlog(float arg) { l1 = log10(arg); }

protected:
    std::vector<float> gen_timeserieslog() const
    {
        std::vector<float> t;

        for(int i = 0; i < ITER; i++)
        {
            t.push_back(pow( 10.0,l0 + i*(l1-l0)/(ITER-1)) );
        }

        return t;
    }
private:
    std::shared_ptr<SYCLbase> m_app;
    std::string config_file;
    bool batch_run = true;
    bool fp64 = true;

    float l0 = 2;
    float l1 = -1;
};



int main(int argc, char **argv)
{

    App app;

    if(argc > 1) { app.set_config_file(argv[1]); }
    if(argc > 2) { app.set_batch_option(atoi(argv[2]));}
    if(argc > 3) { app.set_precision(atoi(argv[3]));}
    if(argc > 4) { app.set_tendlog(atof(argv[4]));}
    if(argc > 5) { REALIZATION = (atoi(argv[5]));}

    app.initialize();

    return app.exec();
}


void App::initialize()
{
    // global configuration
    util::IO::set_output_option_append(false);

    if(fp64)
    {
        m_app.reset(new ode45<double, NEQN>);
        std::cout << "Use 64-bit precision." << std::endl;
    }
    else
    {
        m_app.reset(new ode45<float , NEQN>);
        std::cout << "Use 32-bit precision." << std::endl;
    }
}

int App::exec()
{
    do
    {
        auto start = std::chrono::steady_clock::now();
        if(m_app->load_param(config_file.c_str()) == false) { return 0; };

        for(int it = 0; it < REALIZATION; it++)
        {

            if(m_app->initialize() == false) { return 0; };

            // get timeseries (distribute evenly on logspace)
            auto ts = this->gen_timeserieslog();

            for(int i = 0; i < ITER; i++)
            {
                // getcurrent timepoint
                auto tf = ts.back();
                ts.pop_back();
                char fn[255];
                std::sprintf(fn,"%06d.dat",i);

                // compute
                m_app->calc(tf);
                m_app->finalize(fn);
            };

            m_app->incr_shuffle();

            util::IO::set_output_option_append(true);
        }


        auto end = std::chrono::steady_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Time: " << time << " fps: " << (double)ITER / (time / 1000.0f) << std::endl;

    }while (batch_run);

    return 0;
}

