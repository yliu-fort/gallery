#include "syclode45.hpp"
#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <string>
#include <memory>


constexpr int NEQN = 4;
constexpr int ITER = 16;
class App
{
public:
    App() {}
    void initialize();
    int exec();
    void finalize();

    void set_config_file(const char* in) { if(in) { config_file = in; std::cout << "Set config file from cmd:" << config_file << std::endl;} }
    void set_batch_option(bool flag) { batch_run = flag; }
    void set_precision(bool flag) { fp64 = flag; }

private:
    std::shared_ptr<SYCLbase> m_app;
    std::string config_file;
    bool batch_run = false;
    bool fp64 = true;
};



int main(int argc, char **argv)
{
    App app;

    if(argc > 1) { app.set_config_file(argv[1]); }
    if(argc > 2) {app.set_batch_option(atoi(argv[2]));}
    if(argc > 3) {app.set_precision(atoi(argv[3]));}

    app.initialize();

    return app.exec();
}


void App::initialize()
{
    if(fp64)
    {
        m_app.reset(new ode45<double, NEQN>);
    }
    else
    {
        m_app.reset(new ode45<float , NEQN>);
    }
}

int App::exec()
{
    do
    {
        auto start = std::chrono::steady_clock::now();
        if(m_app->initialize(config_file.c_str()) == false) { return 0; };

        for(int i = 0; i < ITER; i++)
        {
            // SYCL module
            m_app->calc();

        };

        m_app->finalize();

        auto end = std::chrono::steady_clock::now();
        auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
        std::cout << "Time: " << time << " fps: " << (double)ITER / (time / 1000.0f) << std::endl;

    }while (batch_run);

    return 0;
}

