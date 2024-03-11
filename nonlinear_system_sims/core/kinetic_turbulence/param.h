#ifndef PARAMS_H
#define PARAMS_H
#include "defines.h"
#include "linux_utility.h"
#include <ctime>

class Param
{
private:
    // IO
    std::string src_dir = SRC_PATH;
    std::string dir = std::string("_unnamed_result"); // maximum length allowed 255-85 = 170
    uint iter = 0;
    uint dumpInterval = 100;
    uint dumpIntervalShift = 0;

public:
    // simulation type
    algo::simType sim_type = algo::KINETIC_TURB3D;

    // Magic numbers
    const double pi = 3.14159265358979;
    const uint nk = 64; // explicitly defined in kernels
    const double lambda = 0.5;
    const double energy_cascade = -5.0/3.0;
    uint randomSeed = 777777;

    double L = 10.0;
    double eta = L/91.0; // 991: unstable when t > 5
    double E0 = 3.2;

    // Inertial particle
    double tau_inv =100.0;
    double g = -9.81; // assuming +Y direction

    // config
    glm::uvec3 dim = glm::uvec3(2048,16,1);
    double initDist = eta/2.0/pi;
    glm::uvec3 grid = glm::uvec3(16,16,1);

    // runtime
    float t0 = 0.0f;
    float t_end=30.0f;
    float dt = 0.002f;
    float t = t0;

    // numeric
    float tol = 1e-3f;

    Param()
    {
        checkDir(getDir().c_str()); // if not exist, create default dir
    }
    void genericDir()
    {
        time_t now = time(0);

        std::string dname = ctime(&now);
        dname.pop_back(); // remove spurious word

        // Check folder
        std::string path = SRC_PATH + dname;
        checkDir(path.c_str()); // if not exist, create

        // Remove space
        dir = dname;
    }
    std::string getDir() const
    {
        return (this->src_dir + this->dir);
    }

    void printInfo() const
    {
        std::cout << std::endl;
        printf("**Kinetic Simulation Configuration**\n");
        printf("L = %f, L/eta = %f, nk = %d, E0 = %f, seed = %d\n",L,L/eta, nk, E0, randomSeed);
        printf("nx = %d, ny = %d, nz = %d, initial distance = %f\n",dim.x,dim.y,dim.z,initDist);
        printf("tau_inv = %f, g = %f\n",tau_inv, g);
        printf("t0 = %f, t_end = %f, dt = %f, tol = %f, iter = %d\n",t0,t_end,dt,tol,iter);
        printf("output dir: %s\n", dir.c_str());
        printf("**Kinetic Simulation Configuration End**\n");
        std::cout << std::endl;
    }
    void printInfoToDir() const
    {
        FILE* fid;
        fid = fopen ((this->getDir() + "/" + "_config").c_str(), "w");

        //const char* _fid = string(path).append(name).c_str();
        //std::cout << _fid << std::endl;
        //FILE* fid = fopen (_fid,"w");
        if (fid)
        {
            // 1. Current time
            // 2. Compute mean square distance

            fprintf(fid,"OUTPUT_DIR\t%s\n", this->getDir().c_str()); // maximum length allowed 255-85 = 170
            fprintf(fid,"ITER\t%d\n", iter);
            fprintf(fid,"DUMP_INTERVAL\t%d\n", dumpInterval);
            fprintf(fid,"DUMP_INTERVAL_SHIFT\t%d\n", dumpIntervalShift);
            fprintf(fid,"\n");

            // simulation type
            fprintf(fid,"SIMULATION_TYPE\t%d\n", static_cast<int>(sim_type));
            fprintf(fid,"\n");

            // Magic numbers
            fprintf(fid,"PI\t%f\n", pi);
            fprintf(fid,"NK\t%d\n", nk);
            fprintf(fid,"LAMBDA\t%f\n", lambda);
            fprintf(fid,"CASCADE_LAW\t%f\n", energy_cascade);
            fprintf(fid,"RANDOM_SEED\t%d\n", randomSeed);
            fprintf(fid,"\n");

            fprintf(fid,"L\t%f\n", L);
            fprintf(fid,"ETA\t%f\n", eta);
            fprintf(fid,"L/ETA\t%f\n", L/eta);
            fprintf(fid,"E0\t%f\n", E0);
            fprintf(fid,"\n");

            // Inertial particle
            fprintf(fid,"TAU_INV\t%f\n", tau_inv);
            fprintf(fid,"G_Y\t%f\n", g);
            fprintf(fid,"\n");

            // config
            fprintf(fid,"DIM\t%d\t%d\t%d\n", dim.x,dim.y,dim.z);
            fprintf(fid,"INIT_DIST\t%f\n", initDist);
            //fprintf(fid,"GRID\t%f\n", grid.x,grid.y,grid.z);
            fprintf(fid,"\n");

            // runtime
            fprintf(fid,"T_INIT\t%f\n", t0);
            fprintf(fid,"T_END\t%f\n", t_end);
            fprintf(fid,"DT\t%f\n", dt);
            fprintf(fid,"TIME\t%f\n", t);
            fprintf(fid,"\n");

            // numeric
            fprintf(fid,"TOL\t%e\n", tol);

            fclose (fid);
        }else {
            std::cout<<"PARAM_OUTPUT::ERROR"<<std::endl;
        }


    }
    bool needDump() const
    {
        return (iter%dumpInterval == dumpIntervalShift);
    }
    bool timeOver() const
    {
        return !(t < t_end);
    }
    float getCurrentTimestep() const
    {
        return this->t;
    }
    void advCurrentTimestep()
    {
        this->t = getNextTimestep();
    }
    float getNextTimestep() const
    {
        float next_t = this->t + getCurrentTimeInterval();
        if(next_t > this->t_end) { next_t = this->t_end; }
        return (next_t);
    }
    float getCurrentTimeInterval() const
    {
        return (this->dt);
    }
    uint getIter() const
    {
        return (this->iter);
    }
    void operator ++ (int) noexcept
    {
        iter++;
    }

};
#endif
