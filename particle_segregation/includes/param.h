#ifndef PARAMS_H
#define PARAMS_H
#include "filesystem.hpp"
#include <iostream>
#include <string>
#include <map>

namespace algo
{
enum simType
{
    DOUBLE_GYRE=0,
    DOUBLE_GYRE3D=1,
    KINETIC_TURB2D=2,
    KINETIC_TURB3D=3,
    INERTIAL_PARTICLE=4,
    NUM_SIMTYPE=5
};
enum initType
{
    ZERO_VEL=0,
    ZERO_ACC=1,
    FLUID_VEL=2,
    NUM_INITTYPE = 3
};

class Param : public util::IO
{
private:
    // IO
    std::string outputFolder;
    bool reducedOutput = true;
    // simulation type
    algo::simType sim_type = algo::KINETIC_TURB2D;

public:

    // Magic numbers
    const double pi = 3.14159265358979;
    const double energy_cascade = 5.0/3.0;
    unsigned int nk = 64; // explicitly defined in kernels
    double lambda = 0.5;
    unsigned int randomSeed = 777777;
    double xoffset = 0.0;

    double L = 10.0;
    double eta = L/91.0; // 991: unstable when t > 5
    double E0 = 3.2;

    // Inertial particle
    double tau_inv = 10.0;
    double g = -9.81; // assuming +Y direction
    double b = L; // falling distance

    // config
    unsigned int dim[3]{1024,1,1};
    //glm::uvec3 dim = glm::uvec3(2048,16,1);
    double initDist = eta/2.0/pi;
    //glm::uvec3 grid = glm::uvec3(16,16,1);

    // 0: zero acc 1: zero vel
    initType initApproach = ZERO_ACC;

    // runtime
    double t0 = 0.0f;
    double t_end= 30.0f;

    // numeric
    double tol = 1e-3f;
    double absTol = 1e-2f;
    double t_factor = 0.01;

    bool isPassive = false;

    inline int get_chunk_size() const
    {
        return nk/4;
    }

    // todo: remove pre-defined energy cascade 5/3
    // magic number here because we use eqn in Fung 2003
    // to compute E(k), i.e. E(k) = E_0*L*(kL)^(-5/3)
    inline double get_turbulent_rms_velocity() const
    {
        return 0.5419*sqrt(E0);
    }

    // todo: remove pre-defined energy cascade 5/3
    inline double get_kolmogorov_time() const
    {
        return L/get_turbulent_rms_velocity()*pow(eta/L,2.0/3.0);
    }

    // maximum timestep is bounded to 0.01*kolmogorov_time
    inline double get_max_allowed_dt() const
    {
        return this->t_factor*this->get_kolmogorov_time();
    }

    inline void set_t_factor(const double& arg)
    {
        this->t_factor = arg;
    }

    inline bool is_passive() const
    {
        if(std::isinf(tau_inv))
            return true;
        else
            return false;
    }

    inline double get_drift_velocity() const
    {
        if(tau_inv == 0)
            return 0.0;
        return this->g/this->tau_inv;
    }

    // push required parameters into constant buffer
    template<typename T>
    void push_params(std::vector<T>& buffer) const
    {
        buffer.push_back(T(this->get_chunk_size()       ));
        buffer.push_back(T(this->tau_inv                ));
        buffer.push_back(T(this->get_drift_velocity()   ));
        buffer.push_back(T(this->xoffset                ));

        buffer.push_back(T(this->tol                    ));
        buffer.push_back(T(this->absTol                 ));
        buffer.push_back(T(this->get_max_allowed_dt()   ));
        buffer.push_back(0);
    }

    bool config_param(const char* config_file = nullptr)
    {
        if(!util::IO::scan_config(config_file)) {return false;};
        auto db = util::IO::query();

        // Replace default value by queryed configuration if exists
        {
            try { tol  = std::stod(db.at("relTol"));                    } catch(std::exception& e)  { std::cerr << "relTol:\t"     << e.what() << std::endl; }
            try { absTol          = std::stod(db.at("absTol"));         } catch(std::exception& e)  { std::cerr << "absTol:\t"     << e.what() << std::endl; }
            try { L               = std::stod(db.at("L"));              } catch(std::exception& e)  { std::cerr << "L:\t"          << e.what() << std::endl; }
            try { eta             = std::stod(db.at("eta"));            } catch(std::exception& e)  { std::cerr << "eta:\t"        << e.what() << std::endl; }
            try { nk              = std::stoi(db.at("nk"));             } catch(std::exception& e)  { std::cerr << "nk:\t"         << e.what() << std::endl; }
            try { E0              = std::stod(db.at("E0"));             } catch(std::exception& e)  { std::cerr << "E0:\t"         << e.what() << std::endl; }
            try { lambda          = std::stod(db.at("lambda"));         } catch(std::exception& e)  { std::cerr << "lambda:\t"     << e.what() << std::endl; }
            try { randomSeed      = std::stoi(db.at("seed"));           } catch(std::exception& e)  { std::cerr << "seed:\t"       << e.what() << std::endl; }
            try { xoffset         = std::stod(db.at("xoffset"));        } catch(std::exception& e)  { std::cerr << "offset:\t"     << e.what() << std::endl; }
            try { randomSeed     += xoffset;                            } catch(std::exception& e)  { std::cerr << "shuffle:\t"    << e.what() << std::endl; }
            try { tau_inv         = 1.0 / std::stod(db.at("tau"));      } catch(std::exception& e)  { std::cerr << "tau_inv:\t"    << e.what() << std::endl; }
            try { dim[0]          = std::stoi(db.at("numParticles"));   } catch(std::exception& e)  { std::cerr << "dim[0]:\t"     << e.what() << std::endl; }
            try { b               = std::stod(db.at("fallingDistance"));} catch(std::exception& e)  { std::cerr << "falldist:\t"   << e.what() << std::endl; }
            try { g               = std::stod(db.at("g"));              } catch(std::exception& e)  { std::cerr << "grav_acc:\t"   << e.what() << std::endl; }
            try { reducedOutput   = std::stoi(db.at("reducedOutput"));  } catch(std::exception& e)  { std::cerr << "outputFlag:\t" << e.what() << std::endl; }
            try { outputFolder.assign( db.at("outputFolder").c_str() ); } catch(std::exception& e)  { std::cerr << "outputPath:\t" << e.what() << std::endl; }
        }

        // handle passive tracer
        if(std::stod(db.at("tau")) == 0)
            tau_inv = std::numeric_limits<double>::infinity();

        // additional entries
        {
            try { initApproach    = static_cast<initType>(std::stoi(db.at("initApproach"))%NUM_INITTYPE);              } catch(std::exception& e)  { }
        }

        //print_param();

        return true;
    }

    void print_param() const
    {
        std::cout << std::endl;
        printf("**Kinetic Simulation Configuration**\n");
        printf("L = %f, L/eta = %f, nk = %d, E0 = %f, seed = %d\n",L,L/eta, nk, E0, randomSeed);
        printf("nx = %d, ny = %d, nz = %d, initial distance = %f\n",dim[0],dim[1],dim[2],initDist);
        printf("tau_inv = %f, g = %f, fdist = %f, shuffle = %f\n",tau_inv, g, b, xoffset);
        printf("t0 = %f, t_end = %f, dt_max = %f, relTol = %e\n",t0,t_end,get_max_allowed_dt(),tol);
        printf("initial approach = %d\n",initApproach);
        printf("output dir: %s\n", outputFolder.c_str());
        printf("**Kinetic Simulation Configuration End**\n");
        std::cout << std::endl;
    }

};
}

// Generic enum
//#ifndef GENERATE_ENUM_STRINGS
//    #define DECL_ENUM_ELEMENT( element ) element
//    #define BEGIN_ENUM( ENUM_NAME ) typedef enum tag##ENUM_NAME
//    #define END_ENUM( ENUM_NAME ) ENUM_NAME; \
//            char* GetString##ENUM_NAME(enum tag##ENUM_NAME index);
//#else
//    #define DECL_ENUM_ELEMENT( element ) #element
//    #define BEGIN_ENUM( ENUM_NAME ) char* gs_##ENUM_NAME [] =
//    #define END_ENUM( ENUM_NAME ) ; char* GetString##ENUM_NAME(enum \
//            tag##ENUM_NAME index){ return gs_##ENUM_NAME [index]; }
//#endif

#endif
