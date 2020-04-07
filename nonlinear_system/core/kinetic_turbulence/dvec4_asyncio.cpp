#include "cartesian_field.h"
#include <future>
#include <chrono>
#include <string>

namespace async_io {
static std::future<void> _threadHandle[4];
static bool first = true;

void dump_data(const std::string&, const std::string&, const std::vector<glm::vec4>&, const Param&);
void dump_snapshot(const std::string&, const std::string&, const std::vector<glm::vec4>&, const Param&);

void wait()
{
    if(first) { return; }

    // blocking threads
    for(int i = 0; i < 2; i++) { _threadHandle[i].wait(); }
}

void dump(const std::string& fieldname, const std::string& path, const std::vector<glm::vec4>& target, const Param& par)
{
    //std::cout << "enter async call\n" << std::endl;
    if (first || _threadHandle[0].wait_for(std::chrono::seconds(0)) == std::future_status::ready)
    {
        _threadHandle[0] = std::async(std::launch::async, dump_data, fieldname, path, target, par);
    }

    // Dump snapshot
    if (first || _threadHandle[1].wait_for(std::chrono::seconds(0)) == std::future_status::ready)
    {
        _threadHandle[1] = std::async(std::launch::async, dump_snapshot, fieldname, path, target, par);
    }
    first = false;
}

// calculate mean particle distance
// group neighbour particles as pairs
// measure x,y,z direction pairs
void dump_data(const std::string& fieldname, const std::string& path, const std::vector<glm::vec4>& data, const Param& par)
{
    // need to replace by .lock file to prevent data broken
    char name[255];
    sprintf(name, (path + "/" + fieldname + ".mean.dat").c_str());

    FILE* fid;
    fid = fopen (name, "a");

    //const char* _fid = string(path).append(name).c_str();
    //std::cout << _fid << std::endl;
    //FILE* fid = fopen (_fid,"w");
    if (fid)
    {
        // 1. Current time
        // 2. Compute mean square distance
        double mean_distance_sqr = 0;

        // nx 0 -> 1, 1 -> 2, 2 -> 3 ... n-1 -> n
        // ny 0 -> 1, 1 -> 2, 2 -> 3 ... n-1 -> n
        uint numpair = (par.dim.x-1)*par.dim.y*par.dim.z;

        // measure x-pair
        for(uint k = 0; k < par.dim.z; k++)
        {
            for(uint j = 0;j < par.dim.y; j++)
            {
                for(uint i = 0;i < par.dim.x-1; i++)
                {
                    // compute x-dir pairs
                    uint left = i + j*par.dim.x + k*par.dim.x*par.dim.y;
                    uint right = i+1 + j*par.dim.x + k*par.dim.x*par.dim.y;
                    mean_distance_sqr += pow(data[left].x - data[right].x, 2.0)
                            + pow(data[left].y - data[right].y, 2.0)
                            + pow(data[left].z - data[right].z, 2.0);
                    //printf("left = %f,%f,%f\n",data[left].x,data[left].y,data[left].z);
                    //printf("right = %f,%f,%f\n",data[right].x,data[right].y,data[right].z);
                }
            }
        }


        // measure y-pair
        numpair += par.dim.x*(par.dim.y-1)*par.dim.z;
        for(uint k = 0; k < par.dim.z; k++)
        {
            for(uint j = 0;j < par.dim.y-1; j++)
            {
                for(uint i = 0;i < par.dim.x; i++)
                {
                    // compute y-dir pairs
                    uint top = i + j*par.dim.x + k*par.dim.x*par.dim.y;
                    uint bottom = i + (j+1)*par.dim.x + k*par.dim.x*par.dim.y;

                    mean_distance_sqr += pow(data[top].x - data[bottom].x, 2.0)
                            + pow(data[top].y - data[bottom].y, 2.0)
                            + pow(data[top].z - data[bottom].z, 2.0);

                }
            }
        }

        // measure z-pair
        numpair += par.dim.x*par.dim.y*(par.dim.z-1);
        for(uint k = 0; k < par.dim.z-1; k++)
        {
            for(uint j = 0;j < par.dim.y; j++)
            {
                for(uint i = 0;i < par.dim.x; i++)
                {
                    // compute z-dir pairs
                    uint front = i + j*par.dim.x + k*par.dim.x*par.dim.y;
                    uint back = i + j*par.dim.x + (k+1)*par.dim.x*par.dim.y;

                    mean_distance_sqr += pow(data[front].x - data[back].x, 2.0)
                            + pow(data[front].y - data[back].y, 2.0)
                            + pow(data[front].z - data[back].z, 2.0);

                }
            }
        }

        mean_distance_sqr /= (double)numpair;

        fprintf(fid,"%f %f\n",par.t,mean_distance_sqr);

        fclose (fid);
        std::cout << "Dumped data to: " << name << " time: " << par.getCurrentTimestep() << std::endl;
    }else {
        std::cout<<"DUMP_STATS::ERROR"<<std::endl;
    }
}

void dump_snapshot(const std::string& fieldname, const std::string& path, const std::vector<glm::vec4>& data, const Param& par)
{

    char name[255];
    sprintf(name, (path + "/" + fieldname + ".%06d.dat").c_str(), par.getIter());

    FILE* fid;
    fid = fopen (name, "w");

    //const char* _fid = string(path).append(name).c_str();
    //std::cout << _fid << std::endl;
    //FILE* fid = fopen (_fid,"w");
    if (fid)
    {
        // 1. Current time
        // 2. Compute mean square distance
        for(int i = 0; i < data.size();i++)
        {
            fprintf(fid,"%f %f %f\n",data[i].x,data[i].y,data[i].z);
        }

        fclose (fid);
        std::cout << "Dumped snapshot to: " << name << std::endl;
    }else {
        std::cout<<"DUMP_SNAPSHOT::ERROR"<<std::endl;
    }
}

}
