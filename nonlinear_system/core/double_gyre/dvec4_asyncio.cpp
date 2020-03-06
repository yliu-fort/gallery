#include "cartesian_field.h"
#include <future>
#include <chrono>
#include "cmake_source_dir.h"

namespace async_io {
static std::future<void> _threadHandle[4];
static bool first = true;

void dump_data(const std::vector<glm::vec4>&, const float&, const glm::uvec4&);
void dump_ss(const std::vector<glm::vec4>&);

void wait()
{
    if(first) { return; }

    // blocking threads
    for(int i = 0; i < 2; i++) { _threadHandle[i].wait(); }
}

void dump(const std::vector<glm::vec4>& target, const float& label, const glm::uvec4& label2)
{
    //std::cout << "enter async call\n" << std::endl;
    if (first || _threadHandle[0].wait_for(std::chrono::seconds(0)) == std::future_status::ready)
    {
        _threadHandle[0] = std::async(std::launch::async, dump_data, target, label, label2);
    }

    // Dump snapshot
    if (first || _threadHandle[1].wait_for(std::chrono::seconds(0)) == std::future_status::ready)
    {
        _threadHandle[1] = std::async(std::launch::async, dump_ss, target);
    }
    first = false;
}

void dump_data(const std::vector<glm::vec4>& data, const float& timeLabel, const glm::uvec4& dim)
{
    char name[255];
    sprintf(name, FP("result.dat"));

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
        uint numpair = (dim.x-1)*dim.y;
        for(uint i = 0;i < dim.x-1; i++)
        {
            // compute x-dir pairs
            for(uint j = 0;j < dim.y; j++)
            {
                uint left = i*dim.y + j;
                uint right = (i+1)*dim.y + j;
                mean_distance_sqr += pow(data[left].x - data[right].x, 2.0)
                                   + pow(data[left].y - data[right].y, 2.0);
            }
        }
        mean_distance_sqr /= (double)numpair;

        fprintf(fid,"%f %f\n",timeLabel,mean_distance_sqr);

        fclose (fid);
        std::cout << "Dumped data to: " << name << " time: " << timeLabel << std::endl;
    }else {
        std::cout<<"ERROR"<<std::endl;
    }
}

void dump_ss(const std::vector<glm::vec4>& data)
{
    char name[255];
    sprintf(name, FP("snapshot.dat"));

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
            fprintf(fid,"%f %f\n",data[i].x,data[i].y);
        }

        fclose (fid);
        std::cout << "Dumped snapshot to: " << name << std::endl;
    }else {
        std::cout<<"ERROR"<<std::endl;
    }
}

}
