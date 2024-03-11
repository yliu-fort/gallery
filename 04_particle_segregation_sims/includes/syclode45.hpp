#pragma once

#include "syclbase.hpp"
#include "sycl_ode_system.hpp"
#include "param.h"

#define MAX_TRAVEL_TIME (100.0)

template <typename T, int N>
class ode45: public SYCLbase  {

public:
    using vec = cl::sycl::vec<T,N>;

    bool initialize(const char* config_file = nullptr) override;
#ifdef CL_GL_INTEROP
    bool initialize(T* gl_ptr, unsigned int gl_buf, const char* config_file = nullptr );
#endif
    void setup();
    void calc(float tf = MAX_TRAVEL_TIME) override;
    void finalize(const char* filename = nullptr) override;
    void finalize(T* ptr, int stride = 0);

    // temporary solution
    bool load_param(const char* config_file = nullptr) override {
        if(!imported_param.config_param(config_file)) { return false; }
        imported_param.print_param();
        return true;
    }
    void incr_shuffle() override { imported_param.xoffset++; };

protected:
    void configParam();

private:
    size_t array_size;
    algo::Param imported_param;

    cl::sycl::buffer<vec, 1> buf;
    cl::sycl::buffer<T, 1> buf_s;
    cl::sycl::buffer<cl::sycl::vec<T,4>, 1> param;
#ifdef CL_GL_INTEROP
    cl_mem mem_obj;
#endif
};
