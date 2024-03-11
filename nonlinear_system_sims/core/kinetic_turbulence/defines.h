#ifndef DEFINES_H
#define DEFINES_H
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <memory>

//GLEW
#define GLEW_STATIC
#include <GL/glew.h>

#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

#include "cmake_source_dir.h"

// Forward declaration for class Vec4
class Vec4;
class Param;

namespace algo
{
    enum simType
    {
        DOUBLE_GYRE=0,
        DOUBLE_GYRE3D=1,
        KINETIC_TURB2D=2,
        KINETIC_TURB3D=3,
        INERTIAL_PARTICLE=4

    };
    void reload_subroutines();
    void ode45(const Vec4&, const Vec4&, const Param&);
    void ode45_init(uint& ubo, const Param&);
}
namespace renderer
{
    void reload_subroutines();
    template<typename T> void draw(const Vec4&, const T&);
}
namespace async_io
{
    void wait();
    void dump(const std::string&, const std::string&, const std::vector<glm::vec4>&, const Param&);
}

#endif
