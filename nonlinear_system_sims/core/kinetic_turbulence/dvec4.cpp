#include "cartesian_field.h"
#include "shader.h"
#include "cmake_source_dir.h"

static Shader dvec4_subroutine[16];

void Vec4::reload_subroutines()
{
    dvec4_subroutine[0].reload_shader_program_from_files(FP("unary_operator/dvec4_add.glsl"));
    dvec4_subroutine[1].reload_shader_program_from_files(FP("unary_operator/dvec4_minus.glsl"));
    dvec4_subroutine[2].reload_shader_program_from_files(FP("unary_operator/dvec4_multiply.glsl"));
    dvec4_subroutine[3].reload_shader_program_from_files(FP("unary_operator/dvec4_divide.glsl"));

    dvec4_subroutine[4].reload_shader_program_from_files(FP("unary_operator/dvec4_add_by_const.glsl"));
    dvec4_subroutine[5].reload_shader_program_from_files(FP("unary_operator/dvec4_minus_by_const.glsl"));
    dvec4_subroutine[6].reload_shader_program_from_files(FP("unary_operator/dvec4_multiply_by_const.glsl"));
    dvec4_subroutine[7].reload_shader_program_from_files(FP("unary_operator/dvec4_divide_by_const.glsl"));
};

Vec4& Vec4::operator= (const Vec4& other) noexcept
{
    if (this != &other) { // same?
        if (other.mesh != mesh) {         // can not perform copy
            std::cout << "Vec4::operator=:only vectors under the same mesh is operable.\n";
            mesh = other.mesh;
        }
        this->copy( other.d_data );
    }
    return *this;
}

// Increment in GPU
Vec4& Vec4::operator+=(const Vec4& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,rhs.d_data); // rhs
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[0].use();
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}

Vec4& Vec4::operator-=(const Vec4& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,rhs.d_data); // rhs
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[1].use();
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}

Vec4& Vec4::operator*=(const Vec4& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,rhs.d_data); // rhs
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[2].use();
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}

Vec4& Vec4::operator/=(const Vec4& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,rhs.d_data); // rhs
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[3].use();
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}

Vec4& Vec4::operator+=(const glm::vec4& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[4].use();
    dvec4_subroutine[4].setVec4("rhs",rhs);
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}

Vec4& Vec4::operator+=(const float& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[4].use();
    dvec4_subroutine[4].setVec4("rhs",rhs,rhs,rhs,rhs);
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}

Vec4& Vec4::operator-=(const glm::vec4& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[5].use();
    dvec4_subroutine[5].setVec4("rhs",rhs);
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}

Vec4& Vec4::operator-=(const float& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[5].use();
    dvec4_subroutine[5].setVec4("rhs",rhs,rhs,rhs,rhs);
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}

Vec4& Vec4::operator*=(const glm::vec4& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[6].use();
    dvec4_subroutine[6].setVec4("rhs",rhs);
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}

Vec4& Vec4::operator*=(const float& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[6].use();
    dvec4_subroutine[6].setVec4("rhs",rhs,rhs,rhs,rhs);
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}

Vec4& Vec4::operator/=(const glm::vec4& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[7].use();
    dvec4_subroutine[7].setVec4("rhs",rhs);
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}

Vec4& Vec4::operator/=(const float& rhs) noexcept
{
    // Bind arrays
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,d_data); // lhs

    // Deploy kernel
    dvec4_subroutine[7].use();
    dvec4_subroutine[7].setVec4("rhs",rhs,rhs,rhs,rhs);
    glDispatchCompute(mesh->grid<0>(), mesh->grid<1>(), mesh->grid<2>());

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(0);

    return *this;
}
