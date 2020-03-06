#include "cartesian_field.h"
#include "shader.h"
#include "cmake_source_dir.h"
#include "camera.h"
#include <type_traits>

namespace renderer
{

    static Shader renderer;
    static uint VAO=0,VBO=0;

    void reload_subroutines()
    {
        renderer.reload_shader_program_from_files(FP("renderer/points.vert"),FP("renderer/points.frag"));
    };

    template<>
    void draw<Camera>(const Vec4& target, const Camera& camera)
    {
        // Bind to dummy vertex array
        if(VAO == 0)
        {
            glGenVertexArrays(1, &VAO);
            // fill buffer
            glGenBuffers(1, &VBO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            // link vertex attributes
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            glBindVertexArray(0);
        }
        // Draw points
        renderer.use();
        renderer.setMat4("projectionMatrix", camera.GetFrustumMatrix());
        renderer.setInt("numparticles",target.mesh->N());
        glMemoryBarrier(GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,target.d_data);
        glBindVertexArray(VAO);
        glDrawArraysInstanced(GL_POINTS, 0, 1, target.mesh->N());
        glBindVertexArray(0);
    }

}
