#include "slice.h"
#include <iostream>
#include <cmath>
#include <vector>

//GLEW
#define GLEW_STATIC
#include <GL/glew.h>

//GLFW
#include <GLFW/glfw3.h>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

#include "camera.h"
#include "shader.h"
#include "cmake_source_dir.h"

// ------------------------------------------------------------------------
static unsigned int dummyVAO;
static Shader _shaderHandle;

void Slice::ReloadShader()
{
    _shaderHandle.reload_shader_program_from_files(
                FP("glsl/sl.vert"),FP("glsl/mc.frag"),FP("glsl/sl.geom.glsl"));
}

void Slice::Init()
{
    ReloadShader();

    // Dummy VAO
    glGenVertexArrays(1, &dummyVAO);
}

template<typename T, typename S>
void Slice::Draw(int target, float depth,
                      const T& camera,
                      const S& gridSize)
{

    // Set uniform attributes
    _shaderHandle.use();
    _shaderHandle.setMat4("projectionMatrix", camera.GetFrustumMatrix()*glm::scale(glm::mat4(1.0),
                                              glm::vec3(gridSize.x/(float)gridSize.y,1.0f,gridSize.z/(float)gridSize.y)));
    _shaderHandle.setVec3("viewPos",camera.Position);
    _shaderHandle.setInt("volumeTex",target);
    _shaderHandle.setFloat("depth",(float)depth);
    _shaderHandle.setVec3i("gridSize",gridSize);
    _shaderHandle.setVec3("voxelSize",glm::vec3(1.0f/gridSize.x,1.0f/gridSize.y, 1.0f/gridSize.z));

    // Draw
#ifdef FACE_CULLING
    glEnable(GL_CULL_FACE); // Enable to decrease computing load
    glCullFace(GL_BACK);
#endif
    glBindVertexArray(dummyVAO);
    glDrawArraysInstanced(GL_POINTS,0,1,gridSize.x*gridSize.y*1); // for xy plane
    glBindVertexArray(0);
#ifdef FACE_CULLING
    glDisable(GL_CULL_FACE); // Enable to decrease computing load
#endif
}

template void Slice::Draw<Camera, glm::ivec3>(int, float, const Camera&, const glm::ivec3&);
//template void Slice::Draw<double>(int, double, const glm::mat4&, const glm::ivec3&);
//template void Slice::Draw<int>(int, int, const glm::mat4&, const glm::ivec3&);

// Initialize 3d datafield
//Datafield//
static unsigned int _demoVolumeTex;

// an interesting field function
static float torus(float x, float y, float z)
{
    float R2 = 0.4,a2 = 0.05;
    return (x*x + y*y + z*z + R2 - a2)*(x*x + y*y + z*z + R2 - a2) - 4.0*R2*(x*x + z*z)+1.0;
}
static float tangle(float x, float y, float z)
{
    x *= 3.0f;
    y *= 3.0f;
    z *= 3.0f;
    return (x*x*x*x - 5.0f*x*x +y*y*y*y - 5.0f*y*y +z*z*z*z - 5.0f*z*z + 11.8f) * 0.2f + 0.5f;
}
static float surf_of_evolution(float x, float z, float y)
{
    x *= 3.0f;
    y *= 3.0f;
    z *= 3.0f;
    return (x*x + y*y - (log(z+3.2)*log(z+3.2))-0.02 + 1.0);
}

void Slice::Demo(int nx, int ny, int nz)
{
    //Store the volume data to polygonise
    glGenTextures(1, &_demoVolumeTex);
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, _demoVolumeTex);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    //Generate a distance field to the center of the cube
    glm::vec4 *dataField=new glm::vec4[nx*ny*nz];
    for(int k=0; k<nz; k++)
        for(int j=0; j<ny; j++)
            for(int i=0; i<nx; i++){
                //dataField[i+j*nx+k*nx*ny].x = rand() / double(RAND_MAX);
                //dataField[i+j*nx+k*nx*ny].y = rand() / double(RAND_MAX);
                //dataField[i+j*nx+k*nx*ny].z = rand() / double(RAND_MAX);
                dataField[i+j*nx+k*nx*ny].x=glm::distance(glm::vec3(i, j, k),glm::vec3(64))/64.0f;
                dataField[i+j*nx+k*nx*ny].y=glm::distance(glm::vec3(i, j, k),glm::vec3(64))/64.0f;
                dataField[i+j*nx+k*nx*ny].z=glm::distance(glm::vec3(i, j, k),glm::vec3(64))/64.0f;
                //dataField[i+j*nx+k*nx*ny].w = sin(0.4*i)+sin(0.4*j) + sin(0.4*k);
                //dataField[i+j*nx+k*nx*ny].w = sin(0.4*i + 0.4*j + 0.4*k);
                //dataField[i+j*nx+k*nx*ny].w = sin(0.1*i * 0.1*j * 0.1*k);
                //dataField[i+j*nx+k*nx*ny].w = torus(2*i/(float)nx-1,2*j/(float)ny-1,2*k/(float)nz-1);
                dataField[i+j*nx+k*nx*ny].w = tangle(2*i/(float)nx-1,2*j/(float)ny-1,2*k/(float)nz-1);
                //dataField[i+j*nx+k*nx*ny].w = surf_of_evolution(2*i/(float)nx-1,2*j/(float)ny-1,2*k/(float)nz-1);
                //std::cout << dataField[i+j*nx+k*nx*ny] << std::endl;
            }
    glTexImage3D( GL_TEXTURE_3D, 0, GL_RGBA32F, nx, ny, nz, 0,
                  GL_RGBA, GL_FLOAT, dataField);
    delete [] dataField;
    dataField=NULL;
}
