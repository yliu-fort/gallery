#ifndef MARCHING_CUBE_H
#define MARCHING_CUBE_H

#include "glm/glm.hpp"

// triangle table maps same cube vertex index to a list of up to 5 triangles
// which are built from the interpolated edge vertices
// Isosurface require: GLwindow, Camera initialized
// can produce texture with 1-4 channels
namespace IsoSurface {
void Init();
void Demo(int nx, int ny, int nz);
void ReloadShader();
void Finalize();

// Passing 3d texture
template<typename T> void Draw(int, T,const glm::mat4&, const glm::ivec3&);

}

#endif
