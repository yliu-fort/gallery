#ifndef BOUNDING_BOX_H
#define BOUNDING_BOX_H

#include "glm/glm.hpp"

// triangle table maps same cube vertex index to a list of up to 5 triangles
// which are built from the interpolated edge vertices
// Isosurface require: GLwindow, Camera initialized
// can produce texture with 1-4 channels
namespace Boundingbox {
void Init();
void Demo(int nx, int ny, int nz);
void ReloadShader();
void Finalize();

// Passing 3d texture
void Draw(const glm::mat4&);

}

#endif
