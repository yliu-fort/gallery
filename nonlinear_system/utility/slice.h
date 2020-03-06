#ifndef SLICE_H
#define SLICE_H

//#include "glm/glm.hpp"

// 3dtexture slice: x, y or z dir, depth specified
namespace Slice {
void Init();
void Demo(int nx, int ny, int nz);
void ReloadShader();
void Finalize();

// Passing 3d texture
template<typename T, typename S>
void Draw(int, float,const T&, const S&);

}

#endif
