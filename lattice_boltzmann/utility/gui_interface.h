#ifndef GUI_INTERFACE_H
#define GUI_INTERFACE_H

// triangle table maps same cube vertex index to a list of up to 5 triangles
// which are built from the interpolated edge vertices
// Isosurface require: GLwindow, Camera initialized
// can produce texture with 1-4 channels
namespace GuiInterface {
template <typename T> void Init(T* window);
void Demo();
void ReloadShader();
void Finalize();

// Passing 3d texture
void Draw();

}

#endif
