#ifndef IMAGE_IO_H
#define IMAGE_IO_H

// triangle table maps same cube vertex index to a list of up to 5 triangles
// which are built from the interpolated edge vertices
// Isosurface require: GLwindow, Camera initialized
// can produce texture with 1-4 channels
namespace ImageIO {
void Save(int w, int h, int imgIndex,bool verbose=false);
}

#endif
