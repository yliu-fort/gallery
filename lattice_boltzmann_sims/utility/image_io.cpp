#include "image_io.h"
#include <iostream>

//GLFW
#include <GLFW/glfw3.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include "cmake_source_dir.h"

void PPMWriter(unsigned char *in, char *name, int dimx, int dimy)
{
    int i, j;
    FILE *fp = fopen(name, "wb");
    (void) fprintf(fp, "P6\n%d %d\n255\n", dimx, dimy);
    for(int j = 0; j < dimy; ++j)
    {
        for(i = 0; i < dimx; ++i)
        {
            static unsigned char color[3];
            color[0] = in[3*i + 3*(dimy - 1 - j)*dimx]; // Red
            color[1] = in[3*i + 3*(dimy - 1 - j)*dimx + 1]; // Green
            color[2] = in[3*i + 3*(dimy - 1 - j)*dimx + 2]; // Blue
            (void) fwrite(color, 1, 3, fp);
        }
    }
    (void) fclose(fp);
}

static unsigned char* image = nullptr;
void ImageIO::Save(int w, int h, int imgIndex, bool verbose)
{
    image = new unsigned char[4*w*h];
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, image);
    char buffer[255];
    sprintf(buffer, FP("../result/result.%06d.png"), imgIndex);

    //PPMWriter(image, buffer, w, h);
    stbi_flip_vertically_on_write(true);
    stbi_write_png(buffer, w, h, 3, image, w*3);

    if(verbose) std::cout << "ImageIO:: Saved: "<< buffer <<std::endl;

    //imgIndex++;
    delete [] image;
    //free(image);
}
