#include "lbm3dgpu.h"
#include <iostream>
#include <cmath>

//GLEW
#define GLEW_STATIC
#include <GL/glew.h>

#include "glm/glm.hpp"
#include "cmake_source_dir.h"

using namespace std;
using namespace glm;

// Instancing
Lattice3DOpenGLInterface* Lattice3DOpenGLInterface::m_singleton = 0;

// Member function definition
Lattice3DOpenGLInterface::Lattice3DOpenGLInterface(int w, int h, int d):nx(w), ny(h), nz(d), nbstep(0), selector(0)
{
    glGenTextures(2, uvwr);
    glGenTextures(2, occl);
    glGenBuffers( 2, f);

    _setupBuffer(0);
    _setupBuffer(1);

    // Visualization
    glGenTextures( 1, &vis);
    glBindTexture(GL_TEXTURE_3D, vis);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage3D( GL_TEXTURE_3D, 0, GL_RGBA32F, nx, ny, nz, 0,
                  GL_RGBA, GL_FLOAT, NULL);

    // Init param and transfer to device (or pass into, will implement later)
    param = new Lattice3DParam<float>(w, h, d);
    _transferParamToDevice();
}

void Lattice3DOpenGLInterface::_compute(){
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, uvwr[1-selector]);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, occl[1-selector]);
    glBindImageTexture(0, uvwr[selector], 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,f[1 - selector]); // f read
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,2,f[selector]); // f write
    glBindImageTexture(3, occl[selector], 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
    glBindImageTexture(4, vis, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA32F);
    // Deploy kernel
    glm::vec3 grid(8,8,8);
    glDispatchCompute((nx+grid.x-1)/grid.x,(ny+grid.y-1)/grid.y,(nz+grid.z-1)/grid.z);

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
}

void Lattice3DOpenGLInterface::_debug() {cout << "Buffer " << selector << " is using." << endl; }

void Lattice3DOpenGLInterface::_bindTextureRead()
{
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, uvwr[1-selector]);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, occl[1-selector]);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_3D, vis);
}
void Lattice3DOpenGLInterface::_bindTextureWrite()
{
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, uvwr[selector]);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_3D, occl[selector]);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_3D, vis);

}

void Lattice3DOpenGLInterface::_setupBuffer(int i)
{
    GLenum internalFormat = GL_RGBA32F;
    GLenum format = GL_RGBA;
    GLenum datatype = GL_FLOAT;

    glBindTexture(GL_TEXTURE_3D, uvwr[0]);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage3D( GL_TEXTURE_3D, 0, internalFormat, nx, ny, nz, 0,
                  format, datatype, NULL);

    glBindTexture(GL_TEXTURE_3D, uvwr[1]);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage3D( GL_TEXTURE_3D, 0, internalFormat, nx, ny, nz, 0,
                  format, datatype, NULL);

    glBindTexture(GL_TEXTURE_3D, occl[0]);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage3D( GL_TEXTURE_3D, 0, internalFormat, nx, ny, nz, 0,
                  format, datatype, NULL);

    glBindTexture(GL_TEXTURE_3D, occl[1]);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage3D( GL_TEXTURE_3D, 0, internalFormat, nx, ny, nz, 0,
                  format, datatype, NULL);

    //Setting up the Shader Storage Buffer Objects in Your C Program
    glBindBuffer( GL_SHADER_STORAGE_BUFFER, f[0] );
    glBufferData( GL_SHADER_STORAGE_BUFFER, nx*ny*nz * sizeof(float)*28, NULL, GL_STATIC_DRAW );

    //Setting up the Shader Storage Buffer Objects in Your C Program
    glBindBuffer( GL_SHADER_STORAGE_BUFFER, f[1] );
    glBufferData( GL_SHADER_STORAGE_BUFFER, nx*ny*nz * sizeof(float)*28, NULL, GL_STATIC_DRAW );

}

void Lattice3DOpenGLInterface::_transferParamToDevice()
{
    if(ubo != NULL) glDeleteBuffers(1, &ubo);
    glGenBuffers(1, &ubo);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(Lattice3DParam<float>) , param, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    // copy data
    //glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    //GLvoid* p = glMapBuffer(GL_UNIFORM_BUFFER, GL_WRITE_ONLY);
    //memcpy(p, &lc, sizeof(lc));
    //glUnmapBuffer(GL_UNIFORM_BUFFER);
    // binding to UBO0
    unsigned int bindingPoint = 0;
    glBindBufferBase(GL_UNIFORM_BUFFER, bindingPoint, ubo); // Bind to UBO0
    //cout << "Binding parameter inputs to ubo["<< uboBindingPoint << "]" << endl;
}

void Lattice3DOpenGLInterface::_dumpAll()
{
    // IO test
    vec4* p = new vec4[nx*ny*nz];
    // left -> right, x-dir first
    glGetTextureImage(occl[_debugGetReadSelector()], 0, GL_RGBA, GL_FLOAT, nx*ny*nz*sizeof(vec4),p);

    // Write to file
    FILE* fid;
    // Dump alpha field
    fid = fopen (FP("result/alpha.dat"),"w");
    if (fid!=NULL)
    {
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 255; i >= 0; i--) {fprintf(fid,"%1.5f ",p[i+nx*j+nx*ny*k].z);}
                fprintf(fid,"\n");
            }
        }
        fclose (fid);
        cout << "Dumped alpha field to alpha.dat" << endl;
    }

    // Dump vel,rho field
    glGetTextureImage(uvwr[_debugGetReadSelector()], 0, GL_RGBA, GL_FLOAT, nx*ny*nz*sizeof(vec4),p);

    fid = fopen (FP("result/u.dat"),"w");
    if (fid!=NULL)
    {
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 255; i >= 0; i--) {fprintf(fid,"%1.5f ",p[i+nx*j+nx*ny*k].x);}
                fprintf(fid,"\n");
            }
        }
        fclose (fid);
        cout << "Dumped u field to u.dat" << endl;
    }
    fid = fopen (FP("result/v.dat"),"w");
    if (fid!=NULL)
    {
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 255; i >= 0; i--) {fprintf(fid,"%1.5f ",p[i+nx*j+nx*ny*k].y);}
                fprintf(fid,"\n");
            }
        }
        fclose (fid);
        cout << "Dumped v field to v.dat" << endl;
    }
    fid = fopen (FP("result/rho.dat"),"w");
    if (fid!=NULL)
    {
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 255; i >= 0; i--) {fprintf(fid,"%1.5f ",p[i+nx*j+nx*ny*k].z);}
                fprintf(fid,"\n");
            }
        }
        fclose (fid);
        cout << "Dumped rho field to rho.dat" << endl;
    }
    fid = fopen (FP("result/f0.dat"),"w");
    if (fid!=NULL)
    {
        for(int k = 0; k < nz; k++)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 255; i >= 0; i--) {fprintf(fid,"%1.5f ",p[i+nx*j+nx*ny*k].w);}
                fprintf(fid,"\n");
            }
        }
        fclose (fid);
        cout << "Dumped f0 field to f0.dat" << endl;
    }
    delete [] p;
}
