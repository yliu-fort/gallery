#include "lbm2dgpu.h"
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
Lattice2DOpenGLInterface* Lattice2DOpenGLInterface::m_singleton = 0;

// Member function definition
Lattice2DOpenGLInterface::Lattice2DOpenGLInterface(int w, int h):nx(w), ny(h), nbstep(0), selector(0)
{
    glGenFramebuffers(2, buffer);
    glGenTextures(2, uvrf);
    glGenTextures(2, f1_4);
    glGenTextures(2, f5_8);
    glGenTextures(2, occl);
    glGenRenderbuffers(2, rbo);

    _setupBuffer(0);
    _setupBuffer(1);

    // Init param and transfer to device (or pass into, will implement later)
    param = new Lattice2DParam<float>(w, h);
    _transferParamToDevice();
}

void Lattice2DOpenGLInterface::_bind(){glBindFramebuffer(GL_FRAMEBUFFER, buffer[selector]);}
void Lattice2DOpenGLInterface::_unbind(){glBindFramebuffer(GL_FRAMEBUFFER, 0);}

void Lattice2DOpenGLInterface::_debug() {cout << "Buffer " << selector << " is using." << endl; }

void Lattice2DOpenGLInterface::_bindTextureRead()
{
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, uvrf[1-selector]);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, f1_4[1-selector]);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, f5_8[1-selector]);
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, occl[1-selector]);
}
void Lattice2DOpenGLInterface::_bindTextureWrite()
{
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, uvrf[selector]);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, f1_4[selector]);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, f5_8[selector]);
    glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_2D, occl[selector]);
}

void Lattice2DOpenGLInterface::_setupBuffer(int i)
{
    GLenum internalFormat = GL_RGBA32F;
    GLenum format = GL_RGBA;
    GLenum datatype = GL_FLOAT;

    unsigned int gBuffer = buffer[i];

    glBindFramebuffer(GL_FRAMEBUFFER, gBuffer);
    unsigned int bufferTexture1 = uvrf[i];
    unsigned int bufferTexture2 = f1_4[i];
    unsigned int bufferTexture3 = f5_8[i];
    unsigned int bufferTexture4 = occl[i];
    // buffer 1 for u, v, rho, f0
    //glGenTextures(1, &bufferTexture1);
    glBindTexture(GL_TEXTURE_2D, bufferTexture1);
    glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, nx, ny, 0, format, datatype, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, bufferTexture1, 0);
    // buffer 2 for f1~f4
    //glGenTextures(1, &bufferTexture2);
    glBindTexture(GL_TEXTURE_2D, bufferTexture2);
    glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, nx, ny, 0, format, datatype, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, bufferTexture2, 0);
    // buffer 3 for f5~f8
    //glGenTextures(1, &bufferTexture3);
    glBindTexture(GL_TEXTURE_2D, bufferTexture3);
    glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, nx, ny, 0, format, datatype, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, bufferTexture3, 0);
    // buffer 4 for boundary
    //glGenTextures(1, &bufferTexture3);
    glBindTexture(GL_TEXTURE_2D, bufferTexture4);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, nx, ny, 0, GL_RGB, GL_FLOAT, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_2D, bufferTexture4, 0);

    // tell OpenGL which color attachments we'll use (of this framebuffer) for rendering
    unsigned int attachments[4] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3 };
    glDrawBuffers(4, attachments);
    // create and attach depth buffer (renderbuffer)
    unsigned int rboDepth = rbo[i];
    //glGenRenderbuffers(1, &rboDepth);
    glBindRenderbuffer(GL_RENDERBUFFER, rboDepth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, nx, ny);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rboDepth);
    // finally check if framebuffer is complete
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        cout << "Framebuffer not complete!" << endl;
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void Lattice2DOpenGLInterface::_transferParamToDevice()
{
    if(ubo != NULL) glDeleteBuffers(1, &ubo);
    glGenBuffers(1, &ubo);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(Lattice2DParam<float>) , param, GL_DYNAMIC_DRAW);
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

void Lattice2DOpenGLInterface::_dumpAll()
{
        // IO test
        vec4* p = new vec4[nx*ny];
        // left -> right, x-dir first
        glGetTextureImage(occl[_debugGetReadSelector()], 0, GL_RGBA, GL_FLOAT, nx*ny*sizeof(vec4),p);

        // Write to file
        FILE* fid;
        // Dump alpha field
        fid = fopen (FP("result/alpha.dat"),"w");
        if (fid!=NULL)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 255; i >= 0; i--) {fprintf(fid,"%1.5f ",p[i+nx*j].z);}
                fprintf(fid,"\n");
            }
            fclose (fid);
            cout << "Dumped alpha field to alpha.dat" << endl;
        }

        // Dump vel,rho field
        glGetTextureImage(uvrf[_debugGetReadSelector()], 0, GL_RGBA, GL_FLOAT, nx*ny*sizeof(vec4),p);

        fid = fopen (FP("result/u.dat"),"w");
        if (fid!=NULL)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 255; i >= 0; i--) {fprintf(fid,"%1.5f ",p[i+nx*j].x);}
                fprintf(fid,"\n");
            }
            fclose (fid);
            cout << "Dumped u field to u.dat" << endl;
        }
        fid = fopen (FP("result/v.dat"),"w");
        if (fid!=NULL)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 255; i >= 0; i--) {fprintf(fid,"%1.5f ",p[i+nx*j].y);}
                fprintf(fid,"\n");
            }
            fclose (fid);
            cout << "Dumped v field to v.dat" << endl;
        }
        fid = fopen (FP("result/rho.dat"),"w");
        if (fid!=NULL)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 255; i >= 0; i--) {fprintf(fid,"%1.5f ",p[i+nx*j].z);}
                fprintf(fid,"\n");
            }
            fclose (fid);
            cout << "Dumped rho field to rho.dat" << endl;
        }
        fid = fopen (FP("result/f0.dat"),"w");
        if (fid!=NULL)
        {
            for(int j = 0; j < ny; j++)
            {
                for(int i = 255; i >= 0; i--) {fprintf(fid,"%1.5f ",p[i+nx*j].w);}
                fprintf(fid,"\n");
            }
            fclose (fid);
            cout << "Dumped f0 field to f0.dat" << endl;
        }
        delete [] p;
}
