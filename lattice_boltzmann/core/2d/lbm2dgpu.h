#ifndef LBM2DGPU_H
#define LBM2DGPU_H

#include "param2d.h"

class Lattice2DOpenGLInterface
{
    // Allocate and manage stuff on graphic device with opengl framework
    // 4 fixed binding points for textures
    // 0 for uvrf;
    // 1 for f1_4;
    // 2 for f5_8;
    // 3 for occl;
public:
    static void Create(int w, int h)
    {
        m_singleton = new Lattice2DOpenGLInterface(w,h);
    }
    static void Destroy()
    {
        delete m_singleton;
    }
    static void Bind()
    {
        m_singleton->_bind();
    }
    static void Unbind()
    {
        m_singleton->_unbind();
    }
    static void Incr()
    {
        m_singleton->_incr();
    }
    static void BindTexRead()
    {
        m_singleton->_bindTextureRead();
    }
    static void BindTexWrite()
    {
        m_singleton->_bindTextureWrite();
    }
    static int GetTimestep()
    {
        return m_singleton->nbstep;
    }
    static void Reset()
    {
        m_singleton->_reset();
    }
    static void DumpAll()
    {
        m_singleton->_dumpAll();
    }
    static void UpdateParam()
    {
        m_singleton->_transferParamToDevice();
    }
    static int NX()
    {
        return m_singleton->nx;
    }
    static int NY()
    {
        return m_singleton->ny;
    }

private:
    static Lattice2DOpenGLInterface* m_singleton;
    Lattice2DParam<float>* param;

private:
    // Constructor, Destructor, utility functions
    Lattice2DOpenGLInterface(int w, int h);
    ~Lattice2DOpenGLInterface(){}
    void _bind();
    void _unbind();
    void _incr(){ selector = (selector+1)%2;nbstep++; }
    void _reset(){ nbstep = 0; }
    void _debug();
    int _debugGetWriteSelector() { return selector; }
    int _debugGetReadSelector() { return 1-selector; }
    void _bindTextureRead();
    void _bindTextureWrite();

private:
    // Buffers, framebuffer objects management
    unsigned int buffer[2];
    unsigned int uvrf[2];
    unsigned int f1_4[2];
    unsigned int f5_8[2];
    unsigned int occl[2];
    unsigned int rbo[2];
    int nx, ny, nbstep;
    int selector;

    unsigned int ubo = NULL;
    unsigned int padding = NULL;

    void _setupBuffer(int i);
    void _transferParamToDevice();
    void _dumpAll();
};

#endif
