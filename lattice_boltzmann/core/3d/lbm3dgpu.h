#ifndef LBM2DGPU_H
#define LBM2DGPU_H

#include "param3d.h"

class Lattice3DOpenGLInterface
{
    // Allocate and manage stuff on graphic device with opengl framework
    // 4 fixed binding points for textures
    // 0 for uvrf;
    // 1 for f1_4;
    // 2 for f5_8;
    // 3 for occl;
public:
    static void Create(int w, int h, int d)
    {
        m_singleton = new Lattice3DOpenGLInterface(w,h,d);
    }
    static void Destroy()
    {
        delete m_singleton;
    }
    static void Compute()
    {
        m_singleton->_compute();
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
    static int NZ()
    {
        return m_singleton->nz;
    }

private:
    static Lattice3DOpenGLInterface* m_singleton;
    Lattice3DParam<float>* param;

private:
    // Constructor, Destructor, utility functions
    Lattice3DOpenGLInterface(int w, int h, int d);
    ~Lattice3DOpenGLInterface(){}
    void _compute();
    void _incr(){ selector = (selector+1)%2;nbstep++; }
    void _reset(){ nbstep = 0; }
    void _debug();
    int _debugGetWriteSelector() { return selector; }
    int _debugGetReadSelector() { return 1-selector; }
    void _bindTextureRead();
    void _bindTextureWrite();

private:
    // Buffer management
    unsigned int uvwr[2];
    unsigned int f[2];
    unsigned int occl[2];
    unsigned int vis;
    int nx, ny, nz, nbstep;
    int selector;

    unsigned int ubo = NULL;
    unsigned int padding = NULL;

    void _setupBuffer(int i);
    void _transferParamToDevice();
    void _dumpAll();
};

#endif
