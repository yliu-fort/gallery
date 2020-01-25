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
    static void Create(uint w, uint h, uint d, bool restart = false)
    {
        m_singleton = new Lattice3DOpenGLInterface(w,h,d, restart);
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
    static uint GetTimestep()
    {
        return m_singleton->nbstep;
    }
    static void Reset()
    {
        m_singleton->_reset();
    }
    static void DumpAll()
    {
        m_singleton->_dumpAllVTK();
    }
    static bool Isinitializing()
    {
        return m_singleton->_isInitializing();
    }
    static bool IsInInternalCycle()
    {
        return m_singleton->runtime_param->isInInternalCycle();
    }
    static uint NX()
    {
        return m_singleton->nx;
    }
    static uint NY()
    {
        return m_singleton->ny;
    }
    static uint NZ()
    {
        return m_singleton->nz;
    }
    static void Autosave()
    {
        m_singleton->_dumpAllVTKBinary();
    }

private:
    static Lattice3DOpenGLInterface* m_singleton;
    Lattice3DParam<float>* param;
    Lattice3DParamRuntime<float>* runtime_param;

private:
    // Constructor, Destructor, utility functions
    Lattice3DOpenGLInterface(uint w, uint h, uint d, bool restart = false);
    ~Lattice3DOpenGLInterface();
    void _bind(int layer);
    void _unbind();
    void _compute();
    void _swap(){ selector = (selector+1)%2; }
    void _incr(){nbstep++; }
    void _reset();
    void _debug();
    uint _debugGetWriteSelector() const { return selector; }
    uint _debugGetReadSelector() const { return (selector+1)%2; }
    void _bindTextureRead();
    void _bindTextureWrite();
    bool _isInitializing(){return (nbstep == nbstepRef) && isInitializing;}
    size_t _nElem() const {return (static_cast<size_t>(nx)*static_cast<size_t>(ny)*static_cast<size_t>(nz));}

private:
    // Buffer management
    uint uvwr[2];
    uint f[8];
    uint occl[2];
    uint vis;
    uint nbstepRef;
    uint nbstep;
    uint selector;
    uint nx, ny, nz;
    uint ntask_x, ntask_y, ntask_z, tile;
    uint ubo[2];
    float _padding;

    // host data
    float *h_uvwr, *h_occl, **h_f;
    bool isInitializing; // indicating whether we need initialize simulation when start. false when loaded restart file, true when performed a default initialization

    void _setupBuffer(int i);
    void _transferParamToDeviceUBO0();
    void _transferParamToDeviceUBO1();
    void _dumpAll();
    void _dumpAllVTK(const char* path = nullptr);
    void _dumpAllVTKRect(const char* path = nullptr);
    void _dumpAllVTKBinary(const char* path = nullptr) const;
    void _readAllVTKBinary(const char* path = nullptr);
    void _dumpMetadata(const char* path = nullptr) const;
    void _readMetadata(const char* path = nullptr);
};

#endif
