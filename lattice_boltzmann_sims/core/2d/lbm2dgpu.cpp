#include "lbm2dgpu.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

//GLEW
#define GLEW_STATIC
#include <GL/glew.h>

#include "glm/glm.hpp"
#include "cmake_source_dir.h"

using namespace std;
using namespace glm;

// Instancing
Lattice2DOpenGLInterface* Lattice2DOpenGLInterface::m_singleton = nullptr;

#define TILE (1024)

Lattice2DOpenGLInterface::~Lattice2DOpenGLInterface()
{
    if(h_uvwr) delete [] h_uvwr;
    if(h_occl) delete [] h_occl;
    if(h_f) delete [] h_f;
    delete runtime_param;
    delete param;
    glDeleteBuffers(2, ubo);
    glDeleteTextures(1, &vis);
    glDeleteBuffers(2, f);
    glDeleteBuffers(2, occl);
    glDeleteBuffers(2, uvwr);
    std::cout << "Lattice2DOpenGLInterface:: instance has been released.\n";
}

// Member function definition
Lattice2DOpenGLInterface::Lattice2DOpenGLInterface(uint w, uint h, bool restart):
    nbstepRef(0),
    nbstep(0),
    selector(0),
    nx(w), ny(h), nz(1),
    tile(TILE),
    h_uvwr(nullptr),
    h_occl(nullptr),
    h_f(nullptr),
    isInitializing(true)
{
    // Read restart file...
    if(restart) _readAllVTKBinary();
    if(!isInitializing) std::cout << "Lattice2DOpenGLInterface:: Restart from autosave.bin..." << std::endl;
    else std::cout << "Lattice2DOpenGLInterface:: Reading autosave file failed, will start from initial..." << std::endl;

    // Allocate device buffers
    glGenBuffers(2, uvwr);
    glGenBuffers(2, occl);
    glGenBuffers(2, f);

    // setup buffers with host data or null
    _setupBuffer(0);
    _setupBuffer(1);

    // Allocate visualization texture
    glGenTextures(1, &vis);
    glBindTexture(GL_TEXTURE_3D, vis);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexImage3D( GL_TEXTURE_3D, 0, GL_RGBA32F, nx, ny, nz, 0,
                  GL_RGBA, GL_FLOAT, nullptr);

    // Calculate tiling
    ntask_x = (nx+tile-1)/tile;
    ntask_y = (ny+tile-1)/tile;
    ntask_z = (nz+tile-1)/tile;

    // Init param and transfer to device (or pass into, will implement later)
    glGenBuffers(2, ubo); // first initialization
    param = new Lattice2DParam<float>(static_cast<uint>(nx), static_cast<uint>(ny));
    _transferParamToDeviceUBO0();

    runtime_param = new Lattice2DParamRuntime<float>(tile, nbstep);
    _transferParamToDeviceUBO1();

    std::cout << "Lattice2DOpenGLInterface:: instance has been created.\n";
}

void Lattice2DOpenGLInterface::_compute(){

    // Transfer runtime param
    _transferParamToDeviceUBO1();

    // Bind buffers
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,uvwr[1-selector]); // uvwr read
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,occl[1-selector]); // occl read

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,2,uvwr[selector]); // uvwr write
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,3,occl[selector]); // occl write

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,4,f[(1 - selector)]); // f read
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER,5,f[selector]); // f write

    //glBindBufferBase(GL_SHADER_STORAGE_BUFFER,6,vis); // vis write
    glBindImageTexture(0, vis, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA32F);

    // Deploy kernel
    glm::vec3 grid(16,16,1);
    glDispatchCompute((TILE+grid.x-1)/grid.x,(TILE+grid.y-1)/grid.y,1);

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);

    // Swap buffer if internal cycle completed
    if(runtime_param->internalIterationCounter(ntask_x, ntask_y, ntask_z, nbstep))
    {
        _swap();

        // Verbose
        //std::cout << "Timestep: " << nbstep << std::endl;
    }
}

void Lattice2DOpenGLInterface::_debug() {cout << "Buffer " << selector << " is using." << endl; }

void Lattice2DOpenGLInterface::_bindTextureRead()
{
    //glBindBufferBase(GL_SHADER_STORAGE_BUFFER,0,vis); // f read
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_3D, vis);
}
void Lattice2DOpenGLInterface::_setupBuffer(int i)
{
    glBindBuffer( GL_SHADER_STORAGE_BUFFER, uvwr[i] );
    glBufferData( GL_SHADER_STORAGE_BUFFER, _nElem() * sizeof(glm::vec4), h_uvwr, GL_DYNAMIC_DRAW );

    glBindBuffer( GL_SHADER_STORAGE_BUFFER, occl[i] );
    glBufferData( GL_SHADER_STORAGE_BUFFER, _nElem() * sizeof(glm::vec4), h_occl, GL_DYNAMIC_DRAW );

    //Setting up the Shader Storage Buffer Objects in Your C Program
    glBindBuffer( GL_SHADER_STORAGE_BUFFER, f[i] );
    glBufferData( GL_SHADER_STORAGE_BUFFER, _nElem() * sizeof(glm::mat3x4), h_f, GL_DYNAMIC_DRAW );

}
// Will set nbstep to nbstepRef (0 when start from beginning and arbitrary number when restart)
// reset internal iteration counting
// overwrite all buffers with host data
void Lattice2DOpenGLInterface::_reset()
{
    nbstep = nbstepRef;
    runtime_param->reset(nbstep);
    _setupBuffer(0);
    _setupBuffer(1);
}
void Lattice2DOpenGLInterface::_transferParamToDeviceUBO0()
{
    //if(ubo[0]) glDeleteBuffers(1, &ubo[0]);
    //glGenBuffers(1, &ubo[0]);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo[0]);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(Lattice2DParam<float>) , param, GL_STATIC_DRAW);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    // binding to UBO0
    glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo[0]); // Bind to UBO0
    //cout << "Binding parameter inputs to ubo["<< uboBindingPoint << "]" << endl;
}

void Lattice2DOpenGLInterface::_transferParamToDeviceUBO1()
{
    //if(ubo[1]) glDeleteBuffers(1, &ubo[1]);
    //glGenBuffers(1, &ubo[1]);
    glBindBuffer(GL_UNIFORM_BUFFER, ubo[1]);
    glBufferData(GL_UNIFORM_BUFFER, sizeof(Lattice2DParamRuntime<float>) , runtime_param, GL_DYNAMIC_DRAW);
    glBindBuffer(GL_UNIFORM_BUFFER, 0);
    // binding to UBO0
    glBindBufferBase(GL_UNIFORM_BUFFER, 1, ubo[1]); // Bind to UBO1
    //cout << "Binding parameter inputs to ubo["<< uboBindingPoint << "]" << endl;
}

//CSV (Comma Separated Variable) files
//CSV files can be read by ParaView, and are a good quick and dirty format.
//This data can be converted into points or structured grids.
//This data is just a number of rows, each row representing a point in space.
//The columns should include X, Y, Z and any other data. An example follows.
//Cut and paste this block of data into a file named test.csv.
//
//x coord, y coord, z coord, scalar
//0, 0, 0, 0
//1, 0, 0, 1
//0, 1, 0, 2
//1, 1, 0, 3
//-0.5, -0.5, 1, 4
//0.5, -0.5, 1, 5
//-0.5, 0.5, 1, 6
//0.5, 0.5, 1, 7
void Lattice2DOpenGLInterface::_dumpAll(const char* path) const
{
    // Make sure we got complete field
    if(runtime_param->isInInternalCycle()) return;

    // IO test
    vec4* p = new vec4[nx*ny*nz];
    // left -> right, x-dir first
    glGetTextureImage(occl[_debugGetReadSelector()], 0, GL_RGBA, GL_FLOAT,
            _nElem()*sizeof(vec4),p);

    // Write to file
    char name[255];
    if(!path) path = FP("../../result/");
    sprintf(name, "result.%06d.vtk", nbstep);
    const char* _fid = string(path).append(name).c_str();

    // Dump alpha field
    FILE* fid = fopen (_fid,"w");
    if (fid)
    {
        for(uint k = 0; k < nz; k++)
        {
            for(uint j = 0; j < ny; j++)
            {
                for(uint i = 0; i < nx; i++) {fprintf(fid,"%1.5f ", double(p[i+nx*j+nx*ny*k].z));}
                fprintf(fid,"\n");
            }
        }
        fclose (fid);
        cout << "Dumped alpha field to alpha.dat" << endl;
    }

    // Dump vel,rho field
    glGetTextureImage(uvwr[_debugGetReadSelector()], 0, GL_RGBA, GL_FLOAT,
            _nElem()*sizeof(vec4),p);

    fid = fopen (FP("result/u.dat"),"w");
    if (fid)
    {
        for(uint k = 0; k < nz; k++)
        {
            for(uint j = 0; j < ny; j++)
            {
                for(uint i = 0; i < nx; i++) {fprintf(fid,"%1.5f ", double(p[i+nx*j+nx*ny*k].x));}
                fprintf(fid,"\n");
            }
        }
        fclose (fid);
        cout << "Dumped u field to u.dat" << endl;
    }
    fid = fopen (FP("result/v.dat"),"w");
    if (fid)
    {
        for(uint k = 0; k < nz; k++)
        {
            for(uint j = 0; j < ny; j++)
            {
                for(uint i = 0; i < nx; i++) {fprintf(fid,"%1.5f ", double(p[i+nx*j+nx*ny*k].y));}
                fprintf(fid,"\n");
            }
        }
        fclose (fid);
        cout << "Dumped v field to v.dat" << endl;
    }
    fid = fopen (FP("result/w.dat"),"w");
    if (fid)
    {
        for(uint k = 0; k < nz; k++)
        {
            for(uint j = 0; j < ny; j++)
            {
                for(uint i = 0; i < nx; i++) {fprintf(fid,"%1.5f ",double(p[i+nx*j+nx*ny*k].z));}
                fprintf(fid,"\n");
            }
        }
        fclose (fid);
        cout << "Dumped w field to w.dat" << endl;
    }
    fid = fopen (FP("result/rho.dat"),"w");
    if (fid)
    {
        for(uint k = 0; k < nz; k++)
        {
            for(uint j = 0; j < ny; j++)
            {
                for(uint i = 0; i < nx; i++) {fprintf(fid,"%1.5f ",double(p[i+nx*j+nx*ny*k].w));}
                fprintf(fid,"\n");
            }
        }
        fclose (fid);
        cout << "Dumped rho field to rho.dat" << endl;
    }
    delete [] p;
}

// todo: switch all sprintf to snprintf
// VTK file format (legacy, support volume rendering)
void Lattice2DOpenGLInterface::_dumpAllVTK(const char* path) const
{
    // Make sure we got complete field
    if(runtime_param->isInInternalCycle()) return;

    // Verbose
    std::cout << "Dumping...";

    // vtk structured point: x first, y second, z last

    // Dump vel,rho field
    vec4* p = (vec4*)glMapNamedBufferRange(uvwr[_debugGetReadSelector()], 0,
            _nElem()*sizeof(vec4),GL_MAP_READ_BIT);

    char name[255];
    sprintf(name, "result.%06d.vtk", nbstep);

    FILE* fid;
    if(!path) fid = fopen (name, "w");
    else fid = fopen (string(path).append(name).c_str(), "w");

    //const char* _fid = string(path).append(name).c_str();
    //std::cout << _fid << std::endl;
    //FILE* fid = fopen (_fid,"w");
    if (fid)
    {
        // Write header
        (void) fprintf(fid, "# vtk DataFile Version 2.0\nLBM 3D GPU\n");
        (void) fprintf(fid, "ASCII\nDATASET STRUCTURED_POINTS\n");
        (void) fprintf(fid, "DIMENSIONS %d %d %d\nASPECT_RATIO %f %f %f\nORIGIN 0 0 0\n",nx,ny,nz, 1.0/ny, 1.0/ny, 1.0/ny);
        (void) fprintf(fid, "POINT_DATA %d\n",nx*ny*nz);

        // Velocity
        (void) fprintf(fid, "VECTORS vel float \n");
        for(uint k = 0; k < nz; k++)
        {
            for(uint j = 0; j < ny; j++)
            {
                for(uint i = 0; i < nx; i++)
                {
                    fprintf(fid,"%f %f %f\n",
                            double(p[i+nx*j+nx*ny*k].x),
                            double(p[i+nx*j+nx*ny*k].y),
                            double(p[i+nx*j+nx*ny*k].z));
                } // x first

            }
        }

        // density
        (void) fprintf(fid, "SCALARS rho float \nLOOKUP_TABLE default\n");
        for(uint k = 0; k < nz; k++)
        {
            for(uint j = 0; j < ny; j++)
            {
                for(uint i = 0; i < nx; i++)
                {
                    fprintf(fid,"%f\n",double(p[i+nx*j+nx*ny*k].w));
                } // x first

            }
        }

        // Dump visual data
        //vec4* q = new vec4[_nElem()];
        //glGetTextureImage(vis, 0, GL_RGBA, GL_FLOAT, _nElem()*sizeof(vec4),q);
        //delete [] q;


        fclose (fid);
        cout << "Dumped to: " << name << endl;
    }else {
        std::cout<<"ERROR"<<std::endl;
    }

    if(!glUnmapNamedBuffer(uvwr[_debugGetReadSelector()])) std::cout <<"unmapping failed\n";

}

// VTK rectilinear file format (legacy, does not support volume rendering)
void Lattice2DOpenGLInterface::_dumpAllVTKRect(const char* path) const
{
    // Make sure we got complete field
    if(runtime_param->isInInternalCycle()) return;

    // IO test
    vec4* p = new vec4[nx*ny*nz];
    // vtk structured point: x first, y second, z last

    // Dump vel,rho field
    glGetTextureImage(vis, 0, GL_RGBA, GL_FLOAT, _nElem()*sizeof(vec4),p);

    char name[255];
    if(!path) path = FP("../../result/");
    sprintf(name, "result.%06d.vtk", nbstep);
    const char* _fid = string(path).append(name).c_str();

    FILE* fid = fopen (_fid,"w");
    if (fid)
    {
        // Write header
        (void) fprintf(fid, "# vtk DataFile Version 2.0\nLBM 3D GPU\n");
        (void) fprintf(fid, "ASCII\nDATASET RECTILINEAR_GRID\n");
        (void) fprintf(fid, "DIMENSIONS %d %d %d\n",nx,ny,nz);

        // Coordinates
        (void) fprintf(fid, "X_COORDINATES %d float\n",nx);
        for(uint i = 0; i < nx; i++){ fprintf(fid,"%f ",double(i)/ny);}
        (void) fprintf(fid, "\nY_COORDINATES %d float\n",ny);
        for(uint i = 0; i < ny; i++){ fprintf(fid,"%f ",double(i)/ny);}
        (void) fprintf(fid, "\nZ_COORDINATES %d float\n",nz);
        for(uint i = 0; i < nz; i++){ fprintf(fid,"%f ",double(i)/ny);}
        fprintf(fid,"\n");

        (void) fprintf(fid, "POINT_DATA %d\n",nx*ny*nz);
        // Velocity
        (void) fprintf(fid, "VECTORS vel float \n");
        for(uint k = 0; k < nz; k++)
        {
            for(uint j = 0; j < ny; j++)
            {
                for(uint i = 0; i < nx; i++)
                {
                    fprintf(fid,"%f %f %f\n",
                            double(p[i+nx*j+nx*ny*k].x),
                            double(p[i+nx*j+nx*ny*k].y),
                            double(p[i+nx*j+nx*ny*k].z));
                } // x first

            }
        }

        // density
        (void) fprintf(fid, "SCALARS rho float \nLOOKUP_TABLE default\n");
        for(uint k = 0; k < nz; k++)
        {
            for(uint j = 0; j < ny; j++)
            {
                for(uint i = 0; i < nx; i++)
                {
                    fprintf(fid,"%f\n",double(p[i+nx*j+nx*ny*k].w));
                } // x first

            }
        }
        fclose (fid);
        cout << "Dump to: " << name << endl;
    }else {
        std::cout<<"ERROR"<<std::endl;
    }

    delete [] p;

}

// todo: switch all sprintf to snprintf
// VTK file format (legacy, support volume rendering)
template <typename T>
void SwapEnd(T& var)
{
    char* varArray = reinterpret_cast<char*>(&var);
    for(long i = 0; i < static_cast<long>(sizeof(var)/2); i++)
        std::swap(varArray[sizeof(var) - 1 - i],varArray[i]);
}

// isfinite: check nan value, if present, refuse autosave and goto error handling
void Lattice2DOpenGLInterface::_dumpAllVTKBinary(const char* path) const
{
    // Make sure we got complete field
    if(runtime_param->isInInternalCycle()) return;

    // Dump metadata
    _dumpMetadata(path);

    // Verbose
    std::cout << "Dumping binary file...";

    // vtk structured point: x first, y second, z last
    char name[255];
    sprintf(name, "autosave.bin.lock");
    char _name[255];
    sprintf(_name, "autosave.bin");

    std::ofstream vtkstream;
    if(!path) vtkstream.open (name, std::ios::out | std::ios::trunc | std::ios::binary);
    else vtkstream.open (string(path).append(name).c_str(), std::ios::out | std::ios::trunc | std::ios::binary);

    if (!vtkstream) std::cout<<"ERROR"<<std::endl;

    char header[255];
    // Write header
    vtkstream << "# vtk DataFile Version 2.0\n";
    vtkstream << "LBM 2D GPU\n";
    vtkstream << "BINARY\n";
    vtkstream << "DATASET STRUCTURED_POINTS\n";
    sprintf(header, "DIMENSIONS %d %d %d\nASPECT_RATIO %f %f %f\nORIGIN 0 0 0\n",nx,ny,nz, 1.0/ny, 1.0/ny, 1.0/ny);
    vtkstream << header;
    header[0] = '\0'; // reset char
    sprintf(header, "POINT_DATA %d\n",nx*ny*nz);
    vtkstream << header;

    // Dump vel,rho field
    vec4* p;
    mat3x4* q;

    p = static_cast<vec4*>(glMapNamedBufferRange(uvwr[_debugGetReadSelector()], 0,
                           _nElem()*sizeof(vec4),GL_MAP_READ_BIT));

    // Velocity
    vtkstream << "VECTORS vel float \n";
    for(uint k = 0; k < nz; k++)
    {
        for(uint j = 0; j < ny; j++)
        {
            for(uint i = 0; i < nx; i++)
            {
                //fprintf(fid,"%f %f %f\n",
                //        double(p[i+nx*j+nx*ny*k].x),
                //        double(p[i+nx*j+nx*ny*k].y),
                //        double(p[i+nx*j+nx*ny*k].z));
                float u(p[i+nx*j+nx*ny*k].x),v(p[i+nx*j+nx*ny*k].y),w(p[i+nx*j+nx*ny*k].z);
                if((!isfinite(u)) || (!isfinite(v)) || (!isfinite(w))) {goto err;}
                SwapEnd(u);
                SwapEnd(v);
                SwapEnd(w);
                vtkstream.write(reinterpret_cast<char*>(&u), sizeof(float));
                vtkstream.write(reinterpret_cast<char*>(&v), sizeof(float));
                vtkstream.write(reinterpret_cast<char*>(&w), sizeof(float));
            } // x first

        }
    }

    // density
    vtkstream << "SCALARS rho float \nLOOKUP_TABLE default\n";
    for(uint k = 0; k < nz; k++)
    {
        for(uint j = 0; j < ny; j++)
        {
            for(uint i = 0; i < nx; i++)
            {
                float rho = p[i+nx*j+nx*ny*k].w;
                if(!isfinite(rho)) {goto err;}
                SwapEnd(rho);
                vtkstream.write(reinterpret_cast<char*>(&rho), sizeof(float));
            } // x first

        }
    }

    if(!glUnmapNamedBuffer(uvwr[_debugGetReadSelector()])) std::cout <<"unmapping failed\n";
    p = nullptr;

    // Dump occl field
    p = static_cast<vec4*>(glMapNamedBufferRange(occl[_debugGetReadSelector()], 0,
                           _nElem()*sizeof(vec4),GL_MAP_READ_BIT));

    vtkstream << "SCALARS occl float 4 \nLOOKUP_TABLE default\n";
    for(uint k = 0; k < nz; k++)
    {
        for(uint j = 0; j < ny; j++)
        {
            for(uint i = 0; i < nx; i++)
            {
                vec4 bc = p[i+nx*j+nx*ny*k];
                for(int i = 0; i < 4; i++)
                {
                    if(!isfinite(bc[i])) {goto err;}
                    SwapEnd(bc[i]);
                }
                vtkstream.write(reinterpret_cast<char*>(&bc), sizeof(vec4));
            } // x first
        }
    }

    if(!glUnmapNamedBuffer(occl[_debugGetReadSelector()])) std::cout <<"unmapping failed\n";
    p = nullptr;

    // Dump f field (convert to mat3x3)
    q = static_cast<mat3x4*>(glMapNamedBufferRange(f[_debugGetReadSelector()], 0,
                             _nElem()*sizeof(mat3x4),GL_MAP_READ_BIT));

    vtkstream << "TENSORS equilibrium float \n";
    for(uint k = 0; k < nz; k++)
    {
        for(uint j = 0; j < ny; j++)
        {
            for(uint i = 0; i < nx; i++)
            {
                mat3 f(q[i+nx*j+nx*ny*k][0][0],q[i+nx*j+nx*ny*k][0][1],q[i+nx*j+nx*ny*k][0][2],
                        q[i+nx*j+nx*ny*k][1][0],q[i+nx*j+nx*ny*k][1][1],q[i+nx*j+nx*ny*k][1][2],
                        q[i+nx*j+nx*ny*k][2][0],q[i+nx*j+nx*ny*k][2][1],q[i+nx*j+nx*ny*k][2][2]);
                for(int i = 0; i < 3; i++)
                {
                    for(int j = 0; j < 3; j++)
                    {
                        if(!isfinite(f[i][j])) {goto err;}
                        SwapEnd(f[i][j]);
                    }
                }
                vtkstream.write(reinterpret_cast<char*>(&f), sizeof(mat3));
            } // x first
        }
    }

    if(!glUnmapNamedBuffer(f[_debugGetReadSelector()])) std::cout <<"unmapping failed\n";
    q = nullptr;

    // Dump boundary field
    //vec4* q = new vec4[_nElem()];
    //glGetTextureImage(vis, 0, GL_RGBA, GL_FLOAT, _nElem()*sizeof(vec4),q);
    //delete [] q;

    vtkstream.close();

    // Rename to overwrite old file
    if(!path) rename (name, _name);
    else rename (string(path).append(name).c_str(), string(path).append(_name).c_str());

    cout << "to: " << _name << endl;
    return;

err:
    vtkstream.close();
    cout << "autosave interrupted because nan/inf detected in fetched number." << endl;
    return;

}

void Lattice2DOpenGLInterface::_readAllVTKBinary(const char* path)
{
    // Make sure we got complete field
    //if(runtime_param->isInInternalCycle()) return;
    // Will only call in init stage

    // Dump metadata
    _readMetadata(path);

    // Verbose
    std::cout << "Reading binary file...";

    // vtk structured point: x first, y second, z last
    char name[255];
    sprintf(name, "autosave.bin");

    std::ifstream vtkstream;
    if(!path) vtkstream.open (name, std::ios::in | std::ios::binary);
    else vtkstream.open (string(path).append(name).c_str(), std::ios::in | std::ios::binary);

    if (!vtkstream) {std::cout<<"ERROR"<<std::endl;return;}

    char header[255];
    // Write header
    for(int i = 0; i < 8; i++)
    {
        vtkstream.getline(header,255);
        std::cout << header << std::endl;
    }

    // Dump vel,rho field
    vec4* p1 = new vec4[_nElem()];

    // Velocity
    //vtkstream << "VECTORS vel float \n";
    vtkstream.getline(header,255);
    std::cout << header << std::endl;
    for(uint k = 0; k < nz; k++)
    {
        for(uint j = 0; j < ny; j++)
        {
            for(uint i = 0; i < nx; i++)
            {
                float u,v,w;
                vtkstream.read(reinterpret_cast<char*>(&u), sizeof(float));
                vtkstream.read(reinterpret_cast<char*>(&v), sizeof(float));
                vtkstream.read(reinterpret_cast<char*>(&w), sizeof(float));
                SwapEnd(u);
                SwapEnd(v);
                SwapEnd(w);
                p1[i+nx*j+nx*ny*k].x = u;
                p1[i+nx*j+nx*ny*k].y = v;
                p1[i+nx*j+nx*ny*k].z = w;
            } // x first

        }
    }

    // density
    //vtkstream << "SCALARS rho float \nLOOKUP_TABLE default\n";
    vtkstream.getline(header,255);
    std::cout << header << std::endl;
    vtkstream.getline(header,255);
    std::cout << header << std::endl;
    for(uint k = 0; k < nz; k++)
    {
        for(uint j = 0; j < ny; j++)
        {
            for(uint i = 0; i < nx; i++)
            {
                float rho;
                vtkstream.read(reinterpret_cast<char*>(&rho), sizeof(float));
                SwapEnd(rho);
                p1[i+nx*j+nx*ny*k].w = rho;
            } // x first

        }
    }

    // Read occl field
    vec4* p2 = new vec4[_nElem()];

    //vtkstream << "SCALARS occl float 4 \nLOOKUP_TABLE default\n";
    vtkstream.getline(header,255);
    std::cout << header << std::endl;
    vtkstream.getline(header,255);
    std::cout << header << std::endl;
    for(uint k = 0; k < nz; k++)
    {
        for(uint j = 0; j < ny; j++)
        {
            for(uint i = 0; i < nx; i++)
            {
                vec4 bc;
                for(int s = 0; s < 4; s++)
                {
                    vtkstream.read(reinterpret_cast<char*>(&bc[s]), sizeof(float));
                    SwapEnd(bc[s]);
                }
                p2[i+nx*j+nx*ny*k] = bc;
            } // x first
        }
    }

    // Dump f field (convert to mat3x3)
    mat3x4* q = new mat3x4[_nElem()];

    //vtkstream << "TENSORS equilibrium float \n";
    vtkstream.getline(header,255);
    std::cout << header << std::endl;
    for(uint k = 0; k < nz; k++)
    {
        for(uint j = 0; j < ny; j++)
        {
            for(uint i = 0; i < nx; i++)
            {
                mat3x4 f(0.0f);
                for(int s = 0; s < 3; s++)
                {
                    for(int t = 0; t < 3; t++)
                    {
                        vtkstream.read(reinterpret_cast<char*>(&(f[s][t])), sizeof(float));
                        SwapEnd(f[s][t]);
                    }
                }
                //vtkstream.write(reinterpret_cast<char*>(&f), sizeof(mat3));
                q[i+nx*j+nx*ny*k] = f;
            } // x first
        }
    }

    vtkstream.close();

    // Convert type of pointer
    if(h_uvwr) delete [] h_uvwr; // prevent memory leak
    if(h_occl) delete [] h_occl;
    if(h_f) delete [] h_f;

    h_uvwr = (float *)p1;
    h_occl = (float *)p2;
    h_f = (float *)q;

    // set intialization flag to false
    // we no longer need default initialization when start or reset
    isInitializing = false;

    cout << "ok: " << name << endl;

}
// Dumping essential metadata for recovering simulation
// currently include nx, ny, nz, timestep
void Lattice2DOpenGLInterface::_dumpMetadata(const char* path) const
{
    // Verbose
    std::cout << "Dumping meta-data...";

    // vtk structured point: x first, y second, z last
    char name[255];
    sprintf(name, "autosave.meta");

    std::ofstream vtkstream;
    if(!path) vtkstream.open (name, std::ios::out | std::ios::trunc | std::ios::binary);
    else vtkstream.open (string(path).append(name).c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
    if (!vtkstream) std::cout<<"ERROR"<<std::endl;

    // Write header
    vtkstream << "# METADATA v1.0\n";
    uint _nx(nx),_ny(ny),_nz(nz),_nbstep(nbstep);
    SwapEnd(_nx);SwapEnd(_ny);
    SwapEnd(_nz);SwapEnd(_nbstep);
    vtkstream.write(reinterpret_cast<char*>(&_nx), sizeof(float));
    vtkstream.write(reinterpret_cast<char*>(&_ny), sizeof(float));
    vtkstream.write(reinterpret_cast<char*>(&_nz), sizeof(float));
    vtkstream.write(reinterpret_cast<char*>(&_nbstep), sizeof(float));

    vtkstream.close();
    cout << "to: " << name << endl;

    // Verbose
    std::cout << "Meta-data::Dim:"<< nx << "x" << ny << "x" << nz << "::" << nx*ny*nz << std::endl;
    std::cout << "Meta-data::nbstep: " << nbstep << std::endl;

}

void Lattice2DOpenGLInterface::_readMetadata(const char* path)
{
    // Verbose
    std::cout << "Reading meta-data...";

    // vtk structured point: x first, y second, z last
    char name[255];
    sprintf(name, "autosave.meta");

    std::ifstream vtkstream;
    if(!path) vtkstream.open (name, std::ios::in | std::ios::binary);
    else vtkstream.open (string(path).append(name).c_str(), std::ios::in | std::ios::binary);
    if (!vtkstream) {std::cout<<"ERROR"<<std::endl;return;}

    char header[255];
    // Write header
    vtkstream.getline(header,255);
    std::cout << header << std::endl;
    uint _nx,_ny,_nz,_nbstep;
    vtkstream.read(reinterpret_cast<char*>(&_nx), sizeof(uint));
    vtkstream.read(reinterpret_cast<char*>(&_ny), sizeof(uint));
    vtkstream.read(reinterpret_cast<char*>(&_nz), sizeof(uint));
    vtkstream.read(reinterpret_cast<char*>(&_nbstep), sizeof(uint));
    SwapEnd(_nx);SwapEnd(_ny);
    SwapEnd(_nz);SwapEnd(_nbstep);

    // Assign class member...
    nx = _nx; ny = _ny;
    nz = _nz;
    nbstepRef = _nbstep;
    nbstep = _nbstep;

    vtkstream.close();
    cout << "to: " << name << endl;

    // Verbose
    std::cout << "Meta-data::Dim:"<< nx << "x" << ny << "x" << nz << "::" << nx*ny*nz << std::endl;
    std::cout << "Meta-data::nbstep: " << nbstep << std::endl;

}
