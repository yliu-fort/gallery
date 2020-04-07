#ifndef CARTESIAN_FIELD_H
#define CARTESIAN_FIELD_H

// Forward declaration for functions in namespace algo
#include "defines.h"
#include "param.h"

// Consider to write a singleton class to manage pointers

// Allocate geometric field
class Cartesian2d
{
public:
    // Construct: allocate memory in cpu and gpu
    Cartesian2d() : d_data(0), dim(0) {}
    Cartesian2d(uint nx, uint ny, uint nz=1, uint gx=16, uint gy=16, uint gz=1) : d_data(0), dim(0)
    {
        // Dimension of the data
        dim = glm::uvec4(nx,ny,nz,nx*ny*nz);

        // Allocate memory on device side
        glGenBuffers(1, &d_data);
        glBindBuffer( GL_SHADER_STORAGE_BUFFER, d_data );
        glBufferData( GL_SHADER_STORAGE_BUFFER, dim.w * sizeof(glm::vec4), NULL, GL_DYNAMIC_DRAW );
        glBindBuffer( GL_SHADER_STORAGE_BUFFER, 0);

        // Compute necessary grid size
        glm::uvec3 grid(gx,gy,gz);
        gridDim = glm::uvec3((dim.x+grid.x-1)/grid.x,(dim.y+grid.y-1)/grid.y,(dim.z+grid.z-1)/grid.z);

        // UBO, transfer data to device side
        glGenBuffers(1, &ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, ubo);
        glBufferData(GL_UNIFORM_BUFFER, sizeof(glm::uvec4) , glm::value_ptr(dim), GL_STATIC_DRAW);
        glBindBuffer(GL_UNIFORM_BUFFER, 0);
        // binding to UBO0
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo); // Bind to UBO0

    }

    // Destroy: release memory
    ~Cartesian2d()
    {
        glDeleteBuffers(1, &d_data);
        glDeleteBuffers(1, &ubo);
    }
    // Sync: data transfer
    void sync()
    {
        // Copy host data to device
        glBindBuffer( GL_SHADER_STORAGE_BUFFER, d_data );
        glBufferData( GL_SHADER_STORAGE_BUFFER, dim.w * sizeof(glm::vec4), &h_data[0], GL_DYNAMIC_DRAW );
        glBindBuffer( GL_SHADER_STORAGE_BUFFER, 0 );
    }

    // transfer data from device to host
    // must wait previous async execution before get a mapped pointer to deploy new io event
    void gather()
    {
        // Blocking if necessary
        async_io::wait();

        // Get mapped pointer
        h_data.resize(dim.w);
        glm::vec4* data = (glm::vec4*)glMapNamedBufferRange(d_data, 0, dim.w, GL_MAP_READ_BIT);

        // Dump device data
        for(int i = 0; i < dim.w; i++){ h_data[i] = data[i]; }

        if(!glUnmapNamedBuffer(d_data)) std::cout <<"unmapping failed\n";

    }

    void dispatch() const noexcept
    {
        glDispatchCompute(gridDim.x, gridDim.y, gridDim.z);
    }

    uint N() const
    {
        return dim.w;
    }

    template<unsigned int Y>
    uint n() const
    {
        static_assert((Y < 3), "dimension must be 0, 1 or 2!");
        return dim[Y];
    }

    uint memsize() const
    {
        return dim.w*sizeof(glm::vec4);
    }

    template<unsigned int Y>
    uint grid() const
    {
        static_assert((Y < 3), "dimension must be 0, 1 or 2!");
        return gridDim[Y];
    }

private:
    std::vector<glm::vec4> h_data;
    uint d_data;
    glm::uvec4 dim; // {nx, ny, nz, nelem}
    glm::uvec3 gridDim;
    uint ubo; // Bind parameters to device side

};

//4-components vector
class Vec4
{
public:
    static void reload_subroutines();

    // Construct: allocate memory in cpu and gpu
    Vec4(std::shared_ptr<Cartesian2d> inmesh) : mesh(inmesh)
    {
        // Allocate memory on device side
        alloc();
    }

    Vec4(std::vector<glm::vec4>& indata, std::shared_ptr<Cartesian2d> inmesh) : mesh(inmesh)
    {
        // Allocate memory on host side
        h_data.resize(mesh->N());
        h_data.assign(indata.begin(),indata.end());

        // Allocate memory on device side
        alloc(&h_data[0]);
    }

    Vec4(const Vec4 &other) : mesh(other.mesh)
    {
        // Allocate memory on device side
        alloc();

        // Copy gpu data
        this->copy( other.d_data );

        printf("Vec4(Vec4&): Copy constructor is called, this could be very expensive.\n");
    }

    // Destroy: release memory
    ~Vec4()
    {
        glDeleteBuffers(1, &d_data);
    }

    // Sync: data transfer
    void sync()
    {
        // Copy host data to device
        glBindBuffer( GL_SHADER_STORAGE_BUFFER, d_data );
        glBufferData( GL_SHADER_STORAGE_BUFFER, mesh->memsize(), &h_data[0], GL_DYNAMIC_DRAW );
        glBindBuffer( GL_SHADER_STORAGE_BUFFER, 0 );
    }

    void gather()
    {
        // Blocking if necessary
        async_io::wait();

        // Get mapped pointer
        h_data.resize(mesh->N());
        glm::vec4* data = (glm::vec4*)glMapNamedBufferRange(d_data, 0, mesh->N(),GL_MAP_READ_BIT);

        // Dump device data
        for(int i = 0; i < mesh->N(); i++){ h_data[i] = data[i]; }

        if(!glUnmapNamedBuffer(d_data)) {std::cout <<"unmapping failed\n";}
    }

    void dispatch() const noexcept
    {
      mesh->dispatch();
    }

    std::vector<glm::vec4>& get_hptr()
    {
        return h_data;
    }

    const uint& get_dptr() const
    {
        return d_data;
    }

    const std::shared_ptr<Cartesian2d> get_meshptr() const
    {
        return mesh;
    }

    // Computing subroutines ?
    // Unary operators
    Vec4& operator= (const Vec4& other) noexcept;
    Vec4& operator+=(const Vec4& rhs) noexcept;
    Vec4& operator-=(const Vec4& rhs) noexcept;
    Vec4& operator*=(const Vec4& rhs) noexcept;
    Vec4& operator/=(const Vec4& rhs) noexcept;

    Vec4& operator+=(const glm::vec4& rhs) noexcept;
    Vec4& operator+=(const float& rhs) noexcept;
    Vec4& operator-=(const glm::vec4& rhs) noexcept;
    Vec4& operator-=(const float& rhs) noexcept;
    Vec4& operator*=(const glm::vec4& rhs) noexcept;
    Vec4& operator*=(const float& rhs) noexcept;
    Vec4& operator/=(const glm::vec4& rhs) noexcept;
    Vec4& operator/=(const float& rhs) noexcept;
    friend Vec4 operator-(Vec4 lhs) noexcept
    {
        lhs *= -1;
        return lhs;
    }

    // Binary operators
    template<typename T>
    friend Vec4 operator+(Vec4 lhs, const T& rhs) noexcept
    {
      lhs += rhs;
      return lhs;
    }
    template<typename T>
    friend Vec4 operator-(Vec4 lhs, const T& rhs) noexcept
    {
      lhs -= rhs;
      return lhs;
    }
    template<typename T>
    friend Vec4 operator*(Vec4 lhs, const T& rhs) noexcept
    {
      lhs *= rhs;
      return lhs;
    }
    template<typename T>
    friend Vec4 operator/(Vec4 lhs, const T& rhs) noexcept
    {
      lhs /= rhs;
      return lhs;
    }

    friend Vec4 operator+(const float& lhs, Vec4 rhs) noexcept
    {
      rhs += lhs;
      return rhs;
    }

    friend Vec4 operator-(const float& lhs, Vec4 rhs) noexcept
    {
      rhs -= lhs;
      rhs *= -1;
      return rhs;
    }

    friend Vec4 operator*(const float& lhs, Vec4 rhs) noexcept
    {
      rhs *= lhs;
      return rhs;
    }

    // Algorithms need forward declaration & friend declaration here

    void report_minmax()
    {
        gather();
        std::vector<float> px;
        std::vector<float> py;
        std::vector<float> pz;
        std::vector<float> pw;
        for(uint i = 0; i < h_data.size(); i++)
        {
            px.push_back(h_data[i].x);
            py.push_back(h_data[i].y);
            pz.push_back(h_data[i].z);
            pw.push_back(h_data[i].w);
        }
        const auto minx = std::min_element(begin(px), end(px));
        const auto maxx = std::max_element(begin(px), end(px));
        const auto miny = std::min_element(begin(py), end(py));
        const auto maxy = std::max_element(begin(py), end(py));
        const auto minz = std::min_element(begin(pz), end(pz));
        const auto maxz = std::max_element(begin(pz), end(pz));
        const auto minw = std::min_element(begin(pw), end(pw));
        const auto maxw = std::max_element(begin(pw), end(pw));
        printf("min = (%f,%f,%f,%f), max = (%f,%f,%f,%f)\n",*minx,*miny,*minz,*minw,*maxx,*maxy,*maxz,*maxw);
    }

    // Debug functions
    #include "dvec4_debug.h"

private:
    std::vector<glm::vec4> h_data;
    uint d_data;
    std::shared_ptr<Cartesian2d> mesh;

    void alloc(glm::vec4* data = NULL)
    {
        glGenBuffers(1, &d_data);
        glBindBuffer( GL_SHADER_STORAGE_BUFFER, d_data );
        glBufferData( GL_SHADER_STORAGE_BUFFER, mesh->memsize(), data, GL_DYNAMIC_DRAW );
        glBindBuffer( GL_SHADER_STORAGE_BUFFER, 0);
    }

    void copy(uint o_data)
    {
        glBindBuffer( GL_COPY_READ_BUFFER, o_data );
        glBindBuffer( GL_COPY_WRITE_BUFFER, d_data );
        glCopyBufferSubData( GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER ,0, 0, mesh->memsize() );
        glBindBuffer( GL_COPY_READ_BUFFER, 0 );
        glBindBuffer( GL_COPY_WRITE_BUFFER, 0 );
    }
};


#endif
