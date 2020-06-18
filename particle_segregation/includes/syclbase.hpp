#pragma once

#ifdef CL_GL_INTEROP
    #ifdef _WIN32
    #  define WINDOWS_LEAN_AND_MEAN
    #  define NOMINMAX
    #  include <windows.h>
    #endif

    // OpenGL Graphics Includes
    #if defined (__APPLE__) || defined(MACOSX)
    //#include <OpenGL/OpenGL.h>
    //#include <GLUT/glut.h>
    #else
    #if defined(__linux__) || defined(UNIX)
    #include <GL/glx.h>
    #endif
    #endif



    #include <CL/cl_gl.h>
    #if defined (__APPLE__) || defined(MACOSX)
    #define GL_SHARING_EXTENSION "cl_APPLE_gl_sharing"
    #else
    #define GL_SHARING_EXTENSION "cl_khr_gl_sharing"
    #endif
#endif

#include <CL/sycl.hpp>
#include <iostream>


class SYCLbase
{
public:
    virtual ~SYCLbase(){}

    SYCLbase()
    {

        auto async_exception_handler = [](cl::sycl::exception_list l) {
            for (auto ep : l) {
                try {
                    std::rethrow_exception(ep);
                } catch (const cl::sycl::exception& e) {
                    std::cout << "Asynchronous exception caught:\n" << e.what();
                }
            }
        };

        try {
            m_device = m_device_selector.select_device().get();
        } catch (cl::sycl::exception& e) {
            std::cerr << e.what() << std::endl;
        }

        try {
            m_platform = m_device_selector.select_device().get_platform().get();
        } catch (cl::sycl::exception& e) {
            std::cerr << e.what() << std::endl;
        }

        //SYCLprint_info();

#ifdef CL_GL_INTEROP

        // Define OS-specific context properties and create the OpenCL context
#if defined (__APPLE__)
        CGLContextObj kCGLContext = CGLGetCurrentContext();
        CGLShareGroupObj kCGLShareGroup = CGLGetShareGroup(kCGLContext);
        cl_context_properties props[] =
        {
            CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE, (cl_context_properties)kCGLShareGroup,
            0
        };
            m_context = clCreateContext(props, 0,0, NULL, NULL, &err);

#else
#if defined(__linux__) || defined(UNIX)
        cl_context_properties props[] =
        {
            CL_GL_CONTEXT_KHR, (cl_context_properties)glXGetCurrentContext(),
            CL_GLX_DISPLAY_KHR, (cl_context_properties)glXGetCurrentDisplay(),
            CL_CONTEXT_PLATFORM, (cl_context_properties)m_platform,
            0
        };
            m_context = clCreateContext(props, 1, &m_device, NULL, NULL, &err);

#else // Win32
        cl_context_properties props[] =
        {
            CL_GL_CONTEXT_KHR, (cl_context_properties)wglGetCurrentContext(),
            CL_WGL_HDC_KHR, (cl_context_properties)wglGetCurrentDC(),
            CL_CONTEXT_PLATFORM, (cl_context_properties)m_platform,
            0
        };
            m_context = clCreateContext(props, 1, &m_device, NULL, NULL, &err);

#endif
#endif

        if(err) {
            std::cerr << "cl-gl interop initialization failed, "
                         "any program relie on this feature will crash...err code=" << err << "\n";
        }

        try {
            m_q = cl::sycl::queue(cl::sycl::context(m_context),
                                  m_device_selector,
                                  async_exception_handler);
        } catch (cl::sycl::exception& e) {
            std::cerr << "cl::sycl::exception: " << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }


#else
        // construct default queue...
        cl_context_properties props[] = {CL_CONTEXT_PLATFORM, (cl_context_properties)m_platform, 0};
        m_context = clCreateContext(props, 1, &m_device, NULL, NULL, &err);
        m_q = cl::sycl::queue(m_device_selector, async_exception_handler);

        if(err) { std::cerr << "default cl launch error code=" << err << std::endl; exit(EXIT_FAILURE); }

#endif
        // execute user-defined initialization stage...
        //dynamic_cast<derived*>(this)->initialize(args...);
    }

    void SYCLprint_info()
    {

        /* Output device and platform information. */
        auto device = cl::sycl::device(m_device);
        auto deviceName = device.template get_info<cl::sycl::info::device::name>();
        std::cout << " Device Name: " << deviceName << std::endl;
        auto platformName =
                device.get_platform().template get_info<cl::sycl::info::platform::name>();
        std::cout << " Platform Name " << platformName << std::endl;

    }

    // interface methods
    //void SYCLpostinitialize(){}
#ifdef CL_GL_INTEROP
    void create_from_gl_buffer(cl_mem& obj, cl_mem_flags flags, cl_GLuint bufobj)
    {
        obj = clCreateFromGLBuffer(m_context, flags, bufobj , &err);
        if(err) { std::cout << "cl_mem create from gl buffer failed!" << err << "\n";  exit(EXIT_FAILURE); }
    }
#endif
    void wait_and_throw(){ m_q.wait_and_throw(); }
    template <typename Func>
    void submit(Func&& func) { m_q.submit(func); }

    cl::sycl::context getContext() { return m_q.get_context(); }
    cl::sycl::queue getQueue() { return m_q; }

    // derived method
    virtual bool initialize(const char* arg = nullptr) = 0;
    virtual void calc(float arg = 0) = 0;
    virtual void finalize(const char* arg = nullptr) = 0;

    virtual bool load_param(const char* arg = nullptr) = 0;
    virtual void incr_shuffle() = 0; // temporary solution for multiple realizations


protected:

    // global queue
    cl::sycl::default_selector m_device_selector;
    cl::sycl::queue m_q;

    cl_device_id m_device;
    cl_platform_id m_platform;
    cl_context m_context;

    cl_int err;
};
