#include <iostream>
#include <cmath>

//GLEW
#define GLEW_STATIC
#include <GL/glew.h>

//GLFW
#include <GLFW/glfw3.h>
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"

#include "cmake_source_dir.h"
#include "shader.h"
#include "camera.h"
#include "colormap.h"
#include "filesystemmonitor.h"
#include "isosurface.h"

// settings
static int SCR_WIDTH  = 800;
static int SCR_HEIGHT = 600;

static int ny = (32);
static int nx = (32);
static int nz = (32);

// camera
Camera camera = Camera(glm::vec3(0.0f, 0.0f, 5.0f), (float)SCR_WIDTH/SCR_HEIGHT);

static float lastX = SCR_WIDTH / 2.0f;
static float lastY = SCR_HEIGHT / 2.0f;
static bool firstMouse = true;

// timing
static float deltaTime = 0.0f;
static float lastFrame = 0.0f;

// fps recording
static float lastFpsCountFrame = 0;
static int frameCount = 0;

// Pre-declaration
GLFWwindow* initGL(int w, int h);
void processInput(GLFWwindow *window);
void countAndDisplayFps(GLFWwindow*);
unsigned int loadTexture(const char *path, bool gammaCorrection);
void renderQuad();
void renderCube();

// Display texture
static int selectTexture = 0;
static int renderType = 0;
static float slidebar = 0.02f;

// simulation related parameter
static bool initFluid = true;

// Initial customized shaders...
// File path must be modified before compiling the script.
static Shader fluid3DInitShader;
static Shader fluid3DComputeShader;
static Shader postProcessShader;

int main()
{
#if defined(__linux__)
    setenv ("DISPLAY", ":0", 0);
#endif

    // Initialize a window
    GLFWwindow* window = initGL(SCR_WIDTH, SCR_HEIGHT);
    printf("Initial glwindow...\n");
    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // Read shader
    fluid3DInitShader.reload_shader_program_from_files(FP("init3d.compute.glsl"));
   
    // For auto-reloading
    FileSystemMonitor::Init(SRC_PATH);

    // Gen colormap
    Colormap::Viridis();

    // Isosurface
    IsoSurface::Init();
    //IsoSurface::Demo(nx,ny,nz);

    // configure g-buffer framebuffer
    // ------------------------------
    //Lattice3DOpenGLInterface::Create(nx, ny, nz);

    //Store the volume data to polygonise
    unsigned int _demoVolumeTex;
    glGenTextures(1, &_demoVolumeTex);
    glActiveTexture(GL_TEXTURE0);
    glEnable(GL_TEXTURE_3D);
    glBindTexture(GL_TEXTURE_3D, _demoVolumeTex);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    //Generate a distance field to the center of the cube
    glm::vec4 *dataField=new glm::vec4[nx*ny*nz];
    for(int k=0; k<nz; k++)
        for(int j=0; j<ny; j++)
            for(int i=0; i<nx; i++){
                dataField[i+j*nx+k*nx*ny].x = rand() / double(RAND_MAX);
                dataField[i+j*nx+k*nx*ny].y = rand() / double(RAND_MAX);
                dataField[i+j*nx+k*nx*ny].z = rand() / double(RAND_MAX);
                dataField[i+j*nx+k*nx*ny].w=glm::distance(glm::vec3(i, j, k),glm::vec3(16))/16.0f;
                //dataField[i+j*nx+k*nx*ny].w = sin(0.4*i)+sin(0.4*j) + sin(0.4*k);
                //dataField[i+j*nx+k*nx*ny].w = sin(0.4*i + 0.4*j + 0.4*k);
                //dataField[i+j*nx+k*nx*ny].w = sin(0.1*i * 0.1*j * 0.1*k);
                //dataField[i+j*nx+k*nx*ny].w = torus(2*i/(float)nx-1,2*j/(float)ny-1,2*k/(float)nz-1);
                //dataField[i+j*nx+k*nx*ny].w = tangle(2*i/(float)nx-1,2*j/(float)ny-1,2*k/(float)nz-1);
                //dataField[i+j*nx+k*nx*ny].w = surf_of_evolution(2*i/(float)nx-1,2*j/(float)ny-1,2*k/(float)nz-1);
                //std::cout << dataField[i+j*nx+k*nx*ny] << std::endl;
            }
    glTexImage3D( GL_TEXTURE_3D, 0, GL_RGBA32F, nx, ny, nz, 0,
                  GL_RGBA, GL_FLOAT, dataField);
    delete [] dataField;
    dataField=NULL;

    unsigned int fBuf;
    glGenBuffers( 1, &fBuf);
    //Setting up the Shader Storage Buffer Objects in Your C Program
    // This is interleaved storage pattern
    glBindBuffer( GL_SHADER_STORAGE_BUFFER, fBuf );
    glBufferData( GL_SHADER_STORAGE_BUFFER, nx*ny*nz * sizeof(double)*(2), NULL, GL_STATIC_DRAW );
    GLint bufMask = GL_MAP_WRITE_BIT | GL_MAP_INVALIDATE_BUFFER_BIT ; // the invalidate makes a big difference when re-writing
    double* val = (double*)glMapBufferRange(GL_SHADER_STORAGE_BUFFER,0,
                                          nx*ny*nz * sizeof(double)*(2),bufMask);
    for(int k=0; k<nz; k++)
        for(int j=0; j<ny; j++)
            for(int i=0; i<nx; i++){
                int ind = i+j*nx+k*nx*ny;
                int scale = 1;
                val[2*ind]=glm::distance(glm::dvec3(i, j, k),glm::dvec3(8*scale))/(8.0*scale);
                //val[5*ind+0].x=glm::distance(glm::dvec3(i, j, k),glm::dvec3(8*scale))/(8.0*scale);
                //val[5*ind+1].x=glm::distance(glm::dvec3(i, j, k),glm::dvec3(16*scale))/(16.0*scale);
                //val[5*ind+2].x=glm::distance(glm::dvec3(i, j, k),glm::dvec3(32*scale))/(32.0*scale);
                //val[5*ind+3].x=glm::distance(glm::dvec3(i, j, k),glm::dvec3(64*scale))/(64.0*scale);
                //val[5*ind+4].x=glm::distance(glm::vec3(i, j, k),glm::vec3(64*scale))/(64.0f*scale);
            }
    glUnmapBuffer( GL_SHADER_STORAGE_BUFFER );

    while( !glfwWindowShouldClose( window ) )
    {
        // per-frame time logic
        // --------------------
        countAndDisplayFps(window);

        // input
        // -----
        glfwGetFramebufferSize(window, &SCR_WIDTH, &SCR_HEIGHT);
        processInput(window);

        { // launch compute shaders
          if(FileSystemMonitor::Update()) IsoSurface::ReloadShader();
          fluid3DInitShader.use();
          fluid3DInitShader.setInt("img_output", 0);
          fluid3DInitShader.setVec3i("gridSize", glm::ivec3(nx, ny, nz));
          fluid3DInitShader.setFloat("time", glfwGetTime());
          glBindImageTexture(0, _demoVolumeTex, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA32F);
          glBindBufferBase(GL_SHADER_STORAGE_BUFFER,1,fBuf);
          glDispatchCompute(nx,ny,nz);

          // make sure writing to image has finished before read
          glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);
          glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        }

        {
            // 2. rendering pass
            glViewport(0,0,SCR_WIDTH, SCR_HEIGHT);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            //postProcessShader.use();
            //Colormap::Bind();
            //Lattice3DOpenGLInterface::BindTexRead();
            //renderQuad();
            IsoSurface::Draw(0, slidebar+1.0f,
                             camera,
                             glm::ivec3(nx,ny,nz));
        }

        //if(Lattice3DOpenGLInterface::GetTimestep() == 15000) Lattice3DOpenGLInterface::DumpAll();

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate( );

    return 0;
}

void countAndDisplayFps(GLFWwindow* window)
{
    float currentFrame = glfwGetTime();
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;

    frameCount++;
    if(glfwGetTime() - lastFpsCountFrame > 1.0f)
    {
        /*std::cout << "Current fps: "
                  << frameCount/(glfwGetTime() - lastFpsCountFrame)
                  << " runtime:"
                  << glfwGetTime()
                  << std::endl;*/ // deprecated

        char title [256];
        title [255] = '\0';

        snprintf ( title, 255,
                   "FAST+ARB DEMO - FPS: %4.2f | runtime: %.0fs ",
                   frameCount/(glfwGetTime() - lastFpsCountFrame), glfwGetTime() );
        glfwSetWindowTitle(window, title);

        frameCount = 0;
        lastFpsCountFrame = glfwGetTime();
    }
}

// renderCube() renders a 1x1 3D cube in NDC.
// -------------------------------------------------
unsigned int cubeVAO = 0;
unsigned int cubeVBO = 0;
void renderCube()
{
    // initialize (if necessary)
    if (cubeVAO == 0)
    {
        float vertices[] = {
            // back face
            -1.0f, -1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 0.0f, 0.0f, // bottom-left
            1.0f,  1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 1.0f, 1.0f, // top-right
            1.0f, -1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 1.0f, 0.0f, // bottom-right
            1.0f,  1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 1.0f, 1.0f, // top-right
            -1.0f, -1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 0.0f, 0.0f, // bottom-left
            -1.0f,  1.0f, -1.0f,  0.0f,  0.0f, -1.0f, 0.0f, 1.0f, // top-left
            // front face
            -1.0f, -1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 0.0f, 0.0f, // bottom-left
            1.0f, -1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 1.0f, 0.0f, // bottom-right
            1.0f,  1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 1.0f, 1.0f, // top-right
            1.0f,  1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 1.0f, 1.0f, // top-right
            -1.0f,  1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 0.0f, 1.0f, // top-left
            -1.0f, -1.0f,  1.0f,  0.0f,  0.0f,  1.0f, 0.0f, 0.0f, // bottom-left
            // left face
            -1.0f,  1.0f,  1.0f, -1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-right
            -1.0f,  1.0f, -1.0f, -1.0f,  0.0f,  0.0f, 1.0f, 1.0f, // top-left
            -1.0f, -1.0f, -1.0f, -1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-left
            -1.0f, -1.0f, -1.0f, -1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-left
            -1.0f, -1.0f,  1.0f, -1.0f,  0.0f,  0.0f, 0.0f, 0.0f, // bottom-right
            -1.0f,  1.0f,  1.0f, -1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-right
            // right face
            1.0f,  1.0f,  1.0f,  1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-left
            1.0f, -1.0f, -1.0f,  1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-right
            1.0f,  1.0f, -1.0f,  1.0f,  0.0f,  0.0f, 1.0f, 1.0f, // top-right
            1.0f, -1.0f, -1.0f,  1.0f,  0.0f,  0.0f, 0.0f, 1.0f, // bottom-right
            1.0f,  1.0f,  1.0f,  1.0f,  0.0f,  0.0f, 1.0f, 0.0f, // top-left
            1.0f, -1.0f,  1.0f,  1.0f,  0.0f,  0.0f, 0.0f, 0.0f, // bottom-left
            // bottom face
            -1.0f, -1.0f, -1.0f,  0.0f, -1.0f,  0.0f, 0.0f, 1.0f, // top-right
            1.0f, -1.0f, -1.0f,  0.0f, -1.0f,  0.0f, 1.0f, 1.0f, // top-left
            1.0f, -1.0f,  1.0f,  0.0f, -1.0f,  0.0f, 1.0f, 0.0f, // bottom-left
            1.0f, -1.0f,  1.0f,  0.0f, -1.0f,  0.0f, 1.0f, 0.0f, // bottom-left
            -1.0f, -1.0f,  1.0f,  0.0f, -1.0f,  0.0f, 0.0f, 0.0f, // bottom-right
            -1.0f, -1.0f, -1.0f,  0.0f, -1.0f,  0.0f, 0.0f, 1.0f, // top-right
            // top face
            -1.0f,  1.0f, -1.0f,  0.0f,  1.0f,  0.0f, 0.0f, 1.0f, // top-left
            1.0f,  1.0f , 1.0f,  0.0f,  1.0f,  0.0f, 1.0f, 0.0f, // bottom-right
            1.0f,  1.0f, -1.0f,  0.0f,  1.0f,  0.0f, 1.0f, 1.0f, // top-right
            1.0f,  1.0f,  1.0f,  0.0f,  1.0f,  0.0f, 1.0f, 0.0f, // bottom-right
            -1.0f,  1.0f, -1.0f,  0.0f,  1.0f,  0.0f, 0.0f, 1.0f, // top-left
            -1.0f,  1.0f,  1.0f,  0.0f,  1.0f,  0.0f, 0.0f, 0.0f  // bottom-left
        };
        glGenVertexArrays(1, &cubeVAO);
        glGenBuffers(1, &cubeVBO);
        // fill buffer
        glBindBuffer(GL_ARRAY_BUFFER, cubeVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
        // link vertex attributes
        glBindVertexArray(cubeVAO);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }
    // render Cube
    glBindVertexArray(cubeVAO);
    glDrawArrays(GL_TRIANGLES, 0, 36);
    glBindVertexArray(0);
}


// renderQuad() renders a 1x1 XY quad in NDC
// -----------------------------------------
unsigned int quadVAO = 0;
unsigned int quadVBO;
void renderQuad()
{
    if (quadVAO == 0)
    {
        float quadVertices[] = {
            // positions        // texture Coords
            -1.0f,  1.0f, 0.0f, 0.0f, 1.0f,
            -1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
            1.0f,  1.0f, 0.0f, 1.0f, 1.0f,
            1.0f, -1.0f, 0.0f, 1.0f, 0.0f,
        };
        // setup plane VAO
        glGenVertexArrays(1, &quadVAO);
        glGenBuffers(1, &quadVBO);
        glBindVertexArray(quadVAO);
        glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
        glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
    }
    glBindVertexArray(quadVAO);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    glBindVertexArray(0);
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
static int key_1_old_state = GLFW_RELEASE;
static int key_2_old_state = GLFW_RELEASE;
static int key_space_old_state = GLFW_RELEASE;
static int key_i_old_state = GLFW_RELEASE;
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);
    
    if ((glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS) && (key_1_old_state == GLFW_RELEASE))
    {
        selectTexture = (++selectTexture)%4;
        std::cout << "Select texture " << selectTexture << std::endl;
    }
    key_1_old_state = glfwGetKey(window, GLFW_KEY_1);
    if ((glfwGetKey(window, GLFW_KEY_2) == GLFW_PRESS) && (key_2_old_state == GLFW_RELEASE))
    {
        renderType = (++renderType)%3;
        std::cout << "Select render mode " << renderType << std::endl;
    }
    key_2_old_state = glfwGetKey(window, GLFW_KEY_2);
    if ((glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) && (key_space_old_state == GLFW_RELEASE))
    {
        initFluid = true;
    }
    key_space_old_state = glfwGetKey(window, GLFW_KEY_SPACE);

    if (glfwGetKey(window, GLFW_KEY_R) == GLFW_PRESS)
        slidebar += deltaTime;
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS)
        slidebar -= deltaTime;

    if ( (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) && (key_i_old_state == GLFW_RELEASE) ) {
         std::cout << "Shader reloaded." << std::endl;
    }
    key_i_old_state = glfwGetKey(window, GLFW_KEY_I);

}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int w, int h)
{
    glViewport(0, 0, w, h);
    camera.updateAspect((float)w / (float)h);
}

bool mouse_button_right = false;
// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos; // reversed since y-coordinates go from bottom to top

    lastX = xpos;
    lastY = ypos;

    if(mouse_button_right)
        camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
    {mouse_button_right = true;return;}
    mouse_button_right = false;
}

GLFWwindow* initGL(int w, int h)
{
    // Initialise GLFW
    if( !glfwInit() )
    {
        fprintf( stderr, "Failed to initialize GLFW\n" );
        getchar();
        EXIT_FAILURE;
    }
    
    //glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    // Open a window and create its OpenGL context
    GLFWwindow* window = glfwCreateWindow( w, h, "Demo", NULL, NULL);
    if( window == NULL ){
        fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
        getchar();
        glfwTerminate();
        EXIT_FAILURE;
    }
    glfwMakeContextCurrent(window);
    
    // Initialize GLEW
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        getchar();
        glfwTerminate();
        EXIT_FAILURE;
    }
    
    // Initial viewport
    glViewport(0, 0, w, h);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    
    // Query infomation
    int nrAttributes;
    glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &nrAttributes);
    std::cout << "Maximum nr of vertex attributes supported: " << nrAttributes << std::endl;
    std::cout << "Hardware: " <<glGetString(GL_RENDERER) << std::endl;
    glGetIntegerv(GL_MAX_DRAW_BUFFERS, &nrAttributes);
    std::cout << "Maximum nr of color attachments supported: " << nrAttributes << std::endl;
    glGetIntegerv(GL_MAX_FRAGMENT_IMAGE_UNIFORMS, &nrAttributes);
    std::cout << "Maximum nr of image uniforms supported by fragment shader: " << nrAttributes << std::endl;
    glGetIntegerv(GL_MAX_COMPUTE_IMAGE_UNIFORMS, &nrAttributes);
    std::cout << "Maximum nr of image uniforms supported by compute shader: " << nrAttributes << std::endl;

    // Compute Shader Configuration
    int work_grp_cnt[3], work_grp_inv;

    glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &work_grp_cnt[0]);
    glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 1, &work_grp_cnt[1]);
    glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 2, &work_grp_cnt[2]);

    printf("max global (total) work group size x:%i y:%i z:%i\n",
      work_grp_cnt[0], work_grp_cnt[1], work_grp_cnt[2]);

    glGetIntegerv(GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS, &work_grp_inv);
    printf("max local work group invocations %i\n", work_grp_inv);

    // Mouse input mode
    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSwapInterval(0); // No fps constraint
    
    return window;
}
