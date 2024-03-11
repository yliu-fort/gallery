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
#include "slice.h"
#include "boundingbox.h"
#include "image_io.h"

#include "lbm3dgpu.h"

// settings
static int SCR_WIDTH  = 800;
static int SCR_HEIGHT = 600;

// Grid size, must be a integer multiplier of 8 for kernel tiling
static glm::ivec3 scale(glm::ivec3(8,2,2)*4);
static int nx = (8)*scale.x;
static int ny = (8)*scale.y;
static int nz = (8)*scale.z;

// camera
static Camera camera = Camera(glm::vec3(1.0f, 1.0f, 5.0f), (float)SCR_WIDTH/SCR_HEIGHT);

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
static float slidebar = 0.05f;

// simulation related parameter
static bool initFluid = true;
static int timestep = 0;

// Initial customized shaders...
// File path must be modified before compiling the script.
static Shader fluid3DInitShader;
static Shader fluid3DComputeShader;

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
    fluid3DInitShader.reload_shader_program_from_files(FP("fluid3d.init.glsl"));
    fluid3DComputeShader.reload_shader_program_from_files(FP("fluid3d.solve.glsl"));

    // For auto-reloading
    FileSystemMonitor::Init(SRC_PATH);

    // Gen colormap
    Colormap::Plasma();
    Colormap::Bind(15);

    // Isosurface
    IsoSurface::Init();
    Slice::Init();

    // configure g-buffer framebuffer
    // ------------------------------
    Lattice3DOpenGLInterface::Create(nx, ny, nz, true);

    while( !glfwWindowShouldClose( window ) )
    {
        // per-frame time logic
        // --------------------
        timestep = Lattice3DOpenGLInterface::GetTimestep();
        countAndDisplayFps(window);

        // input
        // -----
        glfwGetFramebufferSize(window, &SCR_WIDTH, &SCR_HEIGHT);
        processInput(window);

        // Reset parameters
        if(initFluid)
        {
            initFluid = false;
            slidebar = 0.5f;
            Lattice3DOpenGLInterface::Reset();
        }

        // Shader configuration
        if (FileSystemMonitor::Update())
        {
            fluid3DInitShader.reload_shader_program_from_files(FP("fluid3d.init.glsl"));
            fluid3DComputeShader.reload_shader_program_from_files(FP("fluid3d.solve.glsl"));
            IsoSurface::ReloadShader();
            Slice::ReloadShader();
        }

        // Initialize simulation
        if(Lattice3DOpenGLInterface::Isinitializing()) {fluid3DInitShader.use();}
        else { fluid3DComputeShader.use(); }

        Lattice3DOpenGLInterface::Compute();

        // 2. rendering pass
        if(!Lattice3DOpenGLInterface::IsInInternalCycle())
        {
            glViewport(0,0,SCR_WIDTH, SCR_HEIGHT);
            glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            Lattice3DOpenGLInterface::BindTexRead();
            IsoSurface::Draw(2, 0.5,camera,glm::ivec3(nx,ny,nz));
            Slice::Draw(2, slidebar,camera,glm::ivec3(nx,ny,nz));
            Boundingbox::Draw(camera.GetFrustumMatrix()*glm::translate(glm::scale(glm::mat4(1.0),0.5f*glm::vec3(nx/(float)ny,1.0f, nz/(float)ny)), glm::vec3(1)));
        }
        // For extreme large simulations, need an offscreen version without texture rendering...(isosurface is really heavy when > 16M grid points)
        //Draw bounding box
        //glViewport(0,0,SCR_WIDTH, SCR_HEIGHT);
        //glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        //Boundingbox::Draw(camera.GetFrustumMatrix()*glm::translate(glm::scale(glm::mat4(1.0),glm::vec3(nx/(float)ny,1.0f, nz/(float)ny)), glm::vec3(1)));

        // Autosave
        if(Lattice3DOpenGLInterface::GetTimestep()%2000==1000 ) Lattice3DOpenGLInterface::Autosave();

        // Output vtk result
        //if(
        //        //(Lattice3DOpenGLInterface::GetTimestep()>1000 &&
        //        // Lattice3DOpenGLInterface::GetTimestep()%200==0)
        //        )
        //{
        //    Lattice3DOpenGLInterface::DumpAll();
        //}

        // Auto exit
        if(Lattice3DOpenGLInterface::GetTimestep() > 360000){glfwTerminate();EXIT_SUCCESS;}

        //if(Lattice3DOpenGLInterface::GetTimestep()==50) ImageIO::Save(SCR_WIDTH, SCR_HEIGHT, Lattice3DOpenGLInterface::GetTimestep(),1);

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
                   "FAST+ARB DEMO - FPS: %4.2f | runtime: %.0fs | nstep: %d | Resolution %dx%dx%d: %d",
                   frameCount/(glfwGetTime() - lastFpsCountFrame), glfwGetTime(),timestep, nx,ny,nz,nx*ny*nz );
        glfwSetWindowTitle(window, title);

        frameCount = 0;
        lastFpsCountFrame = glfwGetTime();
    }
    if(deltaTime > 600.0f) {
        std::cout << "No response for 600 sec... exit program." << std::endl;
        glfwTerminate();
        EXIT_FAILURE;
    }
}

// renderCube() renders a 1x1 3D cube in NDC.
// -------------------------------------------------
static unsigned int cubeVAO = 0;
static unsigned int cubeVBO = 0;
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
static unsigned int quadVAO = 0;
static unsigned int quadVBO = 0;
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
        slidebar += 0.1f*deltaTime;
    if (glfwGetKey(window, GLFW_KEY_F) == GLFW_PRESS)
        slidebar -= 0.1f*deltaTime;

    if ( (glfwGetKey(window, GLFW_KEY_I) == GLFW_PRESS) && (key_i_old_state == GLFW_RELEASE) ) {

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
    glGetIntegerv(GL_MAX_SHADER_STORAGE_BLOCK_SIZE, &nrAttributes);
    std::cout << "GL_MAX_SHADER_STORAGE_BLOCK_SIZE is " << nrAttributes << " bytes." << std::endl;
    glGetIntegerv(GL_MAX_SHADER_STORAGE_BUFFER_BINDINGS, &nrAttributes);
    std::cout << "Maximum nr of shader storage buffer binding points is " << nrAttributes << " ." << std::endl;
    
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
