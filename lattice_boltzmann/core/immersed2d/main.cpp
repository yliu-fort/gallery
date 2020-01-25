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
#include "gui_interface.h"
#include "slice.h"
#include "boundingbox.h"
#include "image_io.h"

#include "lbm2dgpu.h"

// settings
static int SCR_WIDTH  = 800;
static int SCR_HEIGHT = 600;

// Grid size, must be a integer multiplier of 8 for kernel tiling
//static glm::uvec2 scale(glm::uvec2(4,1)*4u);
//static uint nx = (16)*scale.x;
//static uint ny = (16)*scale.y;
//static uint nz = (1);

// camera
static Camera camera = Camera(glm::vec3(1.0f, 1.0f, 5.0f), float(SCR_WIDTH)/SCR_HEIGHT);

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

// Display texture
static int selectTexture = 0;
static int renderType = 0;
static float slidebar = 0.5f;

// Dynamic control with gui interface
float tau_gui = 0.5f;

// simulation related parameter
static bool initFluid = true;
static uint timestep = 0;

// Initial customized shaders...
// File path must be modified before compiling the script.
static Shader fluid2DInitShader;
static Shader fluid2DComputeShader;

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
    fluid2DInitShader.reload_shader_program_from_files(FP("fluid2d.init.glsl"));
    fluid2DComputeShader.reload_shader_program_from_files(FP("fluid2d.solve.glsl"));

    // For auto-reloading
    FileSystemMonitor::Init(SRC_PATH);
    GuiInterface::Init(window);

    // Gen colormap
    Colormap::Viridis();
    Colormap::Bind(15);

    // Isosurface
    Slice::Init();

    // configure g-buffer framebuffer
    // ------------------------------
    Lattice2DOpenGLInterface::Init(512*4,512, false);

    while( !glfwWindowShouldClose( window ) )
    {
        // per-frame time logic
        // --------------------
        timestep = Lattice2DOpenGLInterface::GetTimestep();
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
            Lattice2DOpenGLInterface::Reset();
        }
        fluid2DComputeShader.use();
        fluid2DComputeShader.setFloat("tau_gui",tau_gui);

        // Shader configuration
        if (FileSystemMonitor::Update())
        {
            fluid2DInitShader.reload_shader_program_from_files(FP("fluid2d.init.glsl"));
            fluid2DComputeShader.reload_shader_program_from_files(FP("fluid2d.solve.glsl"));
            Slice::ReloadShader();
        }

        // Initialize simulation
        if(Lattice2DOpenGLInterface::Isinitializing()) {fluid2DInitShader.use();}
        else { fluid2DComputeShader.use(); }

        Lattice2DOpenGLInterface::Compute();

        // 2. rendering pass
        if(!Lattice2DOpenGLInterface::IsInInternalCycle()){
            glViewport(0,0,SCR_WIDTH, SCR_HEIGHT);
            glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            Lattice2DOpenGLInterface::BindTexRead();
            Slice::Draw(2, slidebar,camera, glm::ivec3(Lattice2DOpenGLInterface::NX(),Lattice2DOpenGLInterface::NY(),Lattice2DOpenGLInterface::NZ()));
            //Draw bounding box
            Boundingbox::Draw(camera.GetFrustumMatrix()*glm::translate(glm::scale(glm::mat4(1.0), 0.5f*
                                                                                  glm::vec3(Lattice2DOpenGLInterface::NX()/float(Lattice2DOpenGLInterface::NY()), 1.0f,
                                                                                            Lattice2DOpenGLInterface::NZ()/float(Lattice2DOpenGLInterface::NY()))), glm::vec3(1)));
            // Draw gui (on top of other layers so put this to the last order)
            GuiInterface::Draw();
        }

        // Autosave
        //if(Lattice2DOpenGLInterface::GetTimestep()%4000==100 ) Lattice2DOpenGLInterface::Autosave();

        // Output vtk result
        //if(
        //        (Lattice2DOpenGLInterface::GetTimestep()>1000 &&
        //         Lattice2DOpenGLInterface::GetTimestep()%10==0) ||
        //        Lattice2DOpenGLInterface::GetTimestep()==5 // for test
        //        )
        //{
        //    Lattice2DOpenGLInterface::DumpAll();
        //}

        // Output image
        //if(Lattice2DOpenGLInterface::GetTimestep()==50) ImageIO::Save(SCR_WIDTH, SCR_HEIGHT, Lattice2DOpenGLInterface::GetTimestep(),1);

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    Lattice2DOpenGLInterface::Finalize();
    glfwTerminate( );

    return 0;
}

void countAndDisplayFps(GLFWwindow* window)
{
    float currentFrame = float(glfwGetTime());
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;

    frameCount++;
    if(float(glfwGetTime()) - lastFpsCountFrame > 1.0f)
    {
        /*std::cout << "Current fps: "
                  << frameCount/(glfwGetTime() - lastFpsCountFrame)
                  << " runtime:"
                  << glfwGetTime()
                  << std::endl;*/ // deprecated

        char title [256];
        title [255] = '\0';

        snprintf ( title, 255,
                   "FAST+ARB DEMO - FPS: %4.2f | runtime: %.0fs | nstep: %d",
                   double(frameCount/(glfwGetTime() - double(lastFpsCountFrame))), glfwGetTime(),timestep );
        glfwSetWindowTitle(window, title);

        frameCount = 0;
        lastFpsCountFrame = float(glfwGetTime());
    }
    if(deltaTime > 600.0f) {
        std::cout << "No response for 600 sec... exit program." << std::endl;
        glfwTerminate();
        EXIT_FAILURE;
    }
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
    if(!window) return;
    glViewport(0, 0, w, h);
    camera.updateAspect(float(w) / float(h));
}

static bool mouse_button_right = false;
// glfw: whenever the mouse moves, this callback is called
// -------------------------------------------------------
void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if(!window) return;
    if (firstMouse)
    {
        lastX = float(xpos);
        lastY = float(ypos);
        firstMouse = false;
    }

    float xoffset = float(xpos) - lastX;
    float yoffset = lastY - float(ypos); // reversed since y-coordinates go from bottom to top

    lastX = float(xpos);
    lastY = float(ypos);

    if(mouse_button_right)
        camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    if(!window) return;
    camera.ProcessMouseScroll(float(0*xoffset));
    camera.ProcessMouseScroll(float(yoffset));
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if(!window) return;
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
    {mouse_button_right = true;return;}
    mouse_button_right = false;
    mods = 0;
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
    GLFWwindow* window = glfwCreateWindow( w, h, "Demo", nullptr, nullptr);
    if( window == nullptr ){
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
