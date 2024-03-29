#include <iostream>
#include <cmath>

//GLEW
#define GLEW_STATIC
#include <GL/glew.h>

//GLFW
#include <GLFW/glfw3.h>

#include "camera.h"
#include "filesystemmonitor.h"
//#include "gui_interface.h"
#include "image_io.h"
#include "defines.h"
#include "cartesian_field.h"

// settings
static int SCR_WIDTH  = 800;
static int SCR_HEIGHT = 600;

// camera
static Camera camera = Camera(glm::vec3(0.0f, 0.0f, 5.0f), float(SCR_WIDTH)/SCR_HEIGHT);

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
bool countAndDisplayFps(GLFWwindow* window);
void processInput(GLFWwindow *window);

//#define VISUAL
#ifndef VISUAL
    #define BATCH_RUN
    #ifdef BATCH_RUN
        //#define BATCH_TEST
    #endif
#endif
// Shortcut
static bool pause = true;

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

#ifndef BATCH_TEST

    // Read shader
    Vec4::reload_subroutines();
    algo::reload_subroutines();
    renderer::reload_subroutines();

#endif

    // For automatic file reloading
    FileSystemMonitor::Init(SRC_PATH);

#ifdef BATCH_RUN
    // Generate a param pool
    uint numcase = 2*3*5;
    Param pool[numcase];

    // Edit configure pool
    float initDistRatio[3] = {0.5,1,2};
    uint rseeds[3] = {123, 98765432, 38496};
    float time_reso[3] = {0.005,0.002,0.0005};
    float numeric_tol[3] = {1e-2,1e-3,1e-4};
    float total_t[3] = {30,10,5};
    glm::uvec3 initDim[3] = {glm::uvec3(2048,16,1),glm::uvec3(16,2048,1),glm::uvec3(256,256,1)};

    float tau_inv[5] = {100,10,5,2,1};
    float g[2] = {-9.81,0};

    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 3; j++)
        {
            for(int k = 0; k < 5; k++)
            {
                int index = 3*5*i+5*j+k;

                pool[index].dim = glm::uvec3(1024,16,1);
                pool[index].sim_type = algo::INERTIAL_PARTICLE;

                pool[index].g = g[i];
                pool[index].initDist *= initDistRatio[j];
                pool[index].tau_inv = tau_inv[k];


                //pool[index].initDist *= initDistRatio[i];
                //pool[index].tol = numeric_tol[k];
                //pool[index].randomSeed = rseeds[j];
                //pool[index].initDist *= initDistRatio[j];
                //pool[index].dt = time_reso[i];

                //pool[index].t_end = total_t[i];
                //pool[index].t_end = total_t[k] < total_t[i] ? total_t[k] : total_t[i];

            }
        }
    }

    for(uint ipool = 0; ipool < numcase; ipool++)
    {
        pool[ipool].printInfo();
    }

#ifdef BATCH_TEST
    goto exit;
#endif

    for(uint ipool = 0; ipool < numcase; ipool++)
    {
        // configure g-buffer framebuffer
        Param param = pool[ipool];
        param.genericDir();

#else
    Param param;

    // Overwrite default parameters here
    param.dim = glm::uvec3(1024,16,1);
    param.t = 0.0f;param.t_end = 30.0f;
    param.dt = 0.002f;
    param.tau_inv = 10.0;
    param.g = -9.81;
    param.sim_type = algo::INERTIAL_PARTICLE;

#endif
    // Prepare output dir
    param.printInfoToDir();

    std::shared_ptr<Cartesian2d> mesh(new Cartesian2d(param.dim.x,param.dim.y,param.dim.z,param.grid.x,param.grid.y,param.grid.z));

    // Initialize data
    std::vector<glm::vec4> data;

    for(uint k = 0; k < param.dim.z; k++)
    {
        for(uint j = 0; j < param.dim.y; j++)
        {
            for(uint i = 0;i < param.dim.x; i++)
            {
                // initial distance ~ 1/K_nk
                data.push_back(glm::vec4(param.initDist*(i - 0.5f*param.dim.x),
                                         param.initDist*(j - 0.5f*param.dim.y),
                                         param.initDist*(k - 0.5f*param.dim.z),1));
            }
        }
    }

    Vec4 pos(data, mesh);
    Vec4 vel(mesh);

    // Configure kinetic simulation
    uint ks_ubo; glGenBuffers(1, &ks_ubo);
    algo::ode45_init(ks_ubo, param);

    // Dump initial data
    fclose(fopen ((param.getDir() + "/" + "ccd.mean.dat").c_str(), "w")); // Reset file
    pos.gather();
    async_io::dump("ccd", param.getDir(), pos.get_hptr(),param);

    fclose(fopen ((param.getDir() + "/" + "vel.mean.dat").c_str(), "w")); // Reset file
    vel.gather();
    async_io::dump("vel", param.getDir(), vel.get_hptr(),param);

    while( !glfwWindowShouldClose( window ) )
    {
        // per-frame time logic
        // --------------------
        countAndDisplayFps(window);

#ifdef VISUAL
        // input
        glfwGetFramebufferSize(window, &SCR_WIDTH, &SCR_HEIGHT);
        processInput(window);

        // Snapshot
        //if(param.needDump()) {ImageIO::Save(SCR_WIDTH, SCR_HEIGHT, param.getIter(), true);}
        //printf("Current time: %f\n",param.getCurrentTimestep());

        // Compute
        if(!pause)
        {
#endif

            algo::ode45(pos, vel, param);

            param.advCurrentTimestep();
            param++;

#ifndef VISUAL
            // autosave && exit
            if(param.needDump()) {
                pos.gather();
                async_io::dump("ccd", param.getDir(), pos.get_hptr(),param);
                vel.gather();
                async_io::dump("vel", param.getDir(), vel.get_hptr(),param);
            }
            if(param.timeOver()) {break;}
            //pos.report_minmax();
#endif

#ifdef VISUAL
        }

        // Draw points
        glViewport(0,0,SCR_WIDTH, SCR_HEIGHT);
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
        renderer::draw(pos,camera);
#endif

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    async_io::wait();

#ifdef BATCH_RUN
}
#ifdef BATCH_TEST
exit:
#endif
#endif
glfwTerminate( );

return 0;
}

bool countAndDisplayFps(GLFWwindow* window)
{
    float currentFrame = float(glfwGetTime());
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;

    frameCount++;
    if(float(glfwGetTime()) - lastFpsCountFrame > 1.0f)
    {
        std::cout << "Current fps: "
                  << frameCount/(glfwGetTime() - lastFpsCountFrame)
                  << " runtime:"
                  << glfwGetTime()
                  << std::endl; // deprecated

        frameCount = 0;
        lastFpsCountFrame = float(glfwGetTime());
        return true;
    }
    if(deltaTime > 60.0f) {
        std::cout << "No response for 60 sec... exit program." << std::endl;
        glfwTerminate();
        EXIT_FAILURE;
    }
    return false;
}
#ifdef VISUAL
static int key_space_old_state = GLFW_RELEASE;
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

    if ((glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS) && (key_space_old_state == GLFW_RELEASE))
    {
        pause = !pause;
    }
    key_space_old_state = glfwGetKey(window, GLFW_KEY_SPACE);
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
#endif

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
#ifndef VISUAL
    glfwWindowHint(GLFW_VISIBLE, GL_FALSE);
#endif

    // Open a window and create its OpenGL context
    GLFWwindow* window = glfwCreateWindow( w, h, "", nullptr, nullptr);
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

#ifdef VISUAL
    // framebuffer mode
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    // Mouse input mode
    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSwapInterval(1); // 60 fps constraint
#else
    glfwSwapInterval(0); // No fps constraint
#endif
    
    return window;
}
