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

#include "shader.h"
#include "camera.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

// settings
static int SCR_WIDTH  = 1024;
static int SCR_HEIGHT = 1024;

static int nx = 8192;
static int ny = 8192;

static int coherent_nstep = (48*ny)/1024;

// camera
Camera camera(glm::vec3(0.0f, 0.0f, 5.0f), (float)SCR_WIDTH/(float)SCR_HEIGHT);
float lastX = (float)SCR_WIDTH / 2.0;
float lastY = (float)SCR_HEIGHT / 2.0;
bool firstMouse = true;

// timing
float deltaTime = 0.0f;
float lastFrame = 0.0f;

// fps recording
static float lastFpsCountFrame = 0;
static int frameCount = 0;

// Pre-declaration
GLFWwindow* initGL(int w, int h);
void framebuffer_size_callback(GLFWwindow* window, int w, int h);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow *window);
unsigned int loadTexture(const char *path, bool gammaCorrection);
void renderQuad();
void renderCube();

// Display texture
int selectTexture = 0;
int renderType = 1;
bool initFluid = true;

class Lattice2D
{
    // 4 fixed bind points for textures
    // 0 for uvrf;
    // 1 for f1_4;
    // 2 for f5_8;
    // 3 for occl;
public:
    unsigned int buffer[2];
    unsigned int uvrf[2];
    unsigned int f1_4[2];
    unsigned int f5_8[2];
    unsigned int occl[2];
    unsigned int rbo[2];
    int nx, ny, nbstep;


    Lattice2D(int w, int h):nx(w), ny(h), selector(0), nbstep(0)
    {
        glGenFramebuffers(2, buffer);
        glGenTextures(2, uvrf);
        glGenTextures(2, f1_4);
        glGenTextures(2, f5_8);
        glGenTextures(2, occl);
        glGenRenderbuffers(2, rbo);

        setupBuffer(0);
        setupBuffer(1);

    }

    void bind(){glBindFramebuffer(GL_FRAMEBUFFER, buffer[selector]);};
    void unbind(){glBindFramebuffer(GL_FRAMEBUFFER, 0);};
    void incr(){ selector = (selector+1)%2;nbstep++; };
    void debug() {std::cout << "Buffer " << selector << " is using." << std::endl; };
    void bindTextureRead()
    {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, uvrf[1-selector]);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, f1_4[1-selector]);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, f5_8[1-selector]);
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, occl[1-selector]);
    };
    void bindTextureWrite()
    {
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, uvrf[selector]);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, f1_4[selector]);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, f5_8[selector]);
        glActiveTexture(GL_TEXTURE3);
        glBindTexture(GL_TEXTURE_2D, occl[selector]);
    };
    int getOcclusionBindingPoint() const
    {
        return 3;
    };

private:
    int selector;

    void setupBuffer(int i)
    {
        GLenum internalFormat = GL_RGBA32F;
        GLenum format = GL_RGBA;
        GLenum datatype = GL_FLOAT;

        unsigned int gBuffer = buffer[i];

        glBindFramebuffer(GL_FRAMEBUFFER, gBuffer);
        unsigned int bufferTexture1 = uvrf[i];
        unsigned int bufferTexture2 = f1_4[i];
        unsigned int bufferTexture3 = f5_8[i];
        unsigned int bufferTexture4 = occl[i];
        // buffer 1 for u, v, rho, f0
        //glGenTextures(1, &bufferTexture1);
        glBindTexture(GL_TEXTURE_2D, bufferTexture1);
        glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, nx, ny, 0, format, datatype, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, bufferTexture1, 0);
        // buffer 2 for f1~f4
        //glGenTextures(1, &bufferTexture2);
        glBindTexture(GL_TEXTURE_2D, bufferTexture2);
        glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, nx, ny, 0, format, datatype, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, bufferTexture2, 0);
        // buffer 3 for f5~f8
        //glGenTextures(1, &bufferTexture3);
        glBindTexture(GL_TEXTURE_2D, bufferTexture3);
        glTexImage2D(GL_TEXTURE_2D, 0, internalFormat, nx, ny, 0, format, datatype, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT2, GL_TEXTURE_2D, bufferTexture3, 0);
        // buffer 4 for boundary
        //glGenTextures(1, &bufferTexture3);
        glBindTexture(GL_TEXTURE_2D, bufferTexture4);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, nx, ny, 0, GL_RGB, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT3, GL_TEXTURE_2D, bufferTexture4, 0);

        // tell OpenGL which color attachments we'll use (of this framebuffer) for rendering
        unsigned int attachments[4] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3 };
        glDrawBuffers(4, attachments);
        // create and attach depth buffer (renderbuffer)
        unsigned int rboDepth = rbo[i];
        //glGenRenderbuffers(1, &rboDepth);
        glBindRenderbuffer(GL_RENDERBUFFER, rboDepth);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, nx, ny);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rboDepth);
        // finally check if framebuffer is complete
        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
            std::cout << "Framebuffer not complete!" << std::endl;
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }
};

struct Lattice_Constant
{
    glm::mat3 wi;
    glm::mat3 cx;
    glm::mat3 cy;
    float csSqrInv;
    float csSqr;
    float Reg;
    float u_max;
    float tau;


    Lattice_Constant()
    {
        float w0 = 4.0/9.0; // zero weight
        float ws = 1.0/9.0; // adjacent weight
        float wd = 1.0/36.0; // diagonal weight
        wi = glm::mat3(w0,ws,ws,ws,ws,wd,wd,wd,wd);
        cx = glm::mat3(0.0,1.0,0.0,-1.0, 0.0,1.0,-1.0,-1.0, 1.0);
        cy = glm::mat3(0.0,0.0,1.0, 0.0,-1.0,1.0, 1.0,-1.0,-1.0);
        csSqrInv = 3.0;
        csSqr = 1.0/3.0;

        u_max = 0.25;
        Reg = 10.0;
        tau = 0.5 + u_max/(csSqr*Reg);

        // for acoustic simulations, tau -> 0.5 improves accuarcy (also introducing dispersive pattern)
        // for fluid simulations, tau > 0.6 improves to stability

        std::cout << "tau = " << tau << std::endl;
    }
};

class Framebuffer
{
public:

    Framebuffer()
    {
        // configure normal framebuffer
        //unsigned int framebuffer;
        glGenFramebuffers(1, &framebuffer);
        glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
        // render buffer
        //unsigned int rbo;
        glGenRenderbuffers(1, &rbo);
        glBindRenderbuffer(GL_RENDERBUFFER, rbo);
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, SCR_WIDTH, SCR_HEIGHT);
        glBindRenderbuffer(GL_RENDERBUFFER, 0);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rbo);

        // create a color attachment texture
        //unsigned int screenTexture;
        glGenTextures(1, &screenTexture);
        glBindTexture(GL_TEXTURE_2D, screenTexture);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, SCR_WIDTH, SCR_HEIGHT, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, screenTexture, 0);	// we only need a color buffer

        if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
            std::cout << "ERROR::FRAMEBUFFER:: Intermediate framebuffer is not complete!" << std::endl;
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        {
            float quadVertices[] = {
                // positions        // texture Coords
                -1.0f,  1.0f, 0.0f, 0.0f, 1.0f,
                -1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
                1.0f,  1.0f, 0.0f, 1.0f, 1.0f,
                1.0f, -1.0f, 0.0f, 1.0f, 0.0f,
            };
            // setup plane VAO
            glGenVertexArrays(1, &VAO);
            glGenBuffers(1, &VBO);
            glBindVertexArray(VAO);
            glBindBuffer(GL_ARRAY_BUFFER, VBO);
            glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), &quadVertices, GL_STATIC_DRAW);
            glEnableVertexAttribArray(0);
            glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
            glEnableVertexAttribArray(1);
            glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3 * sizeof(float)));
        }
    }
    void bind()
    {
        glBindFramebuffer(GL_FRAMEBUFFER, framebuffer);
    }
    void unbind()
    {
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }
private:
    unsigned int framebuffer;
    unsigned int rbo;
    unsigned int screenTexture;
    unsigned int VAO, VBO;
};

void PPMWriter(unsigned char *in, char *name, int dimx, int dimy)
{
    int i, j;
    FILE *fp = fopen(name, "wb");
    (void) fprintf(fp, "P6\n%d %d\n255\n", dimx, dimy);
    for(int j = 0; j < dimy; ++j)
    {
        for(i = 0; i < dimx; ++i)
        {
            static unsigned char color[3];
            color[0] = in[3*i + 3*(dimy - 1 - j)*dimx]; // Red
            color[1] = in[3*i + 3*(dimy - 1 - j)*dimx + 1]; // Green
            color[2] = in[3*i + 3*(dimy - 1 - j)*dimx + 2]; // Blue
            (void) fwrite(color, 1, 3, fp);
        }
    }
    (void) fclose(fp);
}

static unsigned int imgIndex = 1;
static unsigned char* image = NULL;
void saveImage(float Ttime)
{
    int w, h;
    w = SCR_WIDTH;
    h = SCR_HEIGHT;
    //unsigned char* image = (unsigned char*)malloc(sizeof(unsigned char) * 3* w * h);
    if(image == NULL)
        image = new unsigned char[4*w*h];
    glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, image);
    char buffer[255];
    sprintf(buffer, "/home/yliu/workspace/lbm/capture/result.%06d.png", imgIndex);

    //PPMWriter(image, buffer, w, h);
    stbi_flip_vertically_on_write(true);
    stbi_write_png(buffer, w, h, 3, image, w*3);

    //std::cout << "return from ppmwriter"<<std::endl;

    imgIndex++;
    //delete [] image;
    //free(image);
}


int main()
{
    // Initialize a window
    GLFWwindow* window = initGL(SCR_WIDTH, SCR_HEIGHT);

    // configure global opengl state
    // -----------------------------
    glEnable(GL_DEPTH_TEST);

    // Initial customized shaders...
    Shader shader_INIT("/home/yliu/workspace/lbm/vertex.mrt.shader",
                       "/home/yliu/workspace/lbm/fragment.init.shader");
    Shader shader_MRT("/home/yliu/workspace/lbm/vertex.mrt.shader",
                      "/home/yliu/workspace/lbm/fragment.mrt.shader");
    Shader shader_SCR("/home/yliu/workspace/lbm/vertex.shader",
                      "/home/yliu/workspace/lbm/fragment.shader");

    // Gen colormap
    // generate noise texture
    // ----------------------
    std::vector<glm::vec3> rainbow;
    unsigned int mapResolution = 64;
    for(unsigned int i = 0; i < mapResolution; i++)
    {
        /*plot short rainbow RGB*/
        float f = i/(mapResolution - 1.0f);
        float a=(1-f)/0.25;    //invert and group
        int X=floor(a);    //this is the integer part
        int Y=floor(255*(a-X)); //fractional part from 0 to 255
        float r,g,b;
        switch(X)
        {
        case 0: r=255;g=Y;b=0;break; // RED f = 1
        case 1: r=255-Y;g=255;b=0;break;
        case 2: r=0;g=255;b=Y;break;
        case 3: r=0;g=255-Y;b=255;break;
        case 4: r=0;g=0;b=255;break; // BLUE f = 0
        }

        glm::vec3 color(r/255.0f, g/255.0f, b/255.0f);
        rainbow.push_back(color);
    }

    unsigned int colorMap;
    glGenTextures(1, &colorMap);
    glBindTexture(GL_TEXTURE_1D, colorMap);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGB16, mapResolution, 0, GL_RGB, GL_FLOAT, &rainbow[0]);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    //glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // configure g-buffer framebuffer
    // ------------------------------
    Lattice2D demo(nx, ny);
    Lattice_Constant lc;
    //Framebuffer scr0;

    // Static configuration
    shader_INIT.use();
    shader_INIT.setInt("NX",nx);
    shader_INIT.setInt("NY",ny);
    shader_INIT.setMat3("wi",lc.wi);
    shader_INIT.setMat3("cx",lc.cx);
    shader_INIT.setMat3("cy",lc.cy);
    shader_INIT.setFloat("csSqrInv",lc.csSqrInv);
    shader_INIT.setFloat("csSqr",lc.csSqr);
    shader_INIT.setFloat("Reg",lc.Reg);
    shader_INIT.setFloat("u_max",lc.u_max);
    shader_INIT.setFloat("tau",lc.tau);
    shader_MRT.use();
    shader_MRT.setInt("uvrf_old", 0);
    shader_MRT.setInt("f1_4_old", 1);
    shader_MRT.setInt("f5_8_old", 2);
    shader_MRT.setInt("occl_old", demo.getOcclusionBindingPoint());
    shader_MRT.setInt("NX",nx);
    shader_MRT.setInt("NY",ny);
    shader_MRT.setMat3("wi",lc.wi);
    shader_MRT.setMat3("cx",lc.cx);
    shader_MRT.setMat3("cy",lc.cy);
    shader_MRT.setFloat("csSqrInv",lc.csSqrInv);
    shader_MRT.setFloat("csSqr",lc.csSqr);
    shader_MRT.setFloat("Reg",lc.Reg);
    shader_MRT.setFloat("u_max",lc.u_max);
    shader_MRT.setFloat("tau",lc.tau);
    shader_SCR.use();
    shader_SCR.setInt("occlusion", demo.getOcclusionBindingPoint());
    shader_SCR.setMat4("projectionMatrix", glm::scale(glm::mat4(1.0),
    glm::vec3(1.0/camera.Aspect,ny/(float)nx, 1.0)) );

    // Uniform Buffer (!not implemented)
    //unsigned int ubo;
    //glGenBuffers(1, &ubo);
    //glBindBuffer(GL_UNIFORM_BUFFER, ubo);
    //glBufferData(GL_UNIFORM_BUFFER, sizeof(Lattice_Constant) , &lc, GL_STATIC_DRAW);
    //glBindBuffer(GL_UNIFORM_BUFFER, 0);
    //unsigned int lc_index = glGetUniformBlockIndex(shader_INIT.ID, "Lattice_Constants");
    //glUniformBlockBinding(shader_INIT.ID, lc_index, 0);
    //glUniformBlockBinding(shader_MRT.ID, lc_index, 0);
    //glBindBufferBase(GL_UNIFORM_BUFFER, 0, ubo);

    while( !glfwWindowShouldClose( window ) )
    {
        // per-frame time logic
        // --------------------
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        frameCount++;
        if(glfwGetTime() - lastFpsCountFrame > 1.0f)
        {
            std::cout
                    <<"Current fps: "
                   << frameCount/(glfwGetTime() - lastFpsCountFrame)
                   <<" run step: "
                   << demo.nbstep/coherent_nstep
                   << std::endl;
            frameCount = 0;
            lastFpsCountFrame = glfwGetTime();
        }

        // input
        // -----
        glfwGetFramebufferSize(window, &SCR_WIDTH, &SCR_HEIGHT);
        processInput(window);

        // 渲染
        // 清空颜色缓冲
        glClearColor(0.2f, 0.2f, 0.2f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Initialize fluid
        if(initFluid)
        {
            demo.bind();
            glViewport(0,0,nx,ny);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            shader_INIT.use();
            renderQuad();
            demo.unbind();
            initFluid = false;
            demo.nbstep = 0;
            imgIndex = 1;
        }

        {
            // advance in time
            //demo.debug();
            demo.incr();

            // 1. computing pass (partition to 4 subregions)
            demo.bind();
            shader_MRT.use();
            shader_MRT.setInt("nstep",demo.nbstep);
            demo.bindTextureRead();
            glViewport(0,0,nx,ny);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            renderQuad();

            demo.unbind();
        }

        // 2. rendering pass
        //scr0.bind();
        glViewport(0,0,SCR_WIDTH, SCR_HEIGHT);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        shader_SCR.use();
        //shader_SCR.setMat4("projectionMatrix",camera.GetPerspectiveMatrix()*camera.GetViewMatrix()*glm::scale(glm::mat4(1.0), glm::vec3(1.0,ny/(float)nx, 1.0)));
        shader_SCR.setInt("renderTarget", selectTexture);
        shader_SCR.setInt("renderType", renderType);
        shader_SCR.setInt("colorMap",5);
        glActiveTexture(GL_TEXTURE5);
        glBindTexture(GL_TEXTURE_1D, colorMap);
        demo.bindTextureRead();
        renderQuad();

        // Save output
        if(demo.nbstep%coherent_nstep == 0) saveImage(0);

        //scr0.unbind();

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();


    }

    if(image)
        delete [] image;
    glfwTerminate( );

    return 0;
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

}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int w, int h)
{
    glViewport(0, 0, w, h);
    camera.updateAspect((float)w / (float)h);
}

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

    camera.ProcessMouseMovement(xoffset, yoffset);
}

// glfw: whenever the mouse scroll wheel scrolls, this callback is called
// ----------------------------------------------------------------------
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
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

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    //glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

    // Open a window and create its OpenGL context
    GLFWwindow* window = glfwCreateWindow( w, h, "", NULL, NULL);
    if( window == NULL ){
        fprintf( stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
        getchar();
        glfwTerminate();
        EXIT_FAILURE;
    }
    glfwMakeContextCurrent(window);

    // Initialize GLEW
    if (glewInit() != GLEW_OK) {
        fprintf(stderr, "Failed to initialize GLEW\n");
        getchar();
        glfwTerminate();
        EXIT_FAILURE;
    }

    // Initial viewport
    glViewport(0, 0, w, h);
    //glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // Query infomation
    int nrAttributes;
    glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &nrAttributes);
    std::cout << "Maximum nr of vertex attributes supported: " << nrAttributes << std::endl;
    std::cout << "Hardware: " <<glGetString(GL_RENDERER) << std::endl;
    glGetIntegerv(GL_MAX_DRAW_BUFFERS, &nrAttributes);
    std::cout << "Maximum nr of color attachments supported: " << nrAttributes << std::endl;

    // Mouse input mode
    //glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    //glfwSetCursorPosCallback(window, mouse_callback);
    //glfwSetScrollCallback(window, scroll_callback);
    glfwSwapInterval(0); // No fps constraint

    return window;
}
