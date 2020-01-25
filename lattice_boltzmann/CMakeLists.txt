cmake_minimum_required (VERSION 2.8)
cmake_policy(VERSION 2.8)

project (fast_arb)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules")

link_directories(${CMAKE_SOURCE_DIR}/lib)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")

### SET SOURCE PATH (todo: deprecate)
set(CONFIGURE_SRC_PATH_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/cmake_source_dir.h.in")
message(STATUS "CONFIGURE_SRC_PATH_TEMPLATE at ${CMAKE_CURRENT_SOURCE_DIR}/cmake_source_dir.h.in")

### GLOBAL CONFIGURATIONS
# project libraries
#include_directories("external")
#First make sure you have CMake, Git, and GCC by typing as root (sudo) apt-get install g++ cmake git
#and then get the required packages: Using root (sudo) and type apt-get install libsoil-dev libglm-dev libassimp-dev libglew-dev libglfw3-dev libxinerama-dev libxcursor-dev libxi-dev
# find the required packages
find_package(GLM REQUIRED)
message(STATUS "GLM included at ${GLM_INCLUDE_DIR}")

find_package(GLFW3 REQUIRED)
message(STATUS "Found GLFW3 in ${GLFW3_INCLUDE_DIR}")

find_package(ASSIMP REQUIRED)
message(STATUS "Found ASSIMP in ${ASSIMP_INCLUDE_DIR}")

# find_package(SOIL REQUIRED)
# message(STATUS "Found SOIL in ${SOIL_INCLUDE_DIR}")

find_package(GLEW REQUIRED)
message(STATUS "Found GLEW in ${GLEW_INCLUDE_DIR}")

## INTEGRATE ESSENTIAL GRAPHIC LIBRARIES
set(LIBS glfw GL X11 Xrandr Xi dl GLEW assimp)

# IMAGE_IO LIB
include_directories("external/stb")
add_library(STB_IMAGE "external/stb/stb_image.cpp")
set(LIBS ${LIBS} STB_IMAGE)

# GUI LIB
# Copy imgui_impl_glfw & imgui_impl_opengl3 out to root folder first
include_directories("external/imgui-1.74")
include_directories("external/imgui-1.74/examples")
file(GLOB IMGUI_SOURCES "external/imgui-1.74/*.cpp"
    "external/imgui-1.74/examples/imgui_impl_glfw.cpp"
    "external/imgui-1.74/examples/imgui_impl_opengl3.cpp")
add_library(IMGUI ${IMGUI_SOURCES})
set(LIBS ${LIBS} IMGUI)

# MY UTILITY LIB
include_directories("utility")
file(GLOB UTILITY_SOURCES "utility/*" "utility/glsl/*")
add_library(UTILITY ${UTILITY_SOURCES})
target_link_libraries(UTILITY IMGUI)
set(LIBS ${LIBS} UTILITY)
add_subdirectory("utility") # Provide absolute path for shaders

### CASE CONFIGURATION
add_subdirectory("core/2d")
#add_subdirectory("core/immersed2d")
add_subdirectory("core/3d")
add_subdirectory("core/2d_tex")
#add_subdirectory("core/3d_surf")
### TESTS&DEMOS
add_subdirectory("test")
