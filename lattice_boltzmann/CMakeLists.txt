cmake_minimum_required (VERSION 2.8)
cmake_policy(VERSION 2.8)

project (fast_arb)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules")

link_directories(${CMAKE_SOURCE_DIR}/lib)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

### GLOBAL CONFIGURATIONS

# project libraries
include_directories("utility")

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

add_library(STB_IMAGE "includes/stb/stb_image.cpp")

set(LIBS glfw GL X11 pthread Xrandr Xi dl STB_IMAGE assimp GLEW)

set(CONFIGURE_SRC_PATH_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/cmake_source_dir.h.in")
message(STATUS "CONFIGURE_SRC_PATH_TEMPLATE at ${CMAKE_CURRENT_SOURCE_DIR}/cmake_source_dir.h.in")

file(GLOB UTILITY
    "utility/*.h"
    "utility/*.cpp"
    "utility/glsl/*"
    )

### CASE CONFIGURATION
add_subdirectory("utility")
add_subdirectory("core/2d")
add_subdirectory("core/3d")

### tests
add_subdirectory("test")
