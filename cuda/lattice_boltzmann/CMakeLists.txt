cmake_minimum_required (VERSION 2.8)
cmake_policy(VERSION 2.8)

project (rtCFD)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules")

link_directories(${CMAKE_SOURCE_DIR}/lib)

#list(APPEND CMAKE_CXX_FLAGS "-O3 -std=c++11")

include_directories("includes")
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

set(LIBS glfw3 GL X11 pthread Xrandr Xi dl STB_IMAGE assimp GLEW)

# then create a project file per tutorial
###
add_executable(rtCFD
    "main.cpp"
    "shader.cpp"
    "camera.cpp"
)

target_link_libraries(rtCFD ${LIBS})
###
add_executable(rtCFD_offscreen
    "offscreen.cpp"
    "shader.cpp"
    "camera.cpp"
)

target_link_libraries(rtCFD_offscreen ${LIBS})
