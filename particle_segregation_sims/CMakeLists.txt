cmake_minimum_required(VERSION 3.5)

project(syclode)

option(NVIDIA_BUILD "Build demo for nvidia device" ON)
option(GUI_BUILD "Build interactive demo" ON)

# C++11 enabled
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")


# Find ComputeCpp library
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/modules)
find_package(ComputeCpp REQUIRED)
add_definitions(-D_ALLOW_COMPILER_AND_STL_VERSION_MISMATCH) # workarounds for msvc2019

################################################
# set denpendencies
# set includes
INCLUDE_DIRECTORIES( "${CMAKE_CURRENT_SOURCE_DIR}/includes")
INCLUDE_DIRECTORIES( SYSTEM "${ComputeCpp_INCLUDE_DIRS}")
#INCLUDE_DIRECTORIES( SYSTEM "${OpenCL_INCLUDE_DIRS}")

set(SYCL_SOURCES
    "src/syclode45.cpp"
    )


################################################
# add library for sycl module
add_library(sycl_module STATIC
    ${SYCL_SOURCES}
    )
target_include_directories(sycl_module
    PRIVATE SYSTEM "${ComputeCpp_INCLUDE_DIRS}"
    )

add_sycl_to_target(
    TARGET sycl_module
    SOURCES ${SYCL_SOURCES}
    )


################################################
# executable cannot be launched from outside of qtcreator
# add QT/bin into PATH variable in windows
# Tell CMake to create the helloworld executable

link_libraries("-lstdc++fs")

add_executable(sycl_ode45 "src/main.cpp" "src/filesystem.cpp")
target_link_libraries(sycl_ode45 ${LIBS} sycl_module)


add_executable(sycl_ks2d "src/main.cpp" "src/kinetic_turbulence.cpp" "src/filesystem.cpp")
add_sycl_to_target( TARGET sycl_ks2d SOURCES "src/kinetic_turbulence.cpp" )



add_executable(sycl_ks3d "src/pair_diffusion3d_main.cpp" "src/pair_diffusion3d.cpp" "src/filesystem.cpp")
add_sycl_to_target( TARGET sycl_ks3d SOURCES "src/pair_diffusion3d.cpp" )


# tests
add_executable(simpleode "tests/rk45.cpp" )
add_sycl_to_target( TARGET simpleode SOURCES "tests/rk45.cpp" )


add_executable(filesystem "tests/fs.cpp" "src/filesystem.cpp" )
file(GLOB TEST_INPUT "resources/*.txt")
foreach(INPUT_FILE IN LISTS TEST_INPUT)
    get_filename_component(INPUT_FILE "${INPUT_FILE}" NAME)
    message(STATUS "COPY INPUT FILE=${INPUT_FILE}")
    configure_file("resources/${INPUT_FILE}" "resources/${INPUT_FILE}")
endforeach()

add_executable( random_engine "tests/rd.cpp" )

add_executable(test_ks2d_host "tests/single_traj.cpp" "src/filesystem.cpp")
add_sycl_to_target( TARGET test_ks2d_host SOURCES "tests/single_traj.cpp" )

add_executable(test_ks3d_host "tests/single_traj3d.cpp" "src/filesystem.cpp")
add_sycl_to_target( TARGET test_ks3d_host SOURCES "tests/single_traj3d.cpp" )

################################################
if(GUI_BUILD)
    ### gui part
    ## if dll not found, adding the QTDIR C:\Qt\5.4\msvc2017_64\ to system environment variable
    # Turn on automatic invocation of the MOC, UIC & RCC
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTORCC ON)

    # There may be a way to tell up front if Qt5 is going to be found, but I haven't found
    # a foolproof way to do it yet, so settle for the default error message for now.
    #if(NOT CMAKE_PREFIX_PATH AND NOT Qt5Widgets_DIR)
    #    message(WARNING "CMAKE_PREFIX_PATH is not defined, so find_package may not work. Set the CMAKE_PREFIX_PATH "
    #            "environment variable to the install prefix of Qt 5, either on the command line as "
    #            "-DCMAKE_PREFIX_PATH=\"path/to/Qt5/lib/cmake\" or with set(CMAKE_PREFIX_PATH path/to/Qt5/lib/cmake)")
    #endif(NOT CMAKE_PREFIX_PATH AND NOT Qt5Widgets_DIR)

    # Find the QtWidgets library
    find_package(Qt5 REQUIRED COMPONENTS Widgets)

    ################################################
    # set denpendencies
    if(WIN32)
        set(LIBS opengl32) # glfw3 if use binary libs
    elseif(UNIX AND NOT APPLE)
        set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
        find_package(OpenGL REQUIRED)
        add_definitions(${OPENGL_DEFINITIONS})
        find_package(X11 REQUIRED)
        # note that the order is important for setting the libs
        # use pkg-config --libs $(pkg-config --print-requires --print-requires-private glfw3) in a terminal to confirm
        set(LIBS X11 Xrandr Xinerama Xi Xxf86vm Xcursor GL dl pthread)
        set (CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -ldl")
    elseif(APPLE)
        INCLUDE_DIRECTORIES(/System/Library/Frameworks)
        FIND_LIBRARY(COCOA_LIBRARY Cocoa)
        FIND_LIBRARY(OpenGL_LIBRARY OpenGL)
        FIND_LIBRARY(IOKit_LIBRARY IOKit)
        FIND_LIBRARY(CoreVideo_LIBRARY CoreVideo)
        MARK_AS_ADVANCED(COCOA_LIBRARY OpenGL_LIBRARY)
        SET(APPLE_LIBS ${COCOA_LIBRARY} ${IOKit_LIBRARY} ${OpenGL_LIBRARY} ${CoreVideo_LIBRARY})
        #SET(APPLE_LIBS ${APPLE_LIBS} ${GLFW3_LIBRARY})
        set(LIBS ${LIBS} ${APPLE_LIBS})
    else()
        set(LIBS )
    endif(WIN32)


    ################################################
    # executable cannot be launched from outside of qtcreator
    # add QT/bin into PATH variable in windows
    # Tell CMake to create the helloworld executable

    include_directories(${OpenCL_INCLUDE_DIRS})
    link_directories(${OpenCL_LIBRARY})


    # set sources 2d
    set(SOURCES2D
        "src/guimain.cpp"
        "src/visual_ks2d.cpp"
        "src/openglwindow.cpp"
        "src/filesystem.cpp"
        )

    add_executable(visual_ks2d ${SOURCES2D})
    if(WIN32)
        set_property(TARGET visual_ks2d PROPERTY WIN32_EXECUTABLE true)
    endif()
    # Add the Qt5 Widgets for linking
    target_link_libraries(visual_ks2d  Qt5::Widgets ${LIBS} )
    add_sycl_to_target( TARGET visual_ks2d SOURCES "src/visual_ks2d.cpp" )


    # set sources 3d
    set(SOURCES3D
        "src/guimain3d.cpp"
        "src/visual_ks3d.cpp"
        "src/openglwindow.cpp"
        "src/filesystem.cpp"
        )

    add_executable(visual_ks3d ${SOURCES3D})
    if(WIN32)
        set_property(TARGET visual_ks3d PROPERTY WIN32_EXECUTABLE true)
    endif()
    # Add the Qt5 Widgets for linking
    target_link_libraries(visual_ks3d  Qt5::Widgets ${LIBS} )
    add_sycl_to_target( TARGET visual_ks3d SOURCES "src/visual_ks3d.cpp" )
endif()
