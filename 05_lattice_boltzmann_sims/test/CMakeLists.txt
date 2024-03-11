# Relocate source dir
set(SRC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/")
configure_file("${CONFIGURE_SRC_PATH_TEMPLATE}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake_source_dir.h" @ONLY)

# OPENGL FUNCTIONALITY
add_executable(hello_triangle "triangle.cpp")
target_link_libraries(hello_triangle ${LIBS})

# UTILITIES
add_executable(macro_expansion "macro_expansion.cpp")
add_executable(inotify "inotify.cpp")

# DEMOS
add_subdirectory("marchingCubes")
add_subdirectory("computeShader")
add_subdirectory("imGUI")
add_subdirectory("IO")
