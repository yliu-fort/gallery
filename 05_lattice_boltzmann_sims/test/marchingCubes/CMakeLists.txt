# MARCHING CUBE
#include_directories("utility")

# Relocate source dir
set(SRC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/")
configure_file("${CONFIGURE_SRC_PATH_TEMPLATE}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake_source_dir.h" @ONLY)

           file(GLOB SOURCE
               "*.h"
               "*.cpp"
               "shader.*"
               #${UTILITY}
               )

           add_executable(marching_cubes ${SOURCE})
target_link_libraries(marching_cubes ${LIBS})

