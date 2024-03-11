# Relocate source dir
set(SRC_PATH "${CMAKE_CURRENT_SOURCE_DIR}/")
configure_file("${CONFIGURE_SRC_PATH_TEMPLATE}" "${CMAKE_CURRENT_SOURCE_DIR}/cmake_source_dir.h" @ONLY)

# 1. Write binary format
add_executable(write_binary "write_binary.cpp")

# 3. Read binary format
add_executable(read_binary "read_binary.cpp")

