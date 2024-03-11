# IMGUI
add_executable(imgui_demo "imgui_demo.cpp")
#target_include_directories(imgui_demo PUBLIC
#  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../external/imgui-1.74/example>
#  $<INSTALL_INTERFACE:external/imgui-1.74>  # <prefix>/include/mylib
#)
target_link_libraries(imgui_demo ${LIBS} IMGUI)
