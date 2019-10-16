# Locate the glew library
#
# This module defines the following variables:
#
# GLEW_LIBRARY the name of the library;
# GLEW_INCLUDE_DIR where to find glew include files.
# GLEW_FOUND true if both the GLEW_LIBRARY and GLEW_INCLUDE_DIR have been found.
#
# To help locate the library and include file, you can define a
# variable called GLFW3_ROOT which points to the root of the glfw library
# installation.
#
# default search dirs
# 
# Cmake file from: https://github.com/daw42/glslcookbook

set( _glew_HEADER_SEARCH_DIRS
"/usr/include"
"/usr/local/include"
"${CMAKE_SOURCE_DIR}/includes"
"C:/Program Files (x86)/glew/include" )
set( _glew_LIB_SEARCH_DIRS
"/usr/lib"
"/usr/lib64"
"/usr/local/lib"
"${CMAKE_SOURCE_DIR}/lib"
"C:/Program Files (x86)/glew/lib-msvc110" )

# Check environment for root search directory
set( _glew_ENV_ROOT $ENV{GLEW_ROOT} )
if( NOT GLEW_ROOT AND _glew_ENV_ROOT )
	set(GLEW_ROOT ${_glew_ENV_ROOT} )
endif()

# Put user specified location at beginning of search
if( GLEW_ROOT )
	list( INSERT _glew_HEADER_SEARCH_DIRS 0 "${GLEW_ROOT}/include" )
	list( INSERT _glew_LIB_SEARCH_DIRS 0 "${GLEW_ROOT}/lib" )
endif()

# Search for the header
FIND_PATH(GLEW_INCLUDE_DIR "GL/glew.h"
PATHS ${_glew_HEADER_SEARCH_DIRS} )

# Search for the library
FIND_LIBRARY(GLEW_LIBRARY NAMES GLEW
PATHS ${_glew_LIB_SEARCH_DIRS} )
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GLEW DEFAULT_MSG
GLEW_LIBRARY GLEW_INCLUDE_DIR)
