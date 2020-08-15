# Created by the script cgal_create_CMakeLists
# This is the CMake script for compiling a set of CGAL applications.

project( executable )


cmake_minimum_required(VERSION 3.1)#(VERSION 2.8.11)
SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
SET(BUILD_SHARED_LIBS OFF)
#SET(CMAKE_EXE_LINKER_FLAGS "-static")
SET(CMAKE_BUILD_TYPE "RELEASE")
SET(CGAL_DIR ${CMAKE_INSTALLED_PREFIX}/lib/CGAL)
# Forward the Ceres QUIET / REQUIRED options.
#Set(Ceres_DIR /home/jason/CUDA/ceres-solver-1.14.0/ceres-bin)
#find_package(Ceres REQUIRED)
#include_directories(${CERES_INCLUDE_DIRS})
#link_directories(${CERES_LIBRARY_DIRS})


# CGAL and its components
add_library(CGAL STATIC IMPORTED COMPONENTS)
find_package( CGAL QUIET COMPONENTS  )

include( ${CGAL_USE_FILE} )
include( CGAL_CreateSingleSourceCGALProgram )

if ( NOT CGAL_FOUND )

  message(STATUS "This project requires the CGAL library, and will not be compiled.")
  return()  

endif()


# Boost and its components
find_package( Boost REQUIRED )

if ( NOT Boost_FOUND )

  message(STATUS "This project requires the Boost library, and will not be compiled.")

  return()  

endif()

# include for local directory

# include for local package
include_directories( BEFORE ../include )


# Creating entries for target: executable
set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
# ############################

#add_executable(pipline_poly "pipline_poly.cpp")
#target_link_libraries(pipline_poly -static)
# ############################

#create_single_source_cgal_program( "multi_polygon.cpp")
add_executable( multi_polygon multi_polygon.cpp polygon.cpp polygon.h mesh.cpp mesh.h )
