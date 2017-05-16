cmake_minimum_required(VERSION 3.4.3 FATAL_ERROR)
project (bebe)

include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${ROOTSYS}/include")
add_definitions (-g -O3 -Wall)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib)

add_subdirectory (src)
add_subdirectory (bin)
