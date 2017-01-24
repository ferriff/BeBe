cmake_minimum_required(VERSION 3.4.3 FATAL_ERROR)
project (BeBe)

include_directories ("${PROJECT_SOURCE_DIR}/include")
include_directories ("${ROOTSYS}/include")
add_definitions (-g -O3 -Wall)

add_subdirectory (src)
add_subdirectory (bin)
