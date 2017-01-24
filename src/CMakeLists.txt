include_directories ("${PROJECT_SOURCE_DIR}/include")

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

add_library (BeBeReader SHARED DataReader.cc)
target_include_directories (BeBeReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library (BeBePulse SHARED Pulse.cc)
target_include_directories (BeBePulse PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
