include_directories ("${PROJECT_SOURCE_DIR}/include")

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

add_library (bebeReader SHARED DataReader.cc)
target_include_directories (BebeReader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library (beBePulse SHARED Pulse.cc)
target_include_directories (bebePulse PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
