include_directories ("${PROJECT_SOURCE_DIR}/include")

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

add_library (bebe_reader SHARED data_reader.cc)
target_include_directories (bebe_reader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library (bebe_pulse SHARED pulse.cc)
target_include_directories (bebe_pulse PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
