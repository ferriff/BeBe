include_directories ("${PROJECT_SOURCE_DIR}/include")

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

# -- test_reading
add_executable (test_reading test_reading.cc)
target_link_libraries (test_reading LINK_PUBLIC bebe_reader Core RIO Tree)

# -- quickana
add_executable (quickana quickana.cc)
target_link_libraries (quickana LINK_PUBLIC bebe_pulse Core RIO Tree Hist)
