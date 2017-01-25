include_directories ("${PROJECT_SOURCE_DIR}/include")

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})

# -- bin2root
add_executable (bin2root bin2root.cc)
target_link_libraries (bin2root LINK_PUBLIC bebe_reader Core RIO Tree)

# -- quickana
add_executable (quickana quickana.cc)
target_link_libraries (quickana LINK_PUBLIC bebe_pulse Core RIO Tree Hist)
