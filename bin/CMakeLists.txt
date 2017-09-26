include_directories ("${PROJECT_SOURCE_DIR}/include")

find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist)
include_directories(${ROOT_INCLUDE_DIRS})

# -- bin2root
add_executable (bin2root bin2root.cc)
target_link_libraries (bin2root bebe_reader ${ROOT_LIBRARIES})

# -- quickana
add_executable (quickana quickana.cc)
target_link_libraries (quickana bebe_pulse ${ROOT_LIBRARIES})

# -- player
add_executable (player player.cc)
target_link_libraries (player bebe_pulse bebe_tree ${ROOT_LIBRARIES})

# -- monito
add_executable (monito monito.cc)
target_link_libraries (monito bebe_pulse bebe_tree ${ROOT_LIBRARIES})
