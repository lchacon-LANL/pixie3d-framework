# create target eigen and add sources in eigen folder
add_library(
  eigen STATIC
  "eigen.F"
)
# set dependencies
target_link_libraries(eigen dependencies)

# do standard actions on target
init_target(eigen)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

