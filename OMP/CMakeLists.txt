# create target OMP and add sources in OMP folder
add_library(
  OMP STATIC
  "OMP_mod.F"
)
# set dependencies
target_link_libraries(OMP dependencies)

# do standard actions on target
init_target(OMP)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

