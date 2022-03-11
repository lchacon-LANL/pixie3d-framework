# create target multigrid and add sources in multigrid folder
add_library(
  multigrid STATIC
  "coupledMG.F"
)
# set dependencies
target_link_libraries(multigrid matrix_solver variables)

# do standard actions on target
init_target(multigrid)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

