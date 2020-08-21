# create target orbit and add sources in orbit folder
add_library(
  orbit STATIC
  "lyapn_mod.F" "orbit_mod.F" "par_int_mod.F" "poincare_io_mod.F"
)
# set dependencies
target_link_libraries(orbit grid driver_setup)

# do standard actions on target
init_target(orbit)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

