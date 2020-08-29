# create target 3dplot and add sources in 3dplot folder
add_library(
  3dplot STATIC
  "3dplot.F" "ADIOS_graphics_mod.F" "graphics_mod.F" "acc_test_mod.F"
  "HDF5_graphics_mod.F" "xdraw_graphics_mod.F"
)
# set dependencies
target_link_libraries(3dplot driver_setup)

# do standard actions on target
init_target(3dplot)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

