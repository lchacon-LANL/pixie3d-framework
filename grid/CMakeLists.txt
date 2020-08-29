# create target grid and add sources in grid folder
add_library(
  grid STATIC
  "grid_anal_map_mod.F"
  "grid_debug_mod.F"
  "grid_metric_mod.F"
  "grid_mpi_mod.F"
  "grid_petsc_mod.F"
  "grid_bc_mod.F"
  "grid_def_mod.F"
  "grid_mg_mod.F"
  "grid_operations_mod.F"
  "grid_spline_mod.F"
  "grid_create_mod.F"
  "grid_diffops_mod.F"
  "grid_mod.F"
  "grid_pack_mod.F"
)
# set dependencies
target_link_libraries(grid 1d_int OMP)

# do standard actions on target
init_target(grid)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

