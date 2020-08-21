# create target driver_setup and add sources in driver_setup folder
add_library(
  driver_setup STATIC
  "proj_int_mod.F"
  "ts_if_mod.F"
  "ts_mod.F"
  "ts_setup_mod.F"
  "ts_base_mod.F"
  "ts_io_mod.F"
  "ts_nk_setup_mod.F"
)
# set dependencies
target_link_libraries(driver_setup variables newtonGMRES eigen)

# do standard actions on target
init_target(driver_setup)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

