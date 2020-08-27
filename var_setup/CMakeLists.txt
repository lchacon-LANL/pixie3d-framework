# create target variables and add sources in var_setup folder
add_library(
    variables STATIC
    "var_io_mod.F" "var_setup_mod.F" "ADIOS_io_mod.F" "ADIOS2_io_mod.F90"
    "var_bc_mod.F" "var_mod.F"
)
# set dependencies
target_link_libraries(variables grid)

# do standard actions on target
init_target(variables)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

