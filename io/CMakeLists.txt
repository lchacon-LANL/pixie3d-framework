# create target io and add sources in io folder
add_library(
  io STATIC
  "graphics_io_mod.F" "io_mod.F" "xdraw_io_mod.F"
)
# set dependencies
target_link_libraries(io dependencies)

# do standard actions on target
init_target(io)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

