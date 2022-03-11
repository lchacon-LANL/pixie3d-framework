# create target polynomial and add sources in polynomial folder
add_library(
  polynomial STATIC
  "DPLINT.f" "DPOLCF.f" "DPOLVL.f"
)
# set dependencies
target_link_libraries(polynomial dependencies)
# set Fortran standard to allow modern compilers to build this
target_compile_options(polynomial PRIVATE -std=legacy)

# do standard actions on target
init_target(polynomial)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

