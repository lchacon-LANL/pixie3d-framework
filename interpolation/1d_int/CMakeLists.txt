# create target 1d_int and add sources in 1d_int folder
add_library(
  1d_int STATIC
  "1dIntDriver_mod.f"
  "DCHFDV.f"
  "DCHFEV.f"
  "DPCHDF.f"
  "DPCHFD.f"
  "DPCHFE.f"
  "DPCHSP.f"
)
# set dependencies
target_link_libraries(1d_int math slatec)

# do standard actions on target
init_target(1d_int)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

