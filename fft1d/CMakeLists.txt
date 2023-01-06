# create target fft and add sources in fft folder
add_library(
  fft STATIC
  "fft.f" "test.f"
)
# set dependencies
target_link_libraries(fft 1d_int)

# do standard actions on target
init_target(fft)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

