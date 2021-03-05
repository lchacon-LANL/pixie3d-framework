# create target math and add sources in math folder
add_library(
  math STATIC
  "math_mod.F"
)
# set dependencies
target_link_libraries(math io)
# set legacy standard because of deleted features FIXME
target_compile_options(math PRIVATE -std=legacy)

# do standard actions on target
init_target(math)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

