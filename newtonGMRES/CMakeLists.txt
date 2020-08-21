# create target newtonGMRES and add sources in newtonGMRES folder
add_library(
  newtonGMRES STATIC
  "nk_mod.F"
)
# set dependencies
target_link_libraries(newtonGMRES grid fpa)

# do standard actions on target
init_target(newtonGMRES)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

