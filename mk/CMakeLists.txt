# create target mk and add sources in mk folder
add_library(mk STATIC "mk_mod.F" "mk_pc_mod.F")
# set dependencies
target_link_libraries(mk driver_setup multigrid)
# not compiled by default
set_target_properties(mk PROPERTIES EXCLUDE_FROM_ALL TRUE)
# add to all target
add_dependencies(framework mk)
# add folder for module files
target_include_directories(mk PUBLIC "${CMAKE_CURRENT_BINARY_DIR}")

