# create target matrix_solver and add sources in matrix_solver folder
add_library(
  matrix_solver STATIC
  "mlSolver.F" "mlsolver_mod.f"
)
# set dependencies
target_link_libraries(matrix_solver grid)

# do standard actions on target
init_target(matrix_solver)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

