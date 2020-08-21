# create target driver and add sources in driver folder
add_library(
  driver STATIC
  "driver.F" "evaluateNonlinearResidual.F" "initializeCalculation.F"
  "PC_ifces.F"
)
# set dependencies
target_link_libraries(driver driver_setup multigrid)

# do standard actions on target
init_target(driver)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

