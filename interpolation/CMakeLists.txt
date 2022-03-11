# include sub folders
add_subdirectory(1d_int)
add_subdirectory(polynomial)

# add this file
list(APPEND cmake_files "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

