# include subfolders
add_subdirectory(2d)
add_subdirectory(3d)
add_subdirectory(hyb3d)

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

