# create target pic_2d and add sources in pic_2d folder
add_library(pic_2d STATIC "Fcc_mod.F" "pic_base_mod.F")
# compile definitions
target_compile_definitions(pic_2d PRIVATE _Npg=4 _prec=8 PIC_CM_DIAG)
# set dependencies
target_link_libraries(pic_2d rng driver_setup multigrid)

# do standard actions on target
init_target(pic_2d)

# because CMake is stupid, we have to do this
set(cmake_files
    ${cmake_files}
    PARENT_SCOPE
)
set(framework_sources
    ${framework_sources}
    PARENT_SCOPE
)

