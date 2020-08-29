# create target pic_3d and add sources in pic_3d folder
add_library(pic_3d STATIC "Fcc_mod.F" "pic3d_mod.F" "pic3d_base_mod.F")
# compile definitions
target_compile_definitions(pic_3d PRIVATE _Npg=4 _prec=8 PIC_CM_DIAG)
# set dependencies
target_link_libraries(pic_3d rng grid)

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

