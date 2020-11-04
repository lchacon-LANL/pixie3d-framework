# file to generate targets for dependencies

# define function to generate targets
# ######################################################################################################################
function(generate package base_dir)
    # check if contrib has been unpacked
    if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${base_dir})
        message(STATUS "Unpackig contrib archive")
        # unpack zip file
        execute_process(
            COMMAND
                ${CMAKE_COMMAND} -E tar xzf ${CMAKE_CURRENT_SOURCE_DIR}/../common_contrib.tgz
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
        )
    endif()

    # define file extensions
    set(file_extensions
        "f"
        "F"
        "f90"
        "F90"
    )
    # inlcude file extensions
    set(include_extensions "fpp" "h")

    foreach(file_ext IN LISTS file_extensions)
        list(APPEND source_globs "${base_dir}/*.${file_ext}")
    endforeach()

    foreach(file_ext IN LISTS include_extensions)
        list(APPEND include_globs "${base_dir}/*.${file_ext}")
    endforeach()

    # find source files
    file(GLOB_RECURSE sources ${source_globs})
    # find include files
    file(GLOB_RECURSE includes ${include_globs})
    # add target
    add_library(${package} STATIC ${sources} ${includes})
    # set Fortran standard to allow modern compilers to build these
    target_compile_options(${package} PRIVATE -std=legacy)
    # do not build automatically, only if linked to other targets or explicitly specified
    set_target_properties(${package} PROPERTIES EXCLUDE_FROM_ALL TRUE)
    # add to the all target
    add_dependencies(framework ${package})
    # add path to include path for modules
    target_include_directories(${package} PUBLIC "${CMAKE_CURRENT_BINARY_DIR}")
endfunction()

# create contrib libs
# ######################################################################################################################

# fixed point iterations
# ######################################################################################################################
generate(fpa "fpa") # find source files and create library target
target_compile_definitions(fpa PRIVATE "NDEBUG") # set special compiler flags for target
target_link_libraries(fpa) # set dependencies

# slatec
# ######################################################################################################################
generate(slatec "slatec") # find source files and create library target
target_link_libraries(slatec) # set dependencies

# rng
# ######################################################################################################################
generate(rng "rng") # find source files and create library target
target_link_libraries(rng) # set dependencies

# ptridiag
# ######################################################################################################################
if(USE_PETSc)
    generate(ptridiag "ptridiag") # find source files and create library
    # target
    target_link_libraries(ptridiag MPI::MPI_Fortran) # set dependencies
endif()

# arpack ! FIXME not tested
# ######################################################################################################################
generate(arpack "arpack") # find source files and create library target
target_link_libraries(arpack) # set dependencies

# btridiag ! FIXME not tested
# ######################################################################################################################
generate(btridiag "btridiag") # find source files and create library
# target
target_link_libraries(btridiag) # set dependencies

# lsode ! FIXME not tested
# ######################################################################################################################
generate(lsode "lsode") # find source files and create library target
target_link_libraries(lsode) # set dependencies

# parareal ! FIXME not tested
# ######################################################################################################################
# generate(parareal "parareal") # find source files and create library target target_link_libraries(parareal) # set
# dependencies

# sdc ! FIXME test folder does not compile
# ######################################################################################################################
# generate(sdc "sdc") # find source files and create library target target_link_libraries(sdc) # set dependencies

# options dependent on contrib targets
# ######################################################################################################################

# arpack
if(USE_ARPACK)
    message(STATUS "Enabling ARPACK")

    # add to interface
    target_link_libraries(dependencies INTERFACE arpack)
    # set compiler flags
    target_compile_definitions(dependencies INTERFACE "arpack")
endif()

