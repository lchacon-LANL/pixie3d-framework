cmake_minimum_required(VERSION 3.13)

# flag if we run as a subproject
set(subproject OFF)
# check if project name is already set
if(DEFINED PROJECT_NAME)
  # switch flag
  set(subproject ON)
endif()

# define project
project(framework)

# enable Fortran compiler
enable_language(Fortran)

message(STATUS "Configuring framework")

# Build type
# ##############################################################################
# Set a default build type if none was specified
set(default_build_type "Release")

# set a default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "${default_build_type}"
      CACHE STRING "Choose the type of build." FORCE
  )
endif()

# Set the possible values of build type for cmake-gui
set_property(
  CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo"
)

# print build type
message(STATUS "Building configuration '${CMAKE_BUILD_TYPE}'.")

# modules
# ##############################################################################
# add modules to find packages, e.g. petsc
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/"
)

# include necessary modules
include(CMakeDependentOption)
include(CheckIPOSupported)

# options
# ##############################################################################
# define interface target, this allows us to set properties for all targets
add_library(
  dependencies INTERFACE
)

# add link time optimazation check for link time optimazation
check_ipo_supported(RESULT lto_supported OUTPUT "IPO / LTO not supported")
# option to use LTO, only available if LTO is available
cmake_dependent_option(
  use_LTO "Use Link time optimazation" ON
  "lto_supported;NOT CMAKE_BUILD_TYPE MATCHES Debug" OFF
)

# function  if option is available and if not, turn off
function(framework_option option_name description default)
  # options are only available, if this is not includes as a subproject,
  # otherwise the app has to set these
  if(NOT subproject)
    # create option
    option(${option_name} ${description} ${default})

    # check if option is defined by main project
  elseif(NOT DEFINED ${option_name})
    # default option to off
    set(${option_name}
        ${default}
        PARENT_SCOPE
    )
  endif()
endfunction()

# function if option is available and if not, turn off
function(framework_dependent_option option_name description default condition disabled_default)
  # options are only available, if this is not includes as a subproject,
  # otherwise the app has to set these
  if(NOT subproject)
    # create option
    cmake_dependent_option(${option_name} ${description} ${default} ${condition} ${disabled_default})

    # check if option is defined by main project
  elseif(NOT DEFINED ${option_name})
    # default option to off
    set(${option_name}
        ${disabled_default}
        PARENT_SCOPE
    )
  endif()
endfunction()

framework_option(USE_OPENMP "Use OpenMP for shared memory parallelization" OFF)
framework_option(USE_PETSc "Enable parallelization using PETSc" OFF)
framework_option(USE_SPRNG "Use SPRNG5 random number generator" OFF)
framework_option(USE_USER_PROVIDED_DT "No idea" OFF)
framework_option(USE_NO_BDF2_XOLD "No idea" OFF)
framework_option(USE_coarse_MG "No idea" OFF)
framework_option(USE_OLD_BC "No idea" OFF)
framework_option(USE_PER_BC_SYNC "No idea" OFF)
framework_option(USE_MK "No idea" OFF)
framework_option(USE_PIC_CUDA "Use CUDA for PIC libs" OFF)
framework_option(USE_HDF5 "Use HDF5 libraries" OFF)
framework_option(USE_ARPACK "Use ARPACK libraries" OFF)
framework_option(USE_SAMRAI "Use SAMRAI" OFF)
framework_option(
    USE_Werror "Turns all warnings into errors, good to fix errors" ON
)

framework_dependent_option(USE_ADIOS2 "Enable ADIOS2 io library" ON "USE_PETSc" OFF)
# adios is only available if PETSc is enabled
framework_dependent_option(USE_ADIOS "Enable ADIOS io library" OFF "USE_PETSc;NOT USE_ADIOS2" OFF)

# apply options to dependencies target
# ##############################################################################
# check if we use LTO
if(use_LTO)
  # echo on screen
  message(STATUS "Link time optimization enabled")
  # set flags
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

if(USE_USER_PROVIDED_DT)
  target_compile_definitions(dependencies INTERFACE "USER_PROVIDED_DT")
endif()

if(USE_NO_BDF2_XOLD)
  target_compile_definitions(dependencies INTERFACE "NO_BDF2_XOLD")
endif()

if(USE_coarse_MG)
  target_compile_definitions(dependencies INTERFACE "coarse_MG")
endif()

if(USE_OLD_BC)
  target_compile_definitions(dependencies INTERFACE "OLD_BC")
endif()

if(USE_PER_BC_SYNC)
  target_compile_definitions(dependencies INTERFACE "PER_BC_SYNC")
endif()

# compiler specific settings
# ##############################################################################
# GNU Fortran
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  # enable preprocessor
  target_compile_options(dependencies INTERFACE "-cpp")
  # set gfortran preprocessor constant
  target_compile_definitions(dependencies INTERFACE "gfortran")

  # FIXME turn on warnings, skip unused dummy argument for now
  # target_compile_options(dependencies INTERFACE -fimplicit-none -Wall -Wno-
  # unused-dummy-argument -Wextra -Wno-compare-reals -Wimplicit-interface
  # -Wimplicit- procedure )

  # check if we want warnings to be errors
  if(USE_Werror)
    # add error flag
    target_compile_options(dependencies INTERFACE -Werror)
  endif()

  # debug flags
  if(CMAKE_BUILD_TYPE MATCHES "Debug")
    # add bracktrace and floating point exception flags
    target_compile_options(
      dependencies
      INTERFACE -fbacktrace -fbounds-check -ffpe-trap=invalid -ffpe-trap=zero
                -ffpe-trap=overflow -ffpe-trap=underflow
                # -ffpe-trap=denormal
    )
  endif()
endif()

# Intel compiler
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  # enable preprocessor
  target_compile_options(dependencies INTERFACE -fpp)
  # set intel compile time constant
  target_compile_definitions(dependencies INTERFACE "ifort")

  # FIXME add warnings and implicit checking

  # set debug flags
  if(CMAKE_BUILD_TYPE MATCHES "Debug")
    target_compile_options(dependencies INTERFACE -check all)
  endif()
endif()

# add debug compiler flags
if(CMAKE_BUILD_TYPE MATCHES "Debug")
  target_compile_definitions(dependencies INTERFACE "DEBUG")
else()
  target_compile_definitions(dependencies INTERFACE "NDEBUG")
endif()

# find required packages
# ##############################################################################
# LAPACK
# ##############################################################################
# define input variable for vendor
if(NOT DEFINED BLA_VENDOR)
  # default is empty
  set(BLA_VENDOR
      ""
      CACHE STRING "BLAS and LAPACK vendor" FORCE
  )
endif()

# define options
set_property(
  CACHE BLA_VENDOR
  PROPERTY STRINGS
           ""
           "Goto"
           "OpenBLAS"
           "FLAME"
           "ATLAS PhiPACK"
           "CXML"
           "DXML"
           "SunPerf"
           "SCSL"
           "SGIMATH"
           "IBMESSL"
           "Intel10_32"
           "Intel10_64lp"
           "Intel10_64lp_seq"
           "Intel10_64ilp"
           "Intel10_64ilp"
           "Intel"
           "ACML"
           "ACML_MP"
           "ACML_GPU"
           "Apple"
           "NAS"
           "Generic"
)

# set(BLA_STATIC ON)
find_package(LAPACK REQUIRED)
# BLAS
find_package(BLAS REQUIRED)

# add LAPACK libraries to interface, this will propagate flags, includes etc
# through the project
target_link_libraries(
  dependencies INTERFACE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES}
)

# PETSc
# ##############################################################################
if(USE_PETSc)
  message(STATUS "Enabling PETSc")

  # MPI
  find_package(MPI REQUIRED)
  # petsc
  find_package(PETSC 3.0.0 EXACT REQUIRED)

  # add to interface
  target_link_libraries(
    dependencies INTERFACE PETSC::PETSC_Fortran MPI::MPI_Fortran
  )

  # set compiler flags for the framework
  target_compile_definitions(dependencies INTERFACE NVAR=8 petsc)

endif()

# add OpenMP support
# ##############################################################################
if(USE_OPENMP)
  message(STATUS "Enabling OpenMP")

  find_package(OpenMP REQUIRED)
  # add libraries to interface
  target_link_libraries(dependencies INTERFACE OpenMP::OpenMP_Fortran)
endif()

# other packages
# ##############################################################################
# adios2
if(USE_ADIOS2)
    message(STATUS "Enabling ADIOS 2")

    find_package(ADIOS2 REQUIRED)
    
    # check for parallel
    if (USE_PETSc)
        # add ADIOS MPI to interface
        target_link_libraries(dependencies INTERFACE adios2::fortran_mpi)
        # add compiler flags
        target_compile_definitions(
            dependencies INTERFACE "ADIOS2" "ADIOS2_VER=${ADIOS2_VERSION}"
                                 "ADIOS2_FILE=\"record-${ADIOS2_VERSION}.bp\""
            )

    else()
        # add ADIOS serial to interface, no flags, since framework does not 
        # support serial version at the moment
        target_link_libraries(dependencies INTERFACE adios2::fortran)

    endif()

endif()

# adios
if(USE_ADIOS)
  message(STATUS "Enabling ADIOS")

  find_package(ADIOS 1.10.0 EXACT REQUIRED)
  # add ADIOS to interface
  target_link_libraries(dependencies INTERFACE ADIOS::ADIOS)
  # add compiler flags
  target_compile_definitions(
    dependencies INTERFACE "adios" "ADIOS_VER=${ADIOS_VERSION}"
                           "ADIOS_FILE=\"record-${ADIOS_VERSION}.bp\""
  )

  # set buffer size
  set(ADIOS_BUFFER_SIZE
      ""
      CACHE STRING "Size of ADIOS buffer in MB"
  )

  # set buffer size
  if(NOT ${ADIOS_BUFFER_SIZE} STREQUAL "")
    target_compile_definitions(
      dependencies INTERFACE "ADIOS_BUFFER_MB=${ADIOS_BUFFER_SIZE}"
    )
  endif()

  # set method
  set(ADIOS_METHOD
      ""
      CACHE STRING "ADIOS method"
  )

  # add option for ADIOS method
  if(NOT ${ADIOS_METHOD} STREQUAL "")
    target_compile_definitions(
      dependencies INTERFACE "ADIOS_METHOD=\"record-${ADIOS_METHOD}.bp\""
    )
  endif()
endif()

# SPRNG5
if(USE_SPRNG)
  message(STATUS "Enabling SPRNG5")

  find_package(SPRNG5 5.0.0 EXACT REQUIRED)
  # add to interface
  target_link_libraries(dependencies INTERFACE SPRNG5::SPRNG5)
  # set compiler flags
  target_compile_definitions(
    dependencies INTERFACE "_SPRNG" "HAVE_CONFIG_H" "POINER_SIZE=8"
                           "DLONG64=long" "INTEGER_STAR_8"
  )
endif()

# HDF5
if(USE_HDF5)
  message(STATUS "Enabling HDF5")

  message(FATAL "HDF5 not yet supported!")

  # FIXME
  find_package(HDF REQUIRED)
  # add to interface
  target_link_libraries(dependencies INTERFACE hdf::hdf)
  # set compiler flags
  target_compile_definitions(dependencies INTERFACE "hdf5" "H5_USE_16_API")
endif()

# SAMRAI
if(USE_SAMRAI)
  message(STATUS "Enabling SAMRAI")

  message(FATAL "SAMRAI not yet supported!")

  # FIXME find_package(SAMRAI REQUIRED) add to interface
  target_link_libraries(dependencies INTERFACE)
  # set compiler flags
  target_compile_definitions(dependencies INTERFACE "samrai")
endif()

# ##############################################################################
# a target to compile all target in this file, allows testing of the framework
add_custom_target(framework)

# ##############################################################################
# contributions
add_subdirectory(contrib)

# Framework targets
# ##############################################################################

# list for CMake files
file(GLOB_RECURSE cmake_files "cmake/**/*.cmake")
list(APPEND cmake_files "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
# list of source files
list(APPEND framework_sources "")

# function for actions on all targets
function(init_target target_name)
  # not compiled by default
  set_target_properties(${target_name} PROPERTIES EXCLUDE_FROM_ALL TRUE)
  # add to all target
  add_dependencies(framework ${target_name})
  # add folder for module files
  target_include_directories(
    ${target_name} PUBLIC "${CMAKE_CURRENT_BINARY_DIR}"
  )
  # add CMake file to list
  list(APPEND cmake_files "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
  # add source files to list
  get_target_property(target_sources ${target_name} SOURCES)

  # turn into absolute path
  foreach(target_source IN LISTS target_sources)
    list(APPEND framework_sources
         "${CMAKE_CURRENT_SOURCE_DIR}/${target_source}"
    )
  endforeach()

  # set these values in outer scope (kind of return value)
  set(cmake_files
      ${cmake_files}
      PARENT_SCOPE
  )
  set(framework_sources
      ${framework_sources}
      PARENT_SCOPE
  )

endfunction()

# targets defined in subfolders
add_subdirectory(3dplot)
add_subdirectory(driver)
add_subdirectory(driver_setup)
add_subdirectory(eigen)
add_subdirectory(fft)
add_subdirectory(grid)
add_subdirectory(interpolation)
add_subdirectory(io)
add_subdirectory(math)
add_subdirectory(matrixLightSolver)
add_subdirectory(MG)
add_subdirectory(newtonGMRES)
add_subdirectory(OMP)
add_subdirectory(orbit)
add_subdirectory(var_setup)

# below are packages that do not yet work
# ##############################################################################
# FIXME missing module spline_field

# add_subdirectory(mk)

# add_subdirectory(pic)

# format target (only if not build as subproject)
# ##############################################################################
if(NOT subproject)

  # find fprettify
  find_program(FPRETTIFY "fprettify")
  # find cmake-format
  find_program(CMAKE_FORMAT "cmake-format")

  # check for clang-format
  if(FPRETTIFY AND CMAKE_FORMAT)

    # define script to remove white spaces
    set(format_whitespaces
        "${CMAKE_CURRENT_SOURCE_DIR}/scripts/utility/remove_whitespaces"
    )

    # add format target
    add_custom_target(
      format
      # format CMake files
      COMMAND ${CMAKE_FORMAT} -c ${CMAKE_CURRENT_SOURCE_DIR}/.cmake-format.py -i
              ${cmake_files}
      COMMAND ${format_whitespaces} ${cmake_files}
      COMMAND ${format_whitespaces} ${framework_sources}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      USES_TERMINAL VERBATIM
      COMMENT "Formatting sources"
    )

  else()
    # print warning
    message(WARNING "Format target not available")
    # check for fprettyfy
    if(NOT FPRETTIFY)
      message(
        WARNING
          "No fprettify executable found, formatting target won't be available! Install with "
          "'pip install --user fprettify' and set variable in cmake"
      )
    endif()
    # check for cmake-format
    if(NOT CMAKE_FORMAT)
      message(
        WARNING
          "No cmake-format executable found, formatting target won't be available! Install with "
          "'pip install --user cmake_format' and set variable in cmake"
      )
    endif()

  endif()
endif()

# print final message
message(STATUS "Finished configuring framework")

