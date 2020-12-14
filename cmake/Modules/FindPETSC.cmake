# module to find a PETSC package

# find PETSC architecture variable
if(NOT PETSC_ARCH)
  set(PETSC_ARCH "$ENV{PETSC_ARCH}")
endif()

# find the PETSC header file
find_path(
  PETSC_DIR
  NAMES include/petsc.h
  HINTS ENV PETSC_DIR
  PATH_SUFFIXES ${PETSC_ARCH}
  DOC "PETSc root folder"
)

if(NOT PETSC_DIR MATCHES "PETSC_DIR-NOTFOUND")

  # get include path
  find_path(
    PETSC_INCLUDE_DIRS
    NAMES petsc.h
    HINTS ${PETSC_DIR}/include
  )

  # get version, set header file
  if(PETSC_INCLUDE_DIRS)
    set(_PETSC_VERSION_H "${PETSC_INCLUDE_DIRS}/petscversion.h")

    # function to extract version
    function(_VERSION_EXTRACT _VERSION_H _VER_COMPONENT _VER_OUTPUT)
      set(CMAKE_MATCH_1 "0")
      set(_expr "^[ \\t]*#define[ \\t]+${_VER_COMPONENT}[ \\t]+([0-9]+)$")
      file(STRINGS "${_VERSION_H}" _ver REGEX "${_expr}")
      string(REGEX MATCH "${_expr}" ver "${_ver}")
      set(${_VER_OUTPUT}
          "${CMAKE_MATCH_1}"
          PARENT_SCOPE
      )
    endfunction()

    # extract version
    _version_extract(
      ${_PETSC_VERSION_H} "PETSC_VERSION_MAJOR" PETSC_VERSION_MAJOR
    )
    _version_extract(
      ${_PETSC_VERSION_H} "PETSC_VERSION_MINOR" PETSC_VERSION_MINOR
    )
    _version_extract(
      ${_PETSC_VERSION_H} "PETSC_VERSION_SUBMINOR" PETSC_VERSION_SUBMINOR
    )
    _version_extract(
      ${_PETSC_VERSION_H} "PETSC_VERSION_PATCH" PETSC_VERSION_PATCH
    )

    # combine
    set(PETSC_VERSION
        "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}.${PETSC_VERSION_PATCH}"
    )

    # find the libraries
    find_library(
      PETSC_LIBRARIES
      NAMES "petsc"
      HINTS ${PETSC_DIR}/lib
    )

    # Determine whether the PETSc layout is old-style (through 2.3.3) or
    # new-style (>= 3.0.0)
    if(EXISTS "${PETSC_DIR}/${PETSC_ARCH}/lib/petsc/conf/petscvariables"
    )# > 3.5
      set(petsc_conf_rules "${PETSC_DIR}/lib/petsc/conf/rules")
      set(petsc_conf_variables "${PETSC_DIR}/lib/petsc/conf/variables")

    elseif(EXISTS "${PETSC_DIR}/${PETSC_ARCH}/include/petscconf.h") # > 2.3.3
      set(petsc_conf_rules "${PETSC_DIR}/conf/rules")
      set(petsc_conf_variables "${PETSC_DIR}/conf/variables")

    elseif(EXISTS "${PETSC_DIR}/bmake/${PETSC_ARCH}/petscconf.h") # <= 2.3.3
      set(petsc_conf_rules "${PETSC_DIR}/bmake/common/rules")
      set(petsc_conf_variables "${PETSC_DIR}/bmake/common/variables")

    elseif(PETSC_DIR)
      message(
        SEND_ERROR
          "The pair PETSC_DIR=${PETSC_DIR} PETSC_ARCH=${PETSC_ARCH} do not specify a valid PETSc installation"
      )
    endif()

    # Put variables into environment since they are needed to get configuration
    # (petscvariables) in the PETSc makefile
    set(ENV{PETSC_DIR} "${PETSC_DIR}")
    set(ENV{PETSC_ARCH} "${PETSC_ARCH}")

    # A temporary makefile to probe the PETSc configuration
    set(petsc_config_makefile "${PROJECT_BINARY_DIR}/Makefile.petsc")
    file(WRITE "${petsc_config_makefile}"
         "## This file was autogenerated by FindPETSc.cmake
# PETSC_DIR  = ${PETSC_DIR}
# PETSC_ARCH = ${PETSC_ARCH}
include ${petsc_conf_rules}
include ${petsc_conf_variables}
show :
\t-@echo -n \${\${VARIABLE}}
"
    )

    # define macro to get PETSc variables
    function(PETSC_GET_VARIABLE _NAME _OUTPUT)
      find_program(MAKE_EXECUTABLE NAMES make gmake)
      mark_as_advanced(FORCE MAKE_EXECUTABLE)
      execute_process(
        COMMAND ${MAKE_EXECUTABLE} --no-print-directory -f
                ${petsc_config_makefile} show VARIABLE=${_NAME}
        OUTPUT_VARIABLE _TMP
        RESULT_VARIABLE petsc_return
      )

      # split into list
      separate_arguments(_LIST UNIX_COMMAND ${_TMP})
      # return value
      set(${_OUTPUT}
          ${_LIST}
          PARENT_SCOPE
      )
    endfunction(PETSC_GET_VARIABLE)

    # get mpi bin folder
    get_filename_component(_MPI_ROOT ${MPI_Fortran_COMPILER} DIRECTORY)
    # get Fortran compiler
    petsc_get_variable(FC _PETSC_Fortran_COMPILER)
    find_program(
      PETSC_Fortran_COMPILER
      NAMES ${_PETSC_Fortran_COMPILER}
      HINTS ${PETSC_DIR} ${_MPI_ROOT} ENV ${MPI_ROOT}
    )

    # get fortran includes
    petsc_get_variable(PETSC_FC_INCLUDES PETSC_Fortran_INCLUDES)
    # get Fortran compiler flags
    petsc_get_variable(PETSC_FCPPFLAGS PETSC_Fortran_FLAGS)
    # get Fortran linker flags
    petsc_get_variable(FC_LINKER_FLAGS PETSC_Fortran_LINKER_FLAGS)

    # get C compiler
    petsc_get_variable(CC _PETSC_C_COMPILER)
    find_program(
      PETSC_C_COMPILER
      NAMES ${_PETSC_C_COMPILER}
      HINTS ${PETSC_DIR} ${_MPI_ROOT} ENV ${MPI_ROOT}
    )

    # get fortran includes
    petsc_get_variable(PETSC_CC_INCLUDES PETSC_C_INCLUDES)
    # get C compiler flags
    petsc_get_variable(PETSC_CCPPFLAGS PETSC_C_FLAGS)
    # get C linker flags
    petsc_get_variable(PCC_LINKER_FLAGS PETSC_C_LINKER_FLAGS)

    # get libraries
    petsc_get_variable(PETSC_LIB PETSC_LIBS)

    # get MPI executable
    petsc_get_variable(MPIEXEC MPIEXEC)
    find_program(
      PETSC_MPIEXEC
      NAMES ${MPIEXEC}
      HINTS ${PETSC_DIR} ${_MPI_ROOT} ENV ${MPI_ROOT}
    )

    # We are done with the temporary Makefile, calling PETSC_GET_VARIABLE after
    # this point is invalid!
    file(REMOVE ${petsc_config_makefile})

  endif()
endif()

# hand over to cmake for processing
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
  PETSC
  FOUND_VAR PETSC_FOUND
  REQUIRED_VARS PETSC_LIBS PETSC_C_COMPILER PETSC_Fortran_COMPILER PETSC_C_FLAGS
                PETSC_Fortran_FLAGS PETSC_MPIEXEC
  VERSION_VAR PETSC_VERSION
)

# hide variables
mark_as_advanced(
  FORCE PETSC_C_COMPILER PETSC_Fortran_COMPILER PETSC_LIBRARIES PETSC_MPIEXEC
  PETSC_INCLUDE_DIRS
)

# define target
if(PETSC_FOUND)
  add_library(PETSC::PETSC_Fortran INTERFACE IMPORTED)
  target_compile_options(
    PETSC::PETSC_Fortran INTERFACE ${PETSC_Fortran_FLAGS}
                                   ${PETSC_Fortan_INCLUDES}
  )
  target_link_options(
    PETSC::PETSC_Fortran INTERFACE ${PETSC_Fortran_LINKER_FLAGS}
  )
  target_include_directories(
    PETSC::PETSC_Fortran INTERFACE ${PETSC_INCLUDE_DIRS}
  )
  target_link_libraries(PETSC::PETSC_Fortran INTERFACE ${PETSC_LIBS})

  add_library(PETSC::PETSC_C INTERFACE IMPORTED)
  target_compile_options(
    PETSC::PETSC_C INTERFACE ${PETSC_C_FLAGS} ${PETSC_C_INCLUDES}
  )
  target_link_options(PETSC::PETSC_C INTERFACE ${PETSC_C_LINKER_FLAGS})
  target_include_directories(PETSC::PETSC_C INTERFACE ${PETSC_INCLUDE_DIRS})
  target_link_libraries(PETSC::PETSC_C INTERFACE ${PETSC_LIBS})
endif()
