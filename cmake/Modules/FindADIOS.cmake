# find the ADIOS I/O package

# search for header file in standard paths
find_path(
  ADIOS_ROOT
  NAMES include/adios.h
  HINTS ${ADIOS_ROOT} ENV ADIOS_ROOT
  DOC "Adios root folder"
)

# get include path
find_path(
  ADIOS_INCLUDE_DIRS
  NAMES adios.h
  HINTS ${ADIOS_ROOT}/include
)

# get version, set header file
if(ADIOS_INCLUDE_DIRS)
  set(_ADIOS_VERSION_H "${ADIOS_INCLUDE_DIRS}/adios_version.h")

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
    ${_ADIOS_VERSION_H} "ADIOS_VERSION_MAJOR" ADIOS_VERSION_MAJOR
  )
  _version_extract(
    ${_ADIOS_VERSION_H} "ADIOS_VERSION_MINOR" ADIOS_VERSION_MINOR
  )
  _version_extract(
    ${_ADIOS_VERSION_H} "ADIOS_VERSION_PATCH" ADIOS_VERSION_PATCH
  )

  # combine
  set(ADIOS_VERSION
      "${ADIOS_VERSION_MAJOR}.${ADIOS_VERSION_MINOR}.${ADIOS_VERSION_PATCH}"
  )

  # find the libraries TODO only Fotran for now
  if(MPI_FOUND)
    message(STATUS "Using ADIOS with MPI binding")

    find_library(
      ADIOS_LIBRARIES
      NAMES "adiosf_v1"
      HINTS ${ADIOS_ROOT}/lib
    )
  else()
    message(STATUS "Using ADIOS without MPI bindings")

    find_library(
      ADIOS_LIBRARIES
      NAMES "adiosf_nompi_v1"
      HINTS ${ADIOS_ROOT}/lib
    )
  endif()

endif()

# hand over to cmake for processing
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
  ADIOS
  FOUND_VAR ADIOS_FOUND
  REQUIRED_VARS ADIOS_INCLUDE_DIRS ADIOS_LIBRARIES
  VERSION_VAR ADIOS_VERSION
)

# hide advanced variables
mark_as_advanced(FORCE ADIOS_INCLUDE_DIRS ADIOS_LIBRARIES)

# define target
if(ADIOS_FOUND)
  add_library(ADIOS::ADIOS INTERFACE IMPORTED)
  set_target_properties(
    ADIOS::ADIOS PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                            "${ADIOS_INCLUDE_DIRS}"
  )
  target_link_libraries(ADIOS::ADIOS INTERFACE ${ADIOS_LIBRARIES} ibverbs m)
endif()

