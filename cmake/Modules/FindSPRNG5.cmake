# find a suitable SPRNG5 install

# search for header file in standard paths
find_path(
  SPRNG5_ROOT
  NAMES include/sprng.h
  HINTS ${SPRNG5_ROOT} ENV SPRNG5_ROOT
  DOC "SPRNG5 root folder"
)

# get include path
find_path(
  SPRNG5_INCLUDE_DIRS
  NAMES "sprng.h"
  HINTS "${SPRNG5_ROOT}/include"
)

# get version, set header file
if(SPRNG5_INCLUDE_DIRS)
  set(_SPRNG5_VERSION_H "${SPRNG5_INCLUDE_DIRS}/sprng.h")

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

  # TODO SPRNG5 does not export any version information
  set(SPRNG5_VERSION_MAJOR 5)
  set(SPRNG5_VERSION_MINOR 0)
  set(SPRNG5_VERSION_PATCH 0)

  # extract version _VERSION_EXTRACT(${_SPRNG5_VERSION_H} "SPRNG5_VERSION_MAJOR"
  # SPRNG5_VERSION_MAJOR) _VERSION_EXTRACT(${_SPRNG5_VERSION_H}
  # "SPRNG5_VERSION_MINOR" SPRNG5_VERSION_MINOR)
  # _VERSION_EXTRACT(${_SPRNG5_VERSION_H} "SPRNG5_VERSION_PATCH"
  # SPRNG5_VERSION_PATCH)

  # combine
  set(SPRNG5_VERSION
      "${SPRNG5_VERSION_MAJOR}.${SPRNG5_VERSION_MINOR}.${SPRNG5_VERSION_PATCH}"
  )

  # find the libraries
  find_library(
    SPRNG5_LIBRARIES
    NAMES "sprng"
    HINTS ${SPRNG5_ROOT}/lib
  )

endif()

# hand over to cmake for processing
include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(
  SPRNG5
  FOUND_VAR SPRNG5_FOUND
  REQUIRED_VARS SPRNG5_INCLUDE_DIRS SPRNG5_LIBRARIES
  VERSION_VAR SPRNG5_VERSION
)

# make variables advanced
mark_as_advanced(FORCE SPRNG5_INCLUDE_DIRS SPRNG5_LIBRARIES)

# define target
if(SPRNG5_FOUND)
  add_library(SPRNG5::SPRNG5 INTERFACE IMPORTED)
  set_target_properties(
    SPRNG5::SPRNG5 PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                              "${SPRNG5_INCLUDE_DIRS}"
  )
  target_link_libraries(SPRNG5::SPRNG5 INTERFACE ${SPRNG5_LIBRARIES})
endif()

