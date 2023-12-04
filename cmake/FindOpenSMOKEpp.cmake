# =============================================================================
#   Author: Timoteo Dinelli <timoteo.dinelli@polimi.it>
# =============================================================================
#
# Find OpenSMOKEpp the most famous library for handling and solving large 
# kinetic schemes. OpenSMOKEpp is header only and developed at CRECK modeling
# group by professor A. Cuoci
#
# To provide the module with a hint about where to find your OpenSMOKEpp installation,
# you can set the environment variable OpenSMOKE_ROOT. The FindOpenSMOKEpp module will
# then look in this path when searching for OpenSMOKEpp paths.
#
# This module will define the following variables:
#  OpenSMOKEpp_FOUND - true if OpenSMOKEpp was found on the system
#  OpenSMOKEpp_INCLUDE_DIRS - Location of the OpenSMOKEpp include directory
#  OpenSMOKEpp_VERSION_MAJOR - Major version
#  OpenSMOKEpp_VERSION_MINOR - Minor version
#  OpenSMOKEpp_VERSION_PATCH - Patch level
#  OpenSMOKEpp_VERSION - Full version string

if(NOT OpenSMOKE_FIND_VERSION)
  if(NOT OpenSMOKE_FIND_VERSION_MAJOR)
    set(OpenSMOKE_FIND_VERSION_MAJOR 0)
  endif()
  if(NOT OpenSMOKE_FIND_VERSION_MINOR)
    set(OpenSMOKE_FIND_VERSION_MINOR 20)
  endif()
  if(NOT OpenSMOKE_FIND_VERSION_PATCH)
    set(OpenSMOKE_FIND_VERSION_PATCH 0)
  endif()

  set(OpenSMOKE_FIND_VERSION "${OpenSMOKE_FIND_VERSION_MAJOR}.${OpenSMOKE_FIND_VERSION_MINOR}.${OpenSMOKE_FIND_VERSION_PATCH}")
endif()

macro(_opensmoke_check_version)
  file(READ "${OPENSMOKE_INCLUDE_DIR}/math/OpenSMOKEStdInclude.h" _opensmoke_version_header)

  string(REGEX MATCH "define[ \t]+__OPENSMOKE_VERSION__[ \t]+\"([0-9.]+)\""  _opensmoke_world_version_match "${_opensmoke_version_header}")
  set(OPENSMOKE_VERSION "${CMAKE_MATCH_1}")
  
  message(STATUS "Found OpenSMOKEpp VERSION: ${OPENSMOKE_VERSION}")

  if(${OPENSMOKE_VERSION} VERSION_LESS ${OpenSMOKE_FIND_VERSION})
    set(OPENSMOKE_VERSION_OK FALSE)
  else()
    set(OPENSMOKE_VERSION_OK TRUE)
  endif()

  if(NOT OPENSMOKE_VERSION_OK)
    message(STATUS "OpenSMOKE++ version ${OPENSMOKE_VERSION} found in ${OPENSMOKE_INCLUDE_DIR}, "
      "but at least version ${OpenSMOKE_FIND_VERSION} is required")
  endif()
endmacro()

if (OPENSMOKE_INCLUDE_DIR)
  # in cache already
  _opensmoke_check_version()
  set(OPENSMOKE_FOUND ${OPENSMOKE_VERSION_OK})
  set(OpenSMOKE_FOUND ${OPENSMOKE_VERSION_OK})

else ()
  
  # search first if an OpenSMOKEConfig.cmake is available in the system,
  # if successful this would set OPENSMOKE_INCLUDE_DIR and the rest of
  # the script will work as usual (For the moment we don't have yet this)
  # OpenSMOKEConfig.cmake
  find_package(OpenSMOKE ${OpenSMOKE_FIND_VERSION} NO_MODULE QUIET)

  if(NOT OPENSMOKE_INCLUDE_DIR)
    find_path(OPENSMOKE_INCLUDE_DIR OpenSMOKEpp NAMES OpenSMOKE_Definitions.h
        HINTS
        ENV OpenSMOKEpp_ROOT
        ENV OpenSMOKEpp_DIR
        ENV OpenSMOKEpp_ROOT_DIR
        PATH_SUFFIXES source
    )
  endif()

  if(OPENSMOKE_INCLUDE_DIR)
    _opensmoke_check_version()
  endif()

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(OpenSMOKEpp DEFAULT_MSG OPENSMOKE_INCLUDE_DIR OPENSMOKE_VERSION_OK)

  mark_as_advanced(OPENSMOKE_INCLUDE_DIR)

endif()

