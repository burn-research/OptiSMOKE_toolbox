# =============================================================================
#   Author: Timoteo Dinelli <timoteo.dinelli@polimi.it>
# =============================================================================
#
# Find BzzMATH library
#
# To provide the module with a hint about where to find your BazzMATH installation,
# you can set the environment variable BzzMATH_ROOT. The FindBazzMATH module will
# then look in this path when searching for  BzzMATH paths.
# 
# This module will define the following variables:
#  BzzMATH_FOUND - true if BzzMATH was found on the system
#  BzzMATH_INCLUDE_DIR - Location of the include
#  BzzMATH_LIBRARIES - Location of the library

#set(LIB_NAME libceq-gcc-${CMAKE_CXX_COMPILER_VERSION}.a)
# Set library name
set(LIB_NAME libBzzMath6_GNU.a)

# Find Include directory
find_path(BzzMATH_INCLUDE_DIR NAMES BzzMath.hpp 
    HINTS
    ENV BzzMATH_ROOT
    PATH_SUFFIXES hpp/release-${CMAKE_CXX_COMPILER_VERSION}
    )

# Find Library
find_library(BzzMATH_LIBRARIES ${LIB_NAME}
    PATH ENV BzzMATH_ROOT
    PATH_SUFFIXES lib/release-${CMAKE_CXX_COMPILER_VERSION})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(BzzMATH DEFAULT_MSG BzzMATH_LIBRARIES LIB_NAME)
find_package_handle_standard_args(BzzMATH DEFAULT_MSG BzzMATH_INCLUDE_DIR)

mark_as_advanced(BzzMATH_LIBRARIES)
mark_as_advanced(BzzMATH_INCLUDE_DIR)
