# =============================================================================
#   Author: Timoteo Dinelli <timoteo.dinelli@polimi.it>
# =============================================================================
#
# Find CEQ library
#
# To provide the module with a hint about where to find your CEQ installation,
# you can set the environment variable CEQ_ROOT. The FindCEQ module will
# then look in this path when searching for CEQ paths.
# 
# At least for us the library will always have a name like 
# libceq-gcc-<gcc-version>.a for instance if <gcc-version> == 11.3.0
# this module will try to locate a library named libceq-gcc-11.3.0.a
# 
# This module will define the following variables:
#  CEQ_FOUND - true if CEQ was found on the system
#  CEQ_library - Location of the CEQ library
# message(STATUS "${LIB_NAME}")

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set(LIB_NAME libceq-clang-${CMAKE_CXX_COMPILER_VERSION}.a)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(LIB_NAME libceq-gcc-${CMAKE_CXX_COMPILER_VERSION}.a)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message(FATAL_ERROR "Contact alberto.cuoci@polimi.it or timoteo.dinelli@polimi.it")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    message(FATAL_ERROR "Contact alberto.cuoci@polimi.it or timoteo.dinelli@polimi.it")
endif()

find_library(CEQ_LIBRARIES ${LIB_NAME} PATH ENV CEQ_ROOT)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CEQ DEFAULT_MSG CEQ_LIBRARIES LIB_NAME)

mark_as_advanced(CEQ_LIBRARIES)
