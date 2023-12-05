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

# search first if an OpenSMOKEConfig.cmake is available in the system,
# if successful this would set OPENSMOKE_INCLUDE_DIR and the rest of
# the script will work as usual (For the moment we don't have yet this)
# OpenSMOKEConfig.cmake

find_package(OpenSMOKEppSolvers NO_MODULE QUIET)
if(NOT OPENSMOKEPPSOLVERS_INCLUDE_DIR)
    find_path(OPENSMOKEPPSOLVERS_INCLUDE_DIR OpenSMOKEppSolvers NAMES batchreactor/BatchReactor.cpp
        HINTS
        ENV OpenSMOKEppSolvers_ROOT
        ENV OpenSMOKEppSolvers_DIR
        ENV OpenSMOKEppSolvers_ROOT_DIR
        PATH_SUFFIXES src
    )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OpenSMOKEppSolvers DEFAULT_MSG OPENSMOKEPPSOLVERS_INCLUDE_DIR)

mark_as_advanced(OPENSMOKEPPSOLVERS_INCLUDE_DIR)
