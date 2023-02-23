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
#  OpenSMOKEpp_INCLUDE_DIRS - Location of the OpenSMOKEpp includes
#  OpenSMOKEpp_VERSION_MAJOR - Major version
#  OpenSMOKEpp_VERSION_MINOR - Minro version
#  OpenSMOKEpp_VERSION_PATCH - Patch level
#  OpenSMOKEpp_VERSION - Full version string

SET(Open_SMOKE_INCLUDE_SEARCH_PATHS
    /usr/include
    /usr/include/OpenSMOKEpp
    /usr/include/opensmokepp
    /usr/local/include
    /usr/local/include/OpenSMOKEpp
    /usr/local/include/opensmokepp
    /opt/OpenSMOKEpp/source
    ${OpenSMOKEpp_USER_PATHS}
    ${OpenSMOKEpp_ROOT}
    ${OpenSMOKEpp_ROOT_DIR}
    ${OpenSMOKEpp_DIR}
)

# Find the OpenSMOKEpp include directories
find_path(OpenSMOKEpp_INCLUDE_DIR OpenSMOKEpp
    PATHS
        ${Open_SMOKE_INCLUDE_SEARCH_PATHS}
    ENV
        OpenSMOKEpp_ROOT
        OpenSMOKEpp_DIR
        OpenSMOKEpp_ROOT_DIR
    PATH_SUFFIXES
        source
    NAMES    
        OpenSMOKE_Definitions.h
)
set(OpenSMOKEpp_INCLUDE_DIR
    "${OpenSMOKEpp_INCLUDE_DIR}"
)

# Exctract the version
# TODO

mark_as_advanced(
    OpenSMOKEpp_INCLUDE_DIR
)
