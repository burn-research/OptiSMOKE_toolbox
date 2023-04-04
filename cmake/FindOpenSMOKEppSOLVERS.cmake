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

SET(Open_SMOKE_SOLVER_INCLUDE_SEARCH_PATHS
    /usr/include
    /usr/include/OpenSMOKEppSolvers
    /usr/local/include/OpenSMOKEppSolvers
    /opt/OpenSMOKEppSolvers/source
    ${OpenSMOKEppSolvers_USER_PATHS}
    ${OpenSMOKEppSOLVERS_ROOT}
    ${OpenSMOKEppSOLVERS_ROOT_DIR}
    ${OpenSMOKEppSOLVERS_DIR}
)

# Find the OpenSMOKEpp include directories
find_path(OpenSMOKEppSOLVERS_INCLUDE_DIR OpenSMOKEppSolvers
    PATHS
        ${Open_SMOKE_SOLVER_INCLUDE_SEARCH_PATHS}
    ENV
		OpenSMOKEppSOLVERS_ROOT
        OpenSMOKEppSOLVERS_DIR
        OpenSMOKEppSOLVERS_ROOT_DIR
    PATH_SUFFIXES
        src
    NAMES
        batchreactor/Grammar_BatchReactor.h
)
set(OpenSMOKEppSOLVERS_INCLUDE_DIR
    "${OpenSMOKEppSOLVERS_INCLUDE_DIR}"
)

# Exctract the version
# TODO

mark_as_advanced(
    OpenSMOKEppSOLVERS_INCLUDE_DIR
)
