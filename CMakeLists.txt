cmake_minimum_required (VERSION 3.16)

project (OptiSMOKEpp
  DESCRIPTION "OptiSMOKEpp Suite"
  LANGUAGES CXX C Fortran
  VERSION 0.3.0
)

# --- Set module path in order to use custom CMake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(PROJECT_SOURCE_DIR "${CMAKE_SOURCE_DIR}/source")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
# set(CMAKE_SKIP_BUILD_RPATH TRUE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# --- Controlling fortran c++ interoperability
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

#-----------------------------------------------------------------------------
# Compulsory libraries
#-----------------------------------------------------------------------------

# --- Boost (Compulsory)
find_package (Boost REQUIRED
  COMPONENTS date_time filesystem program_options system regex timer chrono
  serialization
)

# --- Eigen (Compulsory)
find_package (Eigen3 REQUIRED NO_MODULE)

# --- OpenSMOKE++ (Compulsory)
find_package(OpenSMOKEpp REQUIRED)

# --- OpenSMOKESolvers (Compulsory)
find_package(OpenSMOKEppSolvers REQUIRED)

# --- NLopt++
find_package(NLopt REQUIRED)

# --- Dakota (Optional)
find_package(Dakota REQUIRED
  PATHS $ENV{Dakota_ROOT}
  COMPONENTS dakota_surrogates dakota_util dakota_src dakota_src_fortran nidr
  teuchosremainder teuchosnumerics teuchoscomm teuchosparameterlist
  teuchosparser teuchoscore pecos_util pecos_src lhs lhs_mods lhs_mod dfftpack
  sparsegrid surfpack surfpack surfpack_fortran utilib colin interfaces scolib
  3po pebbl tinyxml approxnn conmin ddace dream fsudace hopspack jega jega_fe
  moga soga eutils utilities ncsuopt cport nomad optpp psuade amplsolver
)
  
message(STATUS "Dakota include:   ${Dakota_INCLUDE_DIRS}")
message(STATUS "Dakota libraries: ${Dakota_LIBRARIES}")
message(STATUS "Dakota binaries:  ${Dakota_EXECUTABLE}")

set(COMPULSORY_LIBRARIES
    ${Boost_LIBRARIES}
    ${NLOPT_LIBRARIES}
    ${Dakota_LIBRARIES}
)

set(COMPULSORY_INCLUDE
    ${Boost_INCLUDE_DIRS}
    ${OPENSMOKE_INCLUDE_DIR}
    ${OPENSMOKEPPSOLVERS_INCLUDE_DIR}
    ${NLOPT_INCLUDE_DIRS}
    ${DAKOTA_INCLUDE_DIRS}
)

#-----------------------------------------------------------------------------
# In principle the following libraries should be optional since for instance
# CEQ and blockstdg in the external subdirectory are only for 1D flames
#-----------------------------------------------------------------------------

option(OPTISMOKE_USE_MKL "Activate intel MKL support for OpenSMOKEpp" FALSE)
option(OPTISMOKE_USE_OPENBLAS "Activate OpenBLAS support for OpenSMOKEpp" FALSE)

find_package(CEQ REQUIRED)
message (STATUS "CEQ library: ${CEQ_LIBRARIES}")

set(LINEAR_ALGEBRA_INCLUDE)
set(LINEAR_ALGEBRA_LIB)
if(OPTISMOKE_USE_MKL AND OPTISMOKE_USE_OPENBLAS)
    message(FATAL_ERROR "Solvers can be compiled with just one BLAS
    distribution at a time! Choose between MKL and OpenBLAS")
endif()

# --- Intel MKL (Optional but strongly rcommended: choose one between MKL=BLAS and OpenBLAS)
if (OPTISMOKE_USE_MKL)
    set (BLA_VENDOR Intel10_64lp_seq)
    find_package (BLAS REQUIRED COMPONENTS mkl_core mkl_sequential)
    message (STATUS "MKL include directory: $ENV{BLAS_ROOT}/include")
    message (STATUS "MKL library: ${BLAS_LIBRARIES}")
    list(APPEND LINEAR_ALGEBRA_INCLUDE $ENV{BLAS_ROOT}/include)
    list(APPEND LINEAR_ALGEBRA_LIB ${BLAS_LIBRARIES})
endif()

# --- OpenBLAS (Optional but strongly Recommended: choose one between MKL and OpenBLAS)
if(OPTISMOKE_USE_OPENBLAS)  
    set (BLA_VENDOR OpenBLAS)
    set (BLA_STATIC ON)
    find_package(BLAS REQUIRED)
    message(STATUS "OpenBLAS include directory: $ENV{BLAS_ROOT}/include")
    message(STATUS "OpenBLAS library: ${BLAS_LIBRARIES}")
    list(APPEND LINEAR_ALGEBRA_INCLUDE  $ENV{BLAS_ROOT}/include)
    list(APPEND LINEAR_ALGEBRA_LIB ${BLAS_LIBRARIES})
endif()

add_subdirectory(external)

#--------------------------------------------------------------------
# Compilation flags and options
#--------------------------------------------------------------------

option(OPTISMOKE_USE_OPENMP "Activate support for OpenMP" FALSE)
if (OPTISMOKE_USE_OPENMP)
    add_compile_options(-fopenmp)
    add_link_options(-fopenmp)
endif()

add_compile_options(-fPIC -m64 -w)
add_link_options(-fPIC -m64 -w)

#--------------------------------------------------------------------
# Targets to compile
#--------------------------------------------------------------------

add_executable (OptiSMOKEpp.sh ${PROJECT_SOURCE_DIR}/OptiSMOKEpp.cpp)

#--------------------------------------------------------------------
# Include flags
#--------------------------------------------------------------------

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(FORTRAN_LIB_NAME -static-libgfortran)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(FORTRAN_LIB_NAME gfortran)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    message(FATAL_ERROR "Contact timoteo.dinelli@polimi.it")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    message(FATAL_ERROR "Contact timoteo.dinelli@polimi.it")
endif()

target_include_directories(OptiSMOKEpp.sh
    PUBLIC
    ${PROJECT_SOURCE_DIR}
    ${COMPULSORY_INCLUDE}
    ${LINEAR_ALGEBRA_INCLUDE}
)

target_link_libraries(OptiSMOKEpp.sh
    PUBLIC
    Eigen3::Eigen
    ${COMPULSORY_LIBRARIES} 
    ${LINEAR_ALGEBRA_LIB}
    ${CEQ_LIBRARIES}
    ${CEQ_LIBRARIES}
    blocktdsge
    gfortran
)
    # ${FORTRAN_LIB_NAME}

target_compile_definitions(OptiSMOKEpp.sh 
    PUBLIC
    OPENSMOKE_USE_MKL=0
    OPENSMOKE_USE_OPENBLAS=1
    OPENSMOKE_USE_SUNDIALS=0
    OPENSMOKE_USE_UMFPACK=0
    OPENSMOKE_USE_BZZMATH=0
    OPENSMOKE_USE_SUPERLU_SERIAL=0
    OPENSMOKE_USE_RADAU=0
    OPENSMOKE_USE_ODEPACK=0
    OPENSMOKE_USE_DVODE=0
    OPENSMOKE_USE_DASPK=0
    OPENSMOKE_USE_MEBDF=0
)

#--------------------------------------------------------------------
# Installing options
#--------------------------------------------------------------------

install(TARGETS OptiSMOKEpp.sh DESTINATION ${INSTALL_BIN_DIR})

#--------------------------------------------------------------------
# Testing options
#--------------------------------------------------------------------

# TODO
