# -----------------------------------------------------------------------------
# Preliminary Settings
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.20)

project(
  OptiSMOKEpp
  DESCRIPTION "OptiSMOKEpp Toolbox"
  LANGUAGES CXX C Fortran
  VERSION 0.3.0)

# -----------------------------------------------------------------------------
# Set module path in order to use custom CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
set(PROJECT_SOURCE_DIR ${CMAKE_SOURCE_DIR}/source)
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# -----------------------------------------------------------------------------
# Compiler options
include(${CMAKE_MODULE_PATH}/CompilerOptions.cmake)

# -----------------------------------------------------------------------------
# Options
include(${CMAKE_MODULE_PATH}/Options.cmake)

# -----------------------------------------------------------------------------
# Include and Control Fortran Interface
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

# -----------------------------------------------------------------------------
# Compulsory libraries
# -----------------------------------------------------------------------------
# Boost
find_package(
  Boost REQUIRED
  COMPONENTS date_time
             filesystem
             program_options
             system
             regex
             timer
             chrono
             serialization)

# -----------------------------------------------------------------------------
# Eigen3
find_package(Eigen3)

# -----------------------------------------------------------------------------
# OpenSMOKE++
find_package(OpenSMOKEpp REQUIRED)

# -----------------------------------------------------------------------------
# OpenSMOKESolvers
find_package(OpenSMOKEppSolvers REQUIRED)

# -----------------------------------------------------------------------------
# Dakota
find_package(
  Dakota REQUIRED
  COMPONENTS dakota_surrogates
             dakota_util
             dakota_src
             dakota_src_fortran
             nidr
             teuchosremainder
             teuchosnumerics
             teuchoscomm
             teuchosparameterlist
             teuchosparser
             teuchoscore
             pecos_util
             pecos_src
             lhs
             lhs_mods
             lhs_mod
             dfftpack
             sparsegrid
             surfpack
             surfpack
             surfpack_fortran
             utilib
             colin
             interfaces
             scolib
             3po
             pebbl
             tinyxml
             approxnn
             conmin
             ddace
             dream
             fsudace
             hopspack
             jega
             jega_fe
             moga
             soga
             eutils
             utilities
             ncsuopt
             cport
             nomad
             optpp
             psuade
             amplsolver)

# -----------------------------------------------------------------------------
# NLopt TODO: find_package(NLopt REQUIRED)

# -----------------------------------------------------------------------------
# Linear Algebra Acceleration (Recommended)
# -----------------------------------------------------------------------------
# TODO: Add MKL

# OpenBLAS
if(OPTISMOKE_USE_OPENBLAS)
  # For OpenBLAS we are not using the cmake standard defined find blas stuff because it doesnt
  # provide correctly the include directory so in order to do that we are using the one
  # automatically generated when building from source the OpenBLAS library usually located in
  # ".../lib/cmake/openblas" Handled in this way because we need to link against the static version
  find_package(OpenBLAS REQUIRED)
  message(STATUS "Found OpenBLAS --> version ${OpenBLAS_VERSION}")
endif()

# --------------------------------------------------------------------
# Targets to compile
# --------------------------------------------------------------------
add_executable(OptiSMOKEpp ${PROJECT_SOURCE_DIR}/OptiSMOKEpp.cpp)

target_include_directories(
  OptiSMOKEpp PUBLIC ${PROJECT_SOURCE_DIR} ${OPENSMOKE_INCLUDE_DIR} ${OpenBLAS_INCLUDE_DIRS}
                     ${OPENSMOKEPPSOLVERS_INCLUDE_DIR})

target_link_libraries(
  OptiSMOKEpp
  PUBLIC Eigen3::Eigen
         Boost::date_time
         Boost::filesystem
         Boost::program_options
         Boost::system
         Boost::regex
         Boost::timer
         Boost::chrono
         ${Dakota_LIBRARIES}
         ${OpenBLAS_INCLUDE_DIRS}/../lib/libopenblas.a)

target_compile_definitions(
  OptiSMOKEpp
  PUBLIC OPENSMOKE_USE_MKL=0
         OPENSMOKE_USE_OPENBLAS=1
         OPENSMOKE_USE_SUNDIALS=0
         OPENSMOKE_USE_UMFPACK=0
         OPENSMOKE_USE_BZZMATH=0
         OPENSMOKE_USE_SUPERLU_SERIAL=0
         OPENSMOKE_USE_RADAU=0
         OPENSMOKE_USE_ODEPACK=0
         OPENSMOKE_USE_DVODE=0
         OPENSMOKE_USE_DASPK=0
         OPENSMOKE_USE_MEBDF=0)

# --------------------------------------------------------------------
# Installing options
# --------------------------------------------------------------------

# TODO:Add installation instructions

# --------------------------------------------------------------------
# Testing options
# --------------------------------------------------------------------

# TODO:Add automatic testing
