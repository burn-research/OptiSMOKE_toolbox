#--------------------------------------------------------------------
# Preliminary Settings
#--------------------------------------------------------------------

cmake_minimum_required (VERSION 3.20)

project (OptiSMOKEpp)

# --- Set module path in order to use custom CMake modules
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} $ENV{EIGEN_LIBRARY_PATH}/cmake)
set (CMAKE_COLOR_MAKEFILE ON)
set (VERSION "3.0.0-Beta")
set (INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation directory for executables")
set (MACOSX_DEPLOYMENT_TARGET 12.6)

if( NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif()

message (STATUS "Compiling in mode: ${CMAKE_BUILD_TYPE}")

#--------------------------------------------------------------------
# Compulsory libraries
#--------------------------------------------------------------------

# --- Source code for the project
message (STATUS "OptiSMOKE++ source: ${PROJECT_SOURCE_DIR}")

# --- Boost (Compulsory)
find_package (Boost REQUIRED COMPONENTS date_time filesystem program_options system regex timer chrono serialization)
message (STATUS "Boost include directory: ${Boost_INCLUDE_DIRS}")
message (STATUS "Boost libraries: ${Boost_LIBRARIES}")

# --- Eigen (Compulsory)
find_package (Eigen3 REQUIRED)
message (STATUS "Eigen libraries: ${EIGEN3_INCLUDE_DIR}")

# --- OpenSMOKE++ (Compulsory)
find_package(OpenSMOKEpp REQUIRED)
message (STATUS "OpenSMOKE++: ${OpenSMOKEpp_INCLUDE_DIR}")

# --- OpenSMOKEppSolvers (Compulsory)
find_package(OpenSMOKEppSOLVERS REQUIRED)
message (STATUS "OpenSMOKEppSOLVERS: ${OpenSMOKEppSOLVERS_INCLUDE_DIR}")

# --- CEQ (Compulsory for OpenSMOKE++)
option (OPTISMOKE_USE_CEQ "Activate CEQ support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_CEQ)
  find_library (CEQ_LIBRARIES libceq.a PATH $ENV{CEQ_ROOT})
  message (STATUS "CEQ library: ${CEQ_LIBRARIES}")
  add_compile_definitions(OPENSMOKE_USE_CEQ=1)
endif()

# --- Dakota (Compulsory)
find_package(Dakota REQUIRED PATHS $ENV{Dakota_ROOT} COMPONENTS dakota_surrogates dakota_util dakota_src dakota_src_fortran nidr 
  teuchosremainder teuchosnumerics teuchoscomm teuchosparameterlist teuchosparser
  teuchoscore pecos_util pecos_src lhs lhs_mods lhs_mod dfftpack sparsegrid
  surfpack surfpack surfpack_fortran utilib colin interfaces scolib 3po
  pebbl tinyxml approxnn conmin ddace dream fsudace hopspack jega jega_fe
  moga soga eutils utilities ncsuopt cport nomad optpp psuade amplsolver)
  
message(STATUS "Dakota include:   ${Dakota_INCLUDE_DIRS}")
message(STATUS "Dakota libraries: ${Dakota_LIBRARIES}")
message(STATUS "Dakota binaries:  ${Dakota_EXECUTABLE}")

set(COMPULSORY_LIBRARIES ${Boost_LIBRARIES} 
                         ${Dakota_LIBRARIES}
                         ${CEQ_LIBRARIES})

set(COMPULSORY_INCLUDE ${Boost_INCLUDE_DIRS} 
                       ${Dakota_INCLUDE_DIRS}
                       ${EIGEN3_INCLUDE_DIR}
                       ${OpenSMOKEpp_INCLUDE_DIR}
                       ${OpenSMOKEppSOLVERS_INCLUDE_DIR})

#--------------------------------------------------------------------
# Semi-Compulsory libraries (desc: TODO)
#--------------------------------------------------------------------

set(SEMI_COMPULSORY_INCLUDE)
set(SEMI_COMPULSORY_LIBRARIES)

# --- Intel MKL (Optional but strongly rcommended: choose one between MKL=BLAS and OpenBLAS)
option (OPTISMOKE_USE_MKL "Activate MKL/BLAS support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_MKL)
  set (BLA_VENDOR Intel10_64lp)
  find_package (BLAS)
  message (STATUS "MKL_LIBRARIES: ${BLAS_LIBRARIES}")
  message (STATUS "MKL_DIR: ${BLAS_DIR}")
  list(APPEND SEMI_COMPULSORY_INCLUDE ${BLAS_DIR})
  list(APPEND SEMI_COMPULSORY_LIBRARIES ${BLAS_LIBRARIES})
  add_compile_definitions(OPENSMOKE_USE_MKL=1)
endif()

# --- OpenBLAS (Optional but strongly Recommended: choose one between MKL and OpenBLAS)
option (OPTISMOKE_USE_OPENBLAS "Activate OpenBLAS support for OpenSMOKEpp" FALSE)
if(OPTISMOKE_USE_OPENBLAS)  
  find_package(OpenBLAS)
  message(STATUS "OpenBLAS_LIBRARIES: ${OpenBLAS_LIBRARIES}")
  message(STATUS "OpenBLAS_DIR: ${OpenBLAS_INCLUDE_DIR}")
  list(APPEND SEMI_COMPULSORY_INCLUDE ${OpenBLAS_INCLUDE_DIR})
  list(APPEND SEMI_COMPULSORY_LIBRARIES ${OpenBLAS_LIBRARIES})
  add_compile_definitions(OPENSMOKE_USE_OPENBLAS=1)
endif()

#--------------------------------------------------------------------
# Optional libraries (desc: TODO)
#--------------------------------------------------------------------

set(OPTIONAL_INCLUDE)
set(OPTIONAL_LIBRARIES)

option (OPTISMOKE_USE_NLOPT "Activate NLOPT++ as optimization library" FALSE)
if(OPTISMOKE_USE_NLOPT)
	if(LINUX)
		find_package(NLopt)
	else()
		find_package(NLOPT)
	endif()
	message(STATUS "NLOPT_DIR: ${NLOPT_INCLUDE_DIRS}")
	message(STATUS "NLOPT_LIBRARIES: ${NLOPT_LIBRARIES}")
	list(APPEND OPTIONAL_INCLUDE ${NLOPT_INCLUDE_DIRS})
	list(APPEND OPTIONAL_LIBRARIES ${NLOPT_LIBRARIES})
	add_compile_definitions(OPTISMOKE_USE_NLOPT=1)
endif()

#--------------------------------------------------------------------
# Compilation flags and options
#--------------------------------------------------------------------

# Suppress warnings
option(OPTISMOKE_USE_OPENMP "Activate support for OpenMP" FALSE)
if (OPTISMOKE_USE_OPENMP)
  add_compile_options(-fopenmp)
  add_link_options(-fopenmp)
endif()

add_compile_options(-std=c++17 -fPIC -m64 -w -lgfortran)
add_link_options(-std=c++17 -fPIC -m64 -w -lgfortran)

#--------------------------------------------------------------------
# Include flags
#--------------------------------------------------------------------

include_directories (${PROJECT_SOURCE_DIR}/source
                     ${COMPULSORY_INCLUDE}
                     ${SEMI_COMPULSORY_INCLUDE}
                     ${OPTIONAL_INCLUDE})

#--------------------------------------------------------------------
# Library flags
#--------------------------------------------------------------------

link_libraries (${COMPULSORY_LIBRARIES}
                ${SEMI_COMPULSORY_LIBRARIES}
                ${OPTIONAL_LIBRARIES})

#--------------------------------------------------------------------
# Targets to compile
#--------------------------------------------------------------------

add_executable (OptiSMOKEpp.sh source/OptiSMOKEpp.cpp)

#--------------------------------------------------------------------
# Installing options
#--------------------------------------------------------------------

install(TARGETS OptiSMOKEpp.sh DESTINATION ${INSTALL_BIN_DIR})

#--------------------------------------------------------------------
# Testing options
#--------------------------------------------------------------------

# enable_testing()
