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

# --- OpenSmoke++ (Compulsory)
find_package(OpenSMOKEpp REQUIRED)
message (STATUS "OpenSmoke++: ${OpenSMOKEpp_INCLUDE_DIR}")

# --- libconfig++ (Compulsory for OpenSMOKE++)
#find_package (CONFIG REQUIRED COMPONENTS libconfig++)
#message (STATUS "libCONFIG++ include directory: ${CONFIG++_INCLUDE_DIR}")
#message (STATUS "libCONFIG++ libraries: ${CONFIG++_LIBRARY}")

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
                       ${OpenSMOKEpp_INCLUDE_DIR})

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

# --- SUNDIALS (Optional)
option (OPTISMOKE_USE_SUNDIALS "Activate SUNDIALS support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_SUNDIALS)
  find_package (SUNDIALS QUIET COMPONENTS sundials_cvode sundials_kinsol sundials_ida sundials_sunlinsollapackband sundials_sunlinsollapackdense)
  message (STATUS "SUNDIALS_INCLUDE_DIRS: ${SUNDIALS_INCLUDE_DIRS}")
  message (STATUS "SUNDIALS_LIBRARIES: ${SUNDIALS_LIBRARIES}")

  list(APPEND OPTIONAL_INCLUDE ${SUNDIALS_INCLUDE_DIRS})
  list(APPEND OPTIONAL_LIBRARIES ${SUNDIALS_LIBRARIES})
  add_compile_definitions(OPENSMOKE_USE_SUNDIALS=1)
endif()

# --- UMFPACK
option (OPTISMOKE_USE_UMFPACK "Activate UMFPACK support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_UMFPACK)
  find_package (UMFPACK QUIET COMPONENTS umfpack)
  message (STATUS "UMFPACK_INCLUDES: ${UMFPACK_INCLUDES}")
  message (STATUS "UMFPACK_LIBRARIES: ${UMFPACK_LIBRARIES}")

  list(APPEND OPTIONAL_INCLUDE ${UMFPACK_INCLUDES})
  list(APPEND OPTIONAL_LIBRARIES ${UMFPACK_LIBRARIES})
  add_compile_definitions(OPENSMOKE_USE_UMFPACK=1)
endif()

# --- BzzMath TODO
option (OPTISMOKE_USE_BZZMATH "Activate BzzMath support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_BZZMATH)
  find_package (BZZMATH QUIET COMPONENTS bzzmath)
  message (STATUS "BZZMATH_INCLUDES: ${BZZMATH_INCLUDES}")
  message (STATUS "BZZMATH_LIBRARIES: ${BZZMATH_LIBRARIES}")

  list(APPEND OPTIONAL_INCLUDE ${BZZMATH_INCLUDES})
  list(APPEND OPTIONAL_LIBRARIES ${BZZMATH_LIBRARIES})
  add_compile_definitions(OPENSMOKE_USE_BZZMATH=1)
endif()

# --- SuperLU-serial (Optional)
option (OPTISMOKE_USE_SUPERLU_SERIAL "Activate SuperLU-serial support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_SUPERLU_SERIAL)
  find_package (SuperLU QUIET)
  message (STATUS "SUPERLU_INCLUDES: ${SUPERLU_INCLUDES}")
  message (STATUS "SUPERLU_LIBRARIES: ${SUPERLU_LIBRARIES}")

  list(APPEND OPTIONAL_INCLUDE ${SUPERLU_INCLUDES})
  list(APPEND OPTIONAL_LIBRARIES ${SUPERLU_LIBRARIES})
  add_compile_definitions(OPENSMOKE_USE_SUPERLU_SERIAL=1)
endif()

# --- RADAU (Optional)
option (OPTISMOKE_USE_RADAU "Activate RADAU support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_RADAU)
  find_library (RADAU_LIBRARIES libRADAU_LAPACK_gcc_64bit.a  PATH $ENV{RADAU_ROOT}/lib)
  message (STATUS "RADAU_LIBRARIES: ${RADAU_LIBRARIES}")
  list(APPEND OPTIONAL_LIBRARIES ${RADAU_LIBRARIES})
  add_compile_definitions(OPENSMOKE_USE_RADAU=1)
endif()

# --- ODEPACK (Optional)
option (OPTISMOKE_USE_ODEPACK "Activate ODEPACK support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_ODEPACK)
  find_library (ODEPACK_LIBRARIES libODEPACK_LAPACK_gcc_64bit.a PATH $ENV{ODEPACK_ROOT}/lib)
  message (STATUS "ODEPACK_LIBRARIES: ${ODEPACK_LIBRARIES}")
  list(APPEND OPTIONAL_LIBRARIES ${ODEPACK_LIBRARIES})
  add_compile_definitions(OPENSMOKE_USE_ODEPACK=1)
  add_compile_definitions(OPENSMOKE_USE_DVODE=1)
endif()

# --- DASPK (Optional)
option (OPTISMOKE_USE_DASPK "Activate DASPK support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_DASPK)
  find_library (DASPK_LIBRARIES libDASPK20_LAPACK_gcc_64bit.a PATH $ENV{DASPK_ROOT}/lib)
  message (STATUS "DASPK_LIBRARIES: ${DASPK_LIBRARIES}")
  list(APPEND OPTIONAL_LIBRARIES ${DASPK_LIBRARIES})
  add_compile_definitions(OPENSMOKE_USE_DASPK=1)
endif()

# --- MEBDF (Optional)
option (OPTISMOKE_USE_MEBDF "Activate MEBDF support for OpenSMOKEpp" FALSE)
if (OPENSMOKE_USE_MEBDF)
  find_library (MEBDF_LIBRARIES libMEBDF_LAPACK_gcc_64bit.a PATH $ENV{MEBDF_ROOT}/lib)
  message (STATUS "MEBDF_LIBRARIES: ${MEBDF_LIBRARIES}")
  list(APPEND OPTIONAL_LIBRARIES ${MEBDF_LIBRARIES})
  add_compile_definitions(OPENSMOKE_USE_MEBDF=1)
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
