#--------------------------------------------------------------------
# Preliminary Settings
#--------------------------------------------------------------------

cmake_minimum_required (VERSION 3.16)

project (OptiSMOKEpp
    DESCRIPTION
    "OptiSMOKEpp Suite"
    LANGUAGES 
    C CXX Fortran)

# --- Set module path in order to use custom CMake modules

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set (PROJECT_SOURCE_DIR "${CMAKE_SOURCE_DIR}/source")

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)
set (CMAKE_COLOR_MAKEFILE ON)
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
set (VERSION 0.3.0)
set (INSTALL_BIN_DIR "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Installation directory for executables")

if( NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release)
endif()

message (STATUS "Compiling in mode: ${CMAKE_BUILD_TYPE}")
message (STATUS "OptiSMOKE++ directory: ${PROJECT_SOURCE_DIR}")

# --- Controlling fortran c++ interoperability
include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

#--------------------------------------------------------------------
# Compulsory libraries
#--------------------------------------------------------------------

# --- Boost (Compulsory)
find_package (Boost REQUIRED COMPONENTS date_time filesystem program_options system regex timer chrono serialization)
message (STATUS "Boost include directory: ${Boost_INCLUDE_DIRS}")
message (STATUS "Boost libraries: ${Boost_LIBRARIES}")

# --- Eigen (Compulsory)
find_package (Eigen3 REQUIRED)

# --- OpenSMOKE++ (Compulsory)
find_package(OpenSMOKEpp REQUIRED)
message (STATUS "OpenSMOKE++ include directory: ${OPENSMOKE_INCLUDE_DIR}")

# --- OpenSMOKESolvers (Compulsory)
find_package(OpenSMOKEppSolvers REQUIRED)
message (STATUS "OpenSMOKESolvers include directory: ${OPENSMOKEPPSOLVERS_INCLUDE_DIR}")

# --- NLopt++
find_package(NLopt REQUIRED)
message(STATUS "NLopt include directory: ${NLOPT_INCLUDE_DIRS}")
message(STATUS "NLopt libraries: ${NLOPT_LIBRARIES}")

# --- Dakota (Optional)
find_package(Dakota REQUIRED PATHS $ENV{Dakota_ROOT} COMPONENTS dakota_surrogates dakota_util dakota_src dakota_src_fortran nidr 
  teuchosremainder teuchosnumerics teuchoscomm teuchosparameterlist teuchosparser
  teuchoscore pecos_util pecos_src lhs lhs_mods lhs_mod dfftpack sparsegrid
  surfpack surfpack surfpack_fortran utilib colin interfaces scolib 3po
  pebbl tinyxml approxnn conmin ddace dream fsudace hopspack jega jega_fe
  moga soga eutils utilities ncsuopt cport nomad optpp psuade amplsolver)
  
message(STATUS "Dakota include:   ${Dakota_INCLUDE_DIRS}")
message(STATUS "Dakota libraries: ${Dakota_LIBRARIES}")
message(STATUS "Dakota binaries:  ${Dakota_EXECUTABLE}")

set(COMPULSORY_LIBRARIES
    ${Boost_LIBRARIES}
    ${NLOPT_LIBRARIES}
    ${Dakota_LIBRARIES})

set(COMPULSORY_INCLUDE
    ${Boost_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
    ${OPENSMOKE_INCLUDE_DIR}
    ${OPENSMOKEPPSOLVERS_INCLUDE_DIR}
    ${NLOPT_INCLUDE_DIRS}
    ${DAKOTA_INCLUDE_DIRS})


# Since at the moment we are not including any flame or what we dont need
# this shit
# --- libblocktdsge.a (Compulsory)
# add_subdirectory(${PROJECT_SOURCE_DIR}/numerics)

#--------------------------------------------------------------------
# Semi-Compulsory libraries (compulsory just for some solvers)
#--------------------------------------------------------------------

# --- CEQ (Compulsory just for some solvers)
option (OPTISMOKE_USE_CEQ "Activate CEQ support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_CEQ)
    find_package(CEQ REQUIRED)
    message (STATUS "CEQ library: ${CEQ_LIBRARIES}")
endif()

set(LINEAR_ALGEBRA_INCLUDE)
set(LINEAR_ALGEBRA_LIB)

if(OPTISMOKE_USE_MKL AND OPTISMOKE_USE_OPENBLAS)
    message(FATAL_ERROR "Solvers can be compiled with just one BLAS distribution at a time! Choose between MKL and OpenBLAS")
endif()

# --- Intel MKL (Optional but strongly rcommended: choose one between MKL=BLAS and OpenBLAS)
option (OPTISMOKE_USE_MKL "Activate MKL/BLAS support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_MKL)
    set (BLA_VENDOR Intel10_64lp_seq)
    find_package (BLAS REQUIRED COMPONENTS mkl_core mkl_sequential)
    message (STATUS "MKL include directory: $ENV{BLAS_ROOT}/include")
    message (STATUS "MKL library: ${BLAS_LIBRARIES}")
    list(APPEND LINEAR_ALGEBRA_INCLUDE $ENV{BLAS_ROOT}/include)
    list(APPEND LINEAR_ALGEBRA_LIB ${BLAS_LIBRARIES})
endif()

# --- OpenBLAS (Optional but strongly Recommended: choose one between MKL and OpenBLAS)
option (OPTISMOKE_USE_OPENBLAS "Activate OpenBLAS support for OpenSMOKEpp" FALSE)
if(OPTISMOKE_USE_OPENBLAS)  
    set (BLA_VENDOR OpenBLAS)
    set (BLA_STATIC ON)
    find_package(BLAS REQUIRED)
    message(STATUS "OpenBLAS include directory: ${BLAS_ROOT}/include")
    message(STATUS "OpenBLAS library: ${BLAS_LIBRARIES}")
    list(APPEND LINEAR_ALGEBRA_INCLUDE  ${BLAS_ROOT}/include)
    list(APPEND LINEAR_ALGEBRA_LIB ${BLAS_LIBRARIES})
endif()

#--------------------------------------------------------------------
# Optional libraries
#--------------------------------------------------------------------

set(OPTIONAL_INCLUDE)
set(OPTIONAL_LIBRARIES)

# --- SUNDIALS (Optional)
option (OPTISMOKE_USE_SUNDIALS "Activate SUNDIALS support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_SUNDIALS)
    find_package (SUNDIALS REQUIRED)
    message (STATUS "SUNDIALS include directory: ${SUNDIALS_INCLUDE_DIRS}")
    message (STATUS "SUNDIALS libraries:  ${SUNDIALS_LIBRARIES}")
    list(APPEND OPTIONAL_INCLUDE ${SUNDIALS_INCLUDE_DIRS})
    list(APPEND OPTIONAL_LIBRARIES ${SUNDIALS_LIBRARIES})
endif()

# --- UMFPACK (Optional) also referred to be SuiteSparse (TODO)
option (OPTISMOKE_USE_UMFPACK "Activate UMFPACK support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_UMFPACK)
    find_package (SuiteSparse REQUIRED COMPONENTS umfpack amd klu cholmod ccolamd camd Config)
    message (STATUS "Found UMFPACK version: ${SuiteSparse_VERSION}")

    message (STATUS "UMFPACK include directory: ${SuiteSparse_umfpack_INCLUDE_DIR}")
    # I leave here the version with separates include cause in the future maybe the behaviour will
    # change however for now the include is a unique directory
    #message (STATUS "UMFPACK include directory: ${SuiteSparse_umfpack_INCLUDE_DIR},
    #${SuiteSparse_amd_INCLUDE_DIR}, ${SuiteSparse_klu_INCLUDE_DIR}, ${SuiteSparse_cholmod_INCLUDE_DIR},
    #${SuiteSparse_ccolamd_INCLUDE_DIR}, ${SuiteSparse_camd_INCLUDE_DIR}, ${SuiteSparse_Config_INCLUDE_DIR}")
    
    message (STATUS "UMFPACK libraries:\n\t${SuiteSparse_umfpack_LIBRARY},\n\t${SuiteSparse_amd_LIBRARY},\n\t${SuiteSparse_klu_LIBRARY},\n\t${SuiteSparse_cholmod_LIBRARY},\n\t${SuiteSparse_ccolamd_LIBRARY},\n\t${SuiteSparse_camd_LIBRARY},\n\t${SuiteSparse_Config_LIBRARY}")
    
    list(APPEND OPTIONAL_INCLUDE ${SuiteSparse_umfpack_INCLUDE_DIR})

    list(APPEND OPTIONAL_LIBRARIES ${SuiteSparse_umfpack_LIBRARY} ${SuiteSparse_amd_LIBRARY}
        ${SuiteSparse_klu_LIBRARY} ${SuiteSparse_cholmod_LIBRARY} ${SuiteSparse_ccolamd_LIBRARY}
        ${SuiteSparse_camd_LIBRARY} ${SuiteSparse_Config_LIBRARY})
endif()

# --- BzzMath
option (OPTISMOKE_USE_BZZMATH "Activate BzzMath support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_BZZMATH)
    find_package (BzzMATH REQUIRED COMPONENTS bzzmath)
    message (STATUS "BzzMATH include directory: ${BzzMATH_INCLUDE_DIR}")
    message (STATUS "BzzMATH library: ${BzzMATH_LIBRARIES}")
    list(APPEND OPTIONAL_INCLUDE ${BzzMATH_INCLUDE_DIR})
    list(APPEND OPTIONAL_LIBRARIES ${BzzMATH_LIBRARIES})
endif()

# --- SuperLU-serial (Optional)
option (OPTISMOKE_USE_SUPERLU_SERIAL "Activate SuperLU-serial support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_SUPERLU_SERIAL)
    find_package (SuperLU REQUIRED)
    message (STATUS "SuperLU include directory: ${SUPERLU_INCLUDE_DIRS}")
    message (STATUS "SuperLU libraries: ${SUPERLU_LIBRARIES}")
    list(APPEND OPTIONAL_INCLUDE ${SUPERLU_INCLUDE_DIRS})
    list(APPEND OPTIONAL_LIBRARIES ${SUPERLU_LIBRARIES})
endif()

# --- RADAU (Optional)
option (OPTISMOKE_USE_RADAU "Activate RADAU support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_RADAU)
    find_library (RADAU_LIBRARIES libRADAU_LAPACK_gcc_64bit.a  PATH $ENV{RADAU_ROOT}/lib)
    message (STATUS "RADAU library: ${RADAU_LIBRARIES}")
    list(APPEND OPTIONAL_LIBRARIES ${RADAU_LIBRARIES})
endif()

# --- ODEPACK (Optional)
option (OPTISMOKE_USE_ODEPACK "Activate ODEPACK support for OpenSMOKEpp" FALSE)
option (OPTISMOKE_USE_DVODE "Activate DVODE support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_ODEPACK)
    find_library (ODEPACK_LIBRARIES libODEPACK_LAPACK_gcc_64bit.a PATH $ENV{ODEPACK_ROOT}/lib)
    message (STATUS "ODEPACK library: ${ODEPACK_LIBRARIES}")
    list(APPEND OPTIONAL_LIBRARIES ${ODEPACK_LIBRARIES})
    list(APPEND USE_OPTIONAL_LIBRARIES OPTISMOKE_USE_ODEPACK)
    option (OPTISMOKE_USE_DVODE "Activate DVODE support for OpenSMOKEpp" TRUE)
endif()

# --- DASPK (Optional)
option (OPTISMOKE_USE_DASPK "Activate DASPK support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_DASPK)
    find_library (DASPK_LIBRARIES libDASPK20_LAPACK_gcc_64bit.a PATH $ENV{DASPK_ROOT}/lib)
    message (STATUS "DASPK library: ${DASPK_LIBRARIES}")
    list(APPEND OPTIONAL_LIBRARIES ${DASPK_LIBRARIES})
endif()

# --- MEBDF (Optional)
option (OPTISMOKE_USE_MEBDF "Activate MEBDF support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_MEBDF)
    find_library (MEBDF_LIBRARIES libMEBDF_LAPACK_gcc_64bit.a PATH $ENV{MEBDF_ROOT}/lib)
    message (STATUS "MEBDF library: ${MEBDF_LIBRARIES}")
    list(APPEND OPTIONAL_LIBRARIES ${MEBDF_LIBRARIES})
endif()

# --- LINPACK  (Optional)
option (OPTISMOKE_USE_LINPACK "Activate LINPACK support for OpenSMOKEpp" FALSE)
if (OPTISMOKE_USE_LINPACK)
    find_library (LINPACK_LIBRARIES libLINPACK_gcc_64bit.a PATH $ENV{LINPACK_ROOT}/lib)
    message (STATUS "LINPACK library: ${LINPACK_LIBRARIES}")
    list(APPEND OPTIONAL_LIBRARIES ${LINPACK_LIBRARIES})
endif()

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
    "${PROJECT_SOURCE_DIR}"
    "${COMPULSORY_INCLUDE}"
    "${LINEAR_ALGEBRA_INCLUDE}"
    "${OPTIONAL_INCLUDE}"
)
target_link_libraries(OptiSMOKEpp.sh
    PUBLIC
    ${COMPULSORY_LIBRARIES} 
    ${LINEAR_ALGEBRA_LIB}
    ${OPTIONAL_LIBRARIES}
    ${FORTRAN_LIB_NAME}
)

target_compile_definitions(OptiSMOKEpp.sh 
    PUBLIC
    OPENSMOKE_USE_MKL=${OPTISMOKE_USE_MKL}
    OPENSMOKE_USE_OPENBLAS=${OPTISMOKE_USE_OPENBLAS}
    OPENSMOKE_USE_SUNDIALS=${OPTISMOKE_USE_SUNDIALS}
    OPENSMOKE_USE_UMFPACK=${OPTISMOKE_USE_UMFPACK}
    OPENSMOKE_USE_BZZMATH=${OPTISMOKE_USE_BZZMATH}
    OPENSMOKE_USE_SUPERLU_SERIAL=${OPTISMOKE_USE_SUPERLU_SERIAL}
    OPENSMOKE_USE_RADAU=${OPTISMOKE_USE_RADAU}
    OPENSMOKE_USE_ODEPACK=${OPTISMOKE_USE_ODEPACK}
    OPENSMOKE_USE_DVODE=${OPTISMOKE_USE_DVODE}
    OPENSMOKE_USE_DASPK=${OPTISMOKE_USE_DASPK}
    OPENSMOKE_USE_MEBDF=${OPTISMOKE_USE_MEBDF}
)

#--------------------------------------------------------------------
# Installing options
#--------------------------------------------------------------------

install(TARGETS OptiSMOKEpp.sh DESTINATION ${INSTALL_BIN_DIR})

#--------------------------------------------------------------------
# Testing options
#--------------------------------------------------------------------

# TODO
