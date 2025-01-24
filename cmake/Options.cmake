# List containing the options activated during the configuration of the CMake project
set(PROGRAM_OPTIONS)

option(OPTISMOKE_USE_MKL "Activate intel MKL support for OpenSMOKEpp" FALSE)
option(OPTISMOKE_USE_OPENBLAS "Activate OpenBLAS support for OpenSMOKEpp" FALSE)

if(OPTISMOKE_USE_MKL AND OPTISMOKE_USE_OPENBLAS)
  message(
    FATAL_ERROR
      "Solvers can be compiled with just one BLAS distribution at a time! Choose between MKL and OpenBLAS"
  )
endif()

# TODO:Testing option(BUILD_TEST "Build the c++ tests for the library. Rquires: GTest" OFF)
