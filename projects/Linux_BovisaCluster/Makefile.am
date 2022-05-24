AUTOMAKE_OPTIONS = foreign

if INTEL_COMPILER
AM_CXXFLAGS = -O3 -m64 -fPIC -std=c++11 @OPENMP_INTEL@
else
AM_CXXFLAGS = -O3 -m64 -fPIC -std=c++11 -Wfatal-errors @OPENMP_GNU@
endif

AM_CPPFLAGS = 	@USE_MKL@ @USE_EXTERNAL_LIBRARIES@ \
		@COMPULSORY_INCLUDE_DIRS@ @EXTERNAL_INCLUDE_DIRS@ @MKL_INCLUDE_DIR@

AM_LDFLAGS  =   @MKL_LIBS@ @BOOST_LIBS@ @DAKOTA_LIBS@ -lgfortran

#List of programs
bin_PROGRAMS = 	OptiSMOKEpp.sh

OptiSMOKEpp_sh_SOURCES = ../../src/OptiSMOKEpp.cpp ../../src/DirectApplicInterface.cpp
OptiSMOKEpp_sh_LDADD = @EXTERNAL_LIBS@
