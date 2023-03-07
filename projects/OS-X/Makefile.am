AUTOMAKE_OPTIONS = foreign

AM_CXXFLAGS = -O3 -m64 -fPIC -std=c++11 -Wfatal-errors -Wno-deprecated-declarations -mmacosx-version-min=12.6 @OPENMP_GNU@

AM_CPPFLAGS = @USE_EXTERNAL_LIBRARIES@ @COMPULSORY_INCLUDE_DIRS@ @EXTERNAL_INCLUDE_DIRS@

AM_LDFLAGS  = @BOOST_LIBS@ @DAKOTA_LIBS@ -lgfortran

#List of programs
bin_PROGRAMS = 	OptiSMOKEpp.sh

OptiSMOKEpp_sh_SOURCES = ../../source/OptiSMOKEpp.cpp
OptiSMOKEpp_sh_LDADD = @EXTERNAL_LIBS@
