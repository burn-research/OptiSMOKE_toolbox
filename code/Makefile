#########################################################################################################################################################################################################
# Specify source directory path
SRCDIR = /home/chimica2/bertolino/OptSpp_versions/Stable_OptiSMOKE_V3
# Specify intallation path
PREFIX = /home/chimica2/bertolino/OptSpp_versions/Stable_OptiSMOKE_V3

# Compiler
CC = g++
CFLAGS = -std=c++11

# Use compiler optimization with specified level
OLVL = -O3

# Use OpenMP
OPENMP = -fopenmp

#########################################################################################################################################################################################################
# Compusolry libraries
PATHBOOST = /software/chimica2/libraries/boost/boost-1.64.0-gcc-6.3.0
PATHDAKOTA = /software/chimica2/libraries/dakota/dakota-6.5
PATHEIGEN = -I/software/chimica2/libraries/eigen/eigen-3.3.4
PATHOPENSMOKE = -IOpenSMOKEpp-master/source
PATHRAPIDXML = -I/software/chimica2/libraries/rapidxml/rapidxml-1.13
# Specific libraries from DAKOTA
# Note that these can change depending on the version of DAKOTA that you use!
LDFLAGSDAKOTA = -l3po -lamplsolver -lapproxnn -lcolin -lconmin -lcport -ldakota_dll_api -ldakota_src_fortran -ldakota_src -lddace -ldfftpack -ldream -leutils -lfsudace -lhopspack -linterfaces -ljega_fe -ljega -llhs_mod -llhs_mods -llhs -lmoga -lncsuopt -lnidr -lnomad -loptpp -lpebbl -lpecos -lpecos_src -lpsuade -lscolib -lsoga -lsurfpack_c_interface -lsurfpack_fortran -lsurfpack_interpreter -lsurfpack -lteuchos -ltinyxml -lutilib -lutilities
LDFLAGSBOOST = -lboost_regex -lboost_timer -lboost_chrono -lboost_date_time -lboost_program_options -lboost_filesystem -lboost_system
#########################################################################################################################################################################################################
# Use MKL libraries
MKL = -DOPENSMOKE_USE_MKL=1
ifdef MKL
# Need to specify path to mkl libraries
PATHINTEL = /software/chimica2/libraries/intel/xe_2018
PATHMKL = -I$(PATHINTEL)/include -L$(PATHINTEL)/lib/intel64
LDFLAGSMKL = -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
endif
#########################################################################################################################################################################################################
# Optional libraries
PATHBZZMATH = /software/chimica2/libraries/bzzmath/bzzmath-6.0
PATHCEQ = /software/chimica2/libraries/ceq
#PATHDASPK = /software/chimica2/libraries/daspk/daspk-gcc-6.3.0/lib
PATHLIBCONFIG = /software/chimica2/libraries/libconfig/libconfig-1.5-gcc-6.3.0
#PATHLINPACK = /software/chimica2/libraries/linpack/linpack-gcc-6.3.0/lib
#PATHMEDF = /software/chimica2/libraries/mebdf/medf-gcc-6.3.0/lib
#PATHODEPACK = /software/chimica2/libraries/odepack/odepack-gcc-6.3.0/lib
#PATHRADAU = /software/chimica2/libraries/radau/radau-gcc-6.3.0/lib
PATHSUITEPARSE = /software/chimica2/libraries/suitesparse/suitesparse-5.1.0-gcc-6.3.0
#PATHSUNDIALS = /software/chimica2/libraries/sundials/sundials-3.1.0-opensmoke++/serial-gcc-4.8.5
PATHSUPERLU = /software/chimica2/libraries/superlu/superlu/superlu-5.2-gcc-6.3.0
#########################################################################################################################################################################################################
# Combine the paths for the optional libraries
OPTIONALPATHS =
ifdef PATHBZZMATH
OPTIONALPATHS += -DOPENSMOKE_USE_BZZMATH=1 -I$(PATHBZZMATH)/hpp/release-6.3.0/ -L$(PATHBZZMATH)/lib/release-6.3.0 -lBzzMath6_GNU
endif
ifdef PATHCEQ
OPTIONALPATHS += $(PATHCEQ)/libceq-gcc-6.3.0.a
endif
ifdef PATHDASPK
OPTIONALPATHS += -DOPENSMOKE_USE_DASPK=1 $(PATHDASPK)/libDASPK20_LAPACK_gcc_64bit.a
endif
ifdef PATHLIBCONFIG
OPTIONALPATHS += -I$(PATHLIBCONFIG)/include $(PATHLIBCONFIG)/lib/libconfig++.a
endif
ifdef PATHLINPACK
OPTIONALPATHS += $(PATHLINPACK)/libLINPACK_gcc_64bit.a
endif
ifdef PATHMEDF
OPTIONALPATHS += -DOPENSMOKE_USE_MEBDF=1 $(PATHMEDF)/libMEBDF_LAPACK_gcc_64bit.a
endif
ifdef PATHODEPACK
OPTIONALPATHS += -DOPENSMOKE_USE_ODEPACK=1 -DOPENSMOKE_USE_DVODE=1 $(PATHODEPACK)/libODEPACK_LAPACK_gcc_64bit.a
endif
ifdef PATHRADAU
OPTIONALPATHS += -DOPENSMOKE_USE_RADAU=1 $(PATHRADAU)/libRADAU_LAPACK_gcc_64bit.a
endif
ifdef PATHSUITEPARSE
OPTIONALPATHS += -DOPENSMOKE_USE_UMFPACK=1 -I$(PATHSUITEPARSE)/include -L$(PATHSUITEPARSE)/lib -lklu -lumfpack -lamd -lcholmod -lccolamd -lcamd -lsuitesparseconfig -lmetis
endif
ifdef PATHSUNDIALS
OPTIONALPATHS += -DOPENSMOKE_USE_SUNDIALS=1 -I$(PATHSUNDIALS)/include -L$(PATHSUNDIALS)/lib -lsundials_ida -lsundials_kinsol -lsundials_cvode -lsundials_sunlinsollapackband -lsundials_sunlinsollapackdense -lsundials_sunmatrixband -lsundials_sunmatrixdense -lsundials_nvecserial -lgfortran
endif
ifdef PATHSUPERLU
OPTIONALPATHS += -DOPENSMOKE_USE_SUPERLU_SERIAL=1 -I$(PATHSUPERLU)/include $(PATHSUPERLU)/lib/libSuperLU52_gcc_64bit.a
endif
#########################################################################################################################################################################################################
all:
	$(CC) $(CFLAGS) $(OLVL) $(OPENMP) $(MKL) $(SRCDIR)/OptiSMOKEpp.cpp $(SRCDIR)/DirectApplicInterface.cpp $(PATHOPENSMOKE) $(PATHRAPIDXML) $(PATHEIGEN) -I$(PATHBOOST)/include -L$(PATHBOOST)/lib -I$(PATHDAKOTA)/include -L$(PATHDAKOTA)/lib $(LDFLAGSDAKOTA) $(LDFLAGSBOOST) $(PATHMKL) $(LDFLAGSMKL) $(OPTIONALPATHS) -o OptiSMOKEpp.sh

install: OptiSMOKEpp.sh
	mkdir -p $(DESTDIR)$(PREFIX)/bin
	cp --remove-destination $< $(DESTDIR)$(PREFIX)/bin/OptiSMOKEpp.sh
clean:
	rm -rf *.o OptiSMOKEpp.sh $(DESTDIR)$(PREFIX)/bin/OptiSMOKEpp.sh
#########################################################################################################################################################################################################
