
AC_INIT([OpenSMOKE++ Suite], [0.15.0], [alberto.cuoci@polimi.it])

AM_INIT_AUTOMAKE

AC_PROG_CC
AC_PROG_CXX

AC_PROG_CC([gcc icc])
AC_PROG_CXX([g++ icpc])

CXXFLAGS=''

AM_CONDITIONAL([INTEL_COMPILER],[test "$CXX" = "icpc"])

# ---------------------------------------------------------------------------------------------------------------------------------------
# Recognize if OpenMP libraries are required by the user
# ---------------------------------------------------------------------------------------------------------------------------------------
OPENMP_INTEL=""
OPENMP_GNU=""

AC_ARG_ENABLE([openmp], AS_HELP_STRING([--enable-openmp], [Enable OpenMP support]))
AS_IF([test "x$enable_openmp" = "xyes"], [OPENMP_INTEL=" -qopenmp"; OPENMP_GNU=" -fopenmp "])

AC_SUBST(OPENMP_INTEL)
AC_SUBST(OPENMP_GNU)

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add OpenSMOKE++ library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------
AC_ARG_WITH([opensmoke], [AS_HELP_STRING([--with-opensmoke=DIR], 
            [installation prefix for OpenSMOKE++ libraries])], 
            [with_opensmoke_include_path=$with_opensmoke/source],
	    AC_MSG_ERROR([Please use --with-opensmokepp to provide the path to OpenSMOKE++. Ask alberto.cuoci@polimi.it for the most recent version of OpenSMOKE++]))

AC_CHECK_FILE([${with_opensmoke_include_path}/OpenSMOKEpp], [], 
              AC_MSG_ERROR([Include file OpenSMOKEpp not found. Ask alberto.cuoci@polimi.it for the most recent version of OpenSMOKE++]))

OPENSMOKEPP_INCLUDE_DIR="-I$with_opensmoke_include_path -I../../src"

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add OpenSMOKE++ solvers (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------
#AC_ARG_WITH([opensmokeppsolvers], [AS_HELP_STRING([--with-opensmokepp-solvers=DIR], 
#            [installation prefix for OpenSMOKE++ solvers])], 
#            [with_opensmokeppsolvers_include_path=$with_opensmokeppsolvers/src],
#	    AC_MSG_ERROR([Please use --with-opensmokepp-solvers to provide the path to OpenSMOKE++ solvers. Ask alberto.cuoci@polimi.it for the most recent version of OpenSMOKE++ solvers]))

#AC_CHECK_FILE([${with_opensmokeppsolvers_include_path}/thermodynamicequilibrium/ThermodynamicEquilibrium.h], [], 
#              AC_MSG_ERROR([Include file ThermodynamicEquilibrium.h not found. Ask alberto.cuoci@polimi.it for the most recent version of OpenSMOKE++ solvers]))

#OPENSMOKEPPSOLVERS_INCLUDE_DIR="-I$with_opensmokeppsolvers_include_path/thermodynamicequilibrium"

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Eigen library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------
AC_ARG_WITH([eigen], [AS_HELP_STRING([--with-eigen=DIR], 
            [installation prefix for Eigen libraries])], 
            [with_eigen_include_path=$with_eigen],
	    AC_MSG_ERROR([Please use --with-eigen to provide the path to Eigen. Try downloading the Eigen libraries from http://eigen.tuxfamily.org/index.php?title=Main_Page]))

AC_CHECK_FILE([${with_eigen_include_path}/Eigen/Eigen], [], 
              AC_MSG_ERROR([Include file Eigen not found. Try downloading the Eigen libraries from http://eigen.tuxfamily.org/index.php?title=Main_Page]))

EIGEN_INCLUDE_DIR="-I$with_eigen_include_path"

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add RapidXML library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------
AC_ARG_WITH([rapidxml], [AS_HELP_STRING([--with-rapidxml=DIR], 
            [installation prefix for RapidXML libraries])], 
            [with_rapidxml_include_path=$with_rapidxml],
	    AC_MSG_ERROR([Please use --with-rapidxml to provide the path to RapidXML. Try downloading the RapidXML libraries from http://rapidxml.sourceforge.net/]))

AC_CHECK_FILE([${with_rapidxml_include_path}/rapidxml.hpp], [], 
              AC_MSG_ERROR([Include file rapidxml.hpp not found. Try downloading the RapidXML libraries from http://rapidxml.sourceforge.net]))

RAPIDXML_INCLUDE_DIR="-I$with_rapidxml_include_path"

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Boost C++ library (compulsory, but system implementation can be used)
# ---------------------------------------------------------------------------------------------------------------------------------------

SYSTEM_BOOST="yes"

AC_ARG_WITH([boost-include], [AS_HELP_STRING([--with-boost-include=DIR], 
            [installation prefix (include) for Boost C++ libraries])], 
            [with_boost_include_path=$with_boost_include; SYSTEM_BOOST="no"])

AC_ARG_WITH([boost-lib], [AS_HELP_STRING([--with-boost-lib=DIR], 
            [installation prefix (lib) for Boost C++ libraries])], 
            [with_boost_lib_path=$with_boost_lib; SYSTEM_BOOST="no"])


AS_IF([test "$SYSTEM_BOOST" = "yes"], [ 

	AC_CHECK_HEADERS([regex.h], [], [Include file regex.h not found. Try downloading the Boost C++ libraries from http://www.boost.org], [])

	BOOST_INCLUDE_DIR=""
	BOOST_LIBS="-lboost_system -lboost_date_time -lboost_chrono -lboost_program_options -lboost_filesystem -lboost_regex -lboost_graph";
])


AS_IF([test "$SYSTEM_BOOST" = "no"], [

	AC_CHECK_FILE([${with_boost_include_path}/boost/regex.h], [], 
	              AC_MSG_ERROR([Include file regex.h not found. Try downloading the Boost C++ libraries from http://www.boost.org]))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_date_time.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_date_time.so], [], 
		      AC_MSG_ERROR([Library file libboost_date_time.a or libboost_date_time.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_filesystem.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_filesystem.so], [], 
		      AC_MSG_ERROR([Library file libboost_filesystem.a or libboost_filesystem.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_program_options.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_program_options.so], [], 
		      AC_MSG_ERROR([Library file libboost_program_options.a or libboost_program_options.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_system.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_system.so], [], 
		      AC_MSG_ERROR([Library file libboost_system.a or libboost_system.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_regex.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_regex.so], [], 
		      AC_MSG_ERROR([Library file libboost_regex.a or libboost_regex.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_timer.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_timer.so], [], 
		      AC_MSG_ERROR([Library file libboost_timer.a or libboost_timer.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_chrono.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_chrono.so], [], 
		      AC_MSG_ERROR([Library file libboost_chrono.a or libboost_chrono.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	AC_CHECK_FILE([${with_boost_lib_path}/libboost_graph.a], [], 
		      AC_CHECK_FILE([${with_boost_lib_path}/libboost_chrono.so], [], 
		      AC_MSG_ERROR([Library file libboost_chrono.a or libboost_chrono.so not found. Try downloading the Boost C++ libraries from http://www.boost.org])))

	BOOST_INCLUDE_DIR="-I$with_boost_include_path"
	BOOST_LIBS="-L$with_boost_lib_path -lboost_regex -lboost_timer -lboost_chrono -lboost_date_time -lboost_program_options -lboost_filesystem -lboost_system"
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Intel MKL library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_MKL=" -DOPENSMOKE_USE_MKL=0 "
MKL_INCLUDE_DIR=""
MKL_LIBS=""

AC_ARG_WITH([mkl], [AS_HELP_STRING([--with-mkl=DIR], 
            [installation prefix for Intel MKL libraries])], 
            [ with_mkl_include_path=$with_mkl/include; 
	      with_mkl_lib_path=$with_mkl/lib/intel64;
              USE_MKL=" -DOPENSMOKE_USE_MKL=1 "; 
	      MKL_INCLUDE_DIR=" -I$with_mkl_include_path "; 
	      MKL_LIBS=" -Wl,--no-as-needed -L$with_mkl_lib_path -lmkl_intel_lp64 -lmkl_core -lmkl_sequential "; ] ,
             AC_MSG_ERROR([Please use --with-mkl to provide the path to Intel MKL. Try downloading the Intel MKL libraries from https://software.intel.com/en-us/intel-mkl]) )

AC_CHECK_FILE([${with_mkl_include_path}/mkl.h], [], 
	      AC_MSG_ERROR([Include file mkl.h not found. Try downloading the Intel MKL libraries from https://software.intel.com/en-us/intel-mkl]))

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add BzzMath library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------

USE_BZZMATH=" -DOPENSMOKE_USE_BZZMATH=0 "
BZZMATH_INCLUDE_DIR=""
BZZMATH_LIBS=""

AC_ARG_WITH([bzzmath-include], [AS_HELP_STRING([--with-bzzmath-include=DIR], 
            [installation prefix (include) for BzzMath libraries])], 
            [ with_bzzmath_include_path=$with_bzzmath_include;USE_BZZMATH=" -DOPENSMOKE_USE_BZZMATH=1 ";BZZMATH_INCLUDE_DIR=" -I$with_bzzmath_include_path ";

AC_CHECK_FILE([${with_bzzmath_include_path}/BzzMath.hpp], [], 
              AC_MSG_ERROR([Include file BzzMath.hpp not found. Try downloading the BzzMath libraries from http://super.chem.polimi.it/download/bzzmath-download])) ])

AC_ARG_WITH([bzzmath-lib], [AS_HELP_STRING([--with-bzzmath-lib=PATH], 
            [path to the BzzMath libraries])], 
            [ with_bzzmath_lib_path=$with_bzzmath_lib; USE_BZZMATH=" -DOPENSMOKE_USE_BZZMATH=1 ";BZZMATH_LIBS=" $with_bzzmath_lib_path ";

AC_CHECK_FILE([${with_bzzmath_lib_path}], [], 
	      AC_MSG_ERROR([Library file ${with_bzzmath_lib_path} not found. Try downloading the BzzMath libraries from http://super.chem.polimi.it/download/bzzmath-download])) ])


# ---------------------------------------------------------------------------------------------------------------------------------------
# Add CEQ library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------

CEQ_LIBS=""

AC_ARG_WITH([ceq], [AS_HELP_STRING([--with-ceq=DIR], 
            [installation prefix (include) for CEQ library])], 
            [ with_ceq_path=$with_ceq; CEQ_LIBS=" $with_ceq_path ";

AC_CHECK_FILE([${with_ceq_path}], [], 
             AC_MSG_ERROR([Library file ${with_ceq_path} not found. Try downloading the CEQ library from http://eccentric.mae.cornell.edu/~pope/CEQ/])) ])

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add Dakota library (compulsory)
# ---------------------------------------------------------------------------------------------------------------------------------------

DAKOTA_INCLUDE_DIR=""
DAKOTA_LIBS=""


AC_ARG_WITH([dakota], [AS_HELP_STRING([--with-dakota=DIR], 
            [installation prefix for DAKOTA library])], 
            [ with_dakota_path=$with_dakota; DAKOTA_LIBS=" $with_dakota_path "],
	    AC_MSG_ERROR([Please use --with-dakota to provide the path to Dakota library. Try downloading the Dakota library from https://dakota.sandia.gov/downloads]))

AC_CHECK_FILE([${with_dakota_path}], [], 
             AC_MSG_ERROR([Library path ${with_dakota_path} not found. Try downloading the Dakota library from https://dakota.sandia.gov/downloads]))

DAKOTA_INCLUDE_DIR=" -I${with_dakota_path}/include "
DAKOTA_LIBS="-L${with_dakota_path}/lib -l3po -lamplsolver -lapproxnn -lcolin -lconmin -lcport -ldakota_dll_api -ldakota_src_fortran -ldakota_src -lddace -ldfftpack -ldream -leutils -lfsudace -lhopspack -linterfaces -ljega_fe -ljega -llhs_mod -llhs_mods -llhs -lmoga -lncsuopt -lnidr -lnomad -loptpp -lpebbl -lpecos_util -lpecos_src -lpsuade -lscolib -lsoga -lsurfpack_c_interface -lsurfpack_fortran -lsurfpack_interpreter -lsurfpack -lteuchoscomm -lteuchosnumerics -lteuchosparameterlist -lteuchosparser -lteuchoscore -lteuchosremainder -ltinyxml -lutilib -lutilities";

#Compulsory libraries
AC_SUBST(OPENSMOKEPP_INCLUDE_DIR)
# AC_SUBST(OPENSMOKEPPSOLVERS_INCLUDE_DIR)
AC_SUBST(EIGEN_INCLUDE_DIR)
AC_SUBST(RAPIDXML_INCLUDE_DIR)
AC_SUBST(BOOST_INCLUDE_DIR)
AC_SUBST(BOOST_LIBS)
AC_SUBST(DAKOTA_INCLUDE_DIR)
AC_SUBST(DAKOTA_LIBS)

#Intel MKL
AC_SUBST(MKL_INCLUDE_DIR)
AC_SUBST(MKL_LIBS)
AC_SUBST(USE_MKL)

#BzzMath libraries
AC_SUBST(BZZMATH_INCLUDE_DIR)
AC_SUBST(BZZMATH_LIBS)

#CEQ
AC_SUBST(CEQ_LIBS)

USE_EXTERNAL_LIBRARIES=""
AC_SUBST(USE_EXTERNAL_LIBRARIES)

COMPULSORY_INCLUDE_DIRS="$OPENSMOKEPP_INCLUDE_DIR $EIGEN_INCLUDE_DIR $RAPIDXML_INCLUDE_DIR $BOOST_INCLUDE_DIR $DAKOTA_INCLUDE_DIR"
AC_SUBST(COMPULSORY_INCLUDE_DIRS)

EXTERNAL_INCLUDE_DIRS="$BZZMATH_INCLUDE_DIR $UMFPACK_INCLUDE_DIR $SUPERLU_INCLUDE_DIR $SUNDIALS_INCLUDE_DIR"
AC_SUBST(EXTERNAL_INCLUDE_DIRS)

EXTERNAL_LIBS="$CEQ_LIBS $BZZMATH_LIBS"
AC_SUBST(EXTERNAL_LIBS)


#
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_FILES([Makefile])

AC_OUTPUT

