
AC_INIT([OptiSMOKE++], [0.3.0], [timoteo.dinelli@polimi.it])

AM_INIT_AUTOMAKE([subdir-objects])

AC_PROG_CC
AC_PROG_CXX

AC_PROG_CC([gcc-12])
AC_PROG_CXX([g++-12])

CXXFLAGS=''

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

OPENSMOKEPP_INCLUDE_DIR="-I$with_opensmoke_include_path -I../../source"

# ---------------------------------------------------------------------------------------------------------------------------------------
# Add OpenSMOKE++ solvers (compulsory) in future they will be compulsory 
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
	BOOST_LIBS="-lboost_regex -lboost_timer -lboost_chrono -lboost_date_time -lboost_program_options -lboost_filesystem -lboost_system";
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
	BOOST_LIBS="-L$with_boost_lib_path -lboost_regex -lboost_timer -lboost_chrono -lboost_date_time -lboost_program_options -lboost_filesystem -lboost_system -lboost_serialization"
])

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
#DAKOTA_LIBS="-L${with_dakota_path}/lib -l3po -lamplsolver -lapproxnn -lcolin -lconmin -lcport -ldakota_dll_api -ldakota_src_fortran -ldakota_src -lddace -ldfftpack -ldream -leutils -lfsudace -lhopspack -linterfaces -ljega_fe -ljega -llhs_mod -llhs_mods -llhs -lmoga -lncsuopt -lnidr -lnomad -loptpp -lpebbl -lpecos -lpecos_src -lpsuade -lscolib -lsoga -lsurfpack_c_interface -lsurfpack_fortran -lsurfpack_interpreter -lsurfpack -lteuchos -ltinyxml -lutilib -lutilities";
DAKOTA_LIBS="-L${with_dakota_path}/lib -lpecos_util -lpecos_src -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchosparser -lteuchoscore -ldakota_src -ldakota_src_fortran -lnidr -llhs -ldfftpack -lsparsegrid -lsurfpack -lsurfpack -lsurfpack_fortran -lutilib -lcolin -linterfaces -lscolib -l3po -lpebbl -ltinyxml -lconmin -lfsudace -lhopspack -ljega -ljega_fe -lmoga -lsoga -leutils -lutilities -lncsuopt -lcport -lnomad -loptpp -lpsuade -lamplsolver";

#Compulsory libraries
AC_SUBST(OPENSMOKEPP_INCLUDE_DIR)
# AC_SUBST(OPENSMOKEPPSOLVERS_INCLUDE_DIR)
AC_SUBST(EIGEN_INCLUDE_DIR)
AC_SUBST(BOOST_INCLUDE_DIR)
AC_SUBST(BOOST_LIBS)
AC_SUBST(DAKOTA_INCLUDE_DIR)
AC_SUBST(DAKOTA_LIBS)

USE_EXTERNAL_LIBRARIES=""
AC_SUBST(USE_EXTERNAL_LIBRARIES)

COMPULSORY_INCLUDE_DIRS="$OPENSMOKEPP_INCLUDE_DIR $EIGEN_INCLUDE_DIR $BOOST_INCLUDE_DIR $DAKOTA_INCLUDE_DIR"
AC_SUBST(COMPULSORY_INCLUDE_DIRS)

EXTERNAL_INCLUDE_DIRS=""
AC_SUBST(EXTERNAL_INCLUDE_DIRS)

EXTERNAL_LIBS=""
AC_SUBST(EXTERNAL_LIBS)

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
