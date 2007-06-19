dnl Process this file with autoconf to produce a configure script.
dnl Fetch the version and timestamp from RELEASE_NOTES
m4_define([ESPRESSO_VERSION],
  [esyscmd([awk '/^\* /\
    { i=2; printf("%s", substr($i,2)); exit}'\
   RELEASE_NOTES]
  )])
m4_define([ESPRESSO_TIMESTAMP],
  [esyscmd([awk '/[(][a-z]*[)]/\
    { i=1; if ($i == "*") { next; }; for(i=2;i<NF;i++) printf("%s ",$i); gsub("[.]", "", $NF); printf $NF; exit}'\
   RELEASE_NOTES]
  )])

# Initialise autoconf
AC_INIT(ESPResSo, [ESPRESSO_VERSION], [espresso@mpip-mainz.mpg.de])

AC_PREREQ(2.50)
AC_CONFIG_SRCDIR([main.c])
AC_CONFIG_AUX_DIR(config)
AC_CONFIG_MACRO_DIR(config)
AC_PREFIX_DEFAULT($HOME/Espresso)
AC_CANONICAL_HOST
AC_CANONICAL_TARGET

# Test if we are in the source directory
AC_MSG_CHECKING([whether configure is running in the source directory])
if test "$srcdir" = "."; then
	objdir="obj-`config/config.guess`"
	AC_MSG_RESULT([yes])
	# create an objdir if necessary
	if test ! -d "$objdir"; then
		mkdir $objdir
	fi
	AC_MSG_NOTICE([Running configure in $objdir])
	# restart configure in the objdir
	cd $objdir
	../configure --enable-chooser "$@"
	ec=$?
	cd ..
	exit $ec
else
	AC_MSG_RESULT(no)
fi

# Initialise automake
AM_INIT_AUTOMAKE([1.9 foreign subdir-objects])
AM_MAINTAINER_MODE
# Make the release notes a dependency of configure
AC_SUBST([CONFIGURE_DEPENDENCIES], ['$(top_srcdir)/RELEASE_NOTES'])

ES_KNOWN_CONFIGS

cat <<EOF
****************************************************************
*                   Check the compiler and MPI                 *
****************************************************************
EOF
save_CFLAGS=$CFLAGS
AC_PROG_CC
AC_PROG_CC_C_O
AM_PROG_CC_C_O
CFLAGS=$save_CFLAGS
ES_CHECK_MPI
ES_CHECK_COMPILER

cat <<EOF
****************************************************************
*                   Check for headers and libraries            *
****************************************************************
EOF
# Checks for libraries.
ES_CHECK_EFENCE
ES_CHECK_TCL
ES_CHECK_TK
ES_CHECK_FFTW

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([arpa/inet.h fcntl.h float.h inttypes.h limits.h netdb.h\
 netinet/in.h stddef.h stdlib.h string.h strings.h sys/file.h sys/socket.h\
 sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MEMCMP
AC_FUNC_REALLOC
AC_FUNC_SELECT_ARGTYPES
AC_TYPE_SIGNAL
AC_TYPE_SIZE_T
AC_FUNC_STRTOD
AC_CHECK_FUNCS([atexit floor getcwd gethostbyname memset pow rint select socket sqrt])

cat <<EOF
****************************************************************
*                   Other settings                             *
****************************************************************
EOF
ES_INIT_CHOOSER
ES_INIT_MYCONFIG
AC_DEFINE(LAST_CHANGE,["ESPRESSO_TIMESTAMP"],[When the last change to the code was done])
AC_ARG_VAR(LIBS,[use this to specify additional libraries to link against, e.g. -lmx])

# Checks for programs
AC_PROG_AWK
AC_PROG_RANLIB
AC_PROG_LN_S
AC_PATH_PROG(DOXYGEN,doxygen,no)
AM_CONDITIONAL(HAVE_DOXYGEN,test .$DOXYGEN != .no)
AC_PATH_PROG(PDFLATEX,pdflatex,no)
AM_CONDITIONAL(HAVE_PDFLATEX,test .$PDFLATEX != .no)

cat <<EOF
****************************************************************
*                   Creating output files                      *
****************************************************************
EOF
# Creating the output
AC_CONFIG_HEADERS([acconfig.h])
AC_CONFIG_LINKS(Makefile)
AC_CONFIG_FILES([
	Makefile-am scripts/Makefile testsuite/Makefile
	doc/Makefile doc/ug/Makefile doc/dg/Makefile
	config/$myconfig:config/myconfig-default.h.in
	chooser/Espresso mpiwrap/Espresso
	])
AC_CONFIG_FILES([Espresso],[chmod 755 Espresso])

AC_OUTPUT

cat <<EOF
****************************************************************
*                   Configuration summary                      *
****************************************************************
platform detected	= $target

Compiler settings:
------------------
MPI			= $with_mpi
compiler		= $CC
linker			= $LD
c compiler flags	= $CFLAGS
c preprocessor flags	= $CPPFLAGS
linker flags		= $LDFLAGS
libraries		= $LIBS
using debug flags	= $enable_debug
using profiling flags	= $enable_profiling
inline			= $MDINLINE
EOF

if test .$with_mpi = .generic; then
   echo "WARNING: a generic MPI version may or may not work, look at the Espresso script in $builddir"
fi
if test "$issue_fake_warning"; then cat <<EOF
********************************************************************
* Could neither detect LAM nor MPICH nor a native MPI environment, *
* using the FAKE implementation for one processor only.            *
* If you have an MPI environment, please specify its type and      *
* the compiler or includes/libraries of your MPI implementation    *
* manually, or, even better, add your MPI environment to           * 
* config/mpi.m4                                                    *
********************************************************************
EOF
fi

cat <<EOF

Libraries:
----------
Tcl version		= $use_tcl
Tk version		= $use_tk
FFTW version		= $use_fftw
efence			= $with_efence
EOF

if test .$use_tcl != .tcl8.5 && test .$use_tcl != .tcl8.4; then
   echo "WARNING: Not all supplied scripts work properly with Tcl version before 8.4, consider upgrading"
fi 

cat <<EOF

Other settings:
---------------
using Espresso chooser	= $enable_chooser
local config header	= $myconfig
pdflatex		= $PDFLATEX
doxygen			= $DOXYGEN
EOF

if test .$pdflatex = .no; then
   echo "To be able to rebuild the user's guide, add pdflatex to your PATH."
fi
if test .$doxygen = .no; then
   echo "To be able to build the developer's guide, add doxygen to your PATH."
fi
