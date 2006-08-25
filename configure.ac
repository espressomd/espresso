dnl Process this file with autoconf to produce a configure script.
dnl Fetch the version and timestamp
m4_define([ESPRESSO_VERSION],
  [esyscmd([awk '/^\* /\
    { i=2; printf("%s", substr($i,2)); exit}'\
   RELEASE_NOTES]
  )])
m4_define([ESPRESSO_TIMESTAMP],
  ['esyscmd([awk '/[(][a-z]*[)]/\
    { i=1; if ($i == "*") { next; }; for(i=2;i<NF;i++) printf("%s ",$i); gsub("[.]", "", $NF); printf $NF; exit}'\
   RELEASE_NOTES]
  )'])

# Initialise autoconf
AC_INIT(ESPResSo, [ESPRESSO_VERSION], [espresso@mpip-mainz.mpg.de])
time_stamp=ESPRESSO_TIMESTAMP

AC_PREREQ(2.50)
AC_CONFIG_SRCDIR([main.c])
AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_HOST
AC_CANONICAL_TARGET

dnl Using the object directory
dnl m4_divert_push([HELP_BEGIN])dnl
dnl AC_HELP_STRING(--objdir=DIR,[where to put the object files])
dnl m4_divert_pop([HELP_BEGIN])dnl

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
AM_INIT_AUTOMAKE

AC_PREFIX_DEFAULT($HOME/Espresso)

# Checks for programs.
AC_PROG_AWK

save_CFLAGS=$CFLAGS
AC_PROG_CC
CFLAGS=$save_CFLAGS

# Checks for libraries.

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
AC_FUNC_STRTOD
AC_CHECK_FUNCS([atexit floor getcwd gethostbyname memset pow rint select socket sqrt])


# Test for (obsolete) option --enable-mode
AC_ARG_ENABLE(mode,[],
	[
	AC_MSG_RESULT([\
********************************************************************************
* The option --enable-mode is obsolete and will be removed in a future version *
* of Espresso! Please use --enable-debug or --enable-profiling instead!        *
********************************************************************************\
	])
	if test .$enable_mode = .debug; then
	   enable_debug=yes
	fi
	if test .$enable_mode = .profiling; then
	   enable_profiling=yes
	fi
	],)


# test whether to enable debugging
AC_ARG_ENABLE(debug,AC_HELP_STRING(--enable-debug,[Turn on debugging]),,)
if test .$enable_debug = .yes; then
   AC_DEFINE(DEBUG,,[Turn on debugging?])
fi

# test whether to enable profiling
AC_ARG_ENABLE(profiling,AC_HELP_STRING(--enable-profiling,[Turn on profiling]),,)
if test .$enable_profiling = .yes; then
   AC_DEFINE(PROFILING,,[Turn on profiling?])
fi

# test whether to enable processor optimisation
AC_ARG_ENABLE(processor-optimization,
	AC_HELP_STRING(--enable-processor-optimization,
	[enable guessing of processor specific optimizations])
	,,enable_processor_optimization=yes)

# ESPRESSO Tests

# TODO: Replace some --enable-FEATURE by --with-PACKAGE
ES_KNOWN_CONFIGS

ES_CHECK_MPI
AM_CONDITIONAL(MPI_FAKE, test .$enable_mpi = .fake)
AC_SUBST(MPI_INVOCATION)

ES_CHECK_EFENCE
ES_CHECK_INLINING
if test .$enable_processor_optimization = .yes; then
        CF_CHECK_OPTFLAGS
fi
ES_CHECK_LINK
ES_CHECK_TCL
ES_CHECK_TK
ES_CHECK_FFTW
AC_CHECK_PROG(doxygen, doxygen, yes, no)

# test whether to use the chooser or not
AC_ARG_ENABLE(chooser,
	AC_HELP_STRING(--enable-chooser,[Create a chooser for the Espresso binary])
	,,)
AM_CONDITIONAL(USE_CHOOSER, test .$enable_chooser = .yes)
if test .$enable_chooser = .yes; then
   execpkglibdir='${pkglibdir}'/obj-`$srcdir/config/config.guess`
   scriptsdir='${pkglibdir}/scripts'
else
   scriptsdir='${pkgdatadir}/scripts'
fi
AC_SUBST(scriptsdir)
AC_SUBST(execpkglibdir)

AC_DEFINE_UNQUOTED(LAST_CHANGE,"$time_stamp",[When the last change to the code was done])
AC_ARG_VAR(LIBS,[use this to specify additional libraries to link against, e.g. -lmx])

dnl Handling the myconfig-header
AC_ARG_WITH(myconfig, 
	AC_HELP_STRING(--with-myconfig=FILE,[name of the local config file [[myconfig.h]]]),
	[ test .$with_myconfig != .no && MYCONFIG_H=$with_myconfig ],
	[ MYCONFIG_H=myconfig.h ])

AC_SUBST(MYCONFIG_H)
AC_DEFINE_UNQUOTED(MYCONFIG_H,"$MYCONFIG_H",[The name of the local configure header])
AH_BOTTOM([
#ifdef MYCONFIG_H
#include MYCONFIG_H
#endif
])

# Creating the output
AC_CONFIG_HEADER([acconfig.h])
AC_CONFIG_LINKS(Makefile)
AC_CONFIG_FILES(
	[Makefile-am scripts/Makefile testsuite/Makefile doc/Makefile] 
	[config/$MYCONFIG_H:config/myconfig-default.h.in]
	[chooser/Espresso mpiwrap/Espresso],
	)
AC_CONFIG_FILES([Espresso],[chmod 755 Espresso])

AC_OUTPUT

