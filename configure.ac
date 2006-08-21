dnl Process this file with autoconf to produce a configure script.
dnl Fetch the version and timestamp
m4_define([ESPRESSO_VERSION],
  [esyscmd([awk '/^\* /\
    { i=2; printf("%s", substr($i,2)); exit}'\
   RELEASE_NOTES]
  )])
m4_define([ESPRESSO_TIMESTAMP],
  ['esyscmd([awk '/[(][a-z]*[)]/\
    { for(i=2;i<NF;i++) printf("%s ",$i); gsub("[.]", "", $NF); printf $NF; exit}'\
   RELEASE_NOTES]
  )'])

# Initialise autoconf
AC_INIT(ESPResSo, [ESPRESSO_VERSION], [espresso@mpip-mainz.mpg.de])
time_stamp=ESPRESSO_TIMESTAMP

AC_PREREQ(2.50)
AC_CONFIG_SRCDIR([main.c])
AC_CONFIG_HEADER([archconfig.h])
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
	# If we are, do the following
	objdir="obj-`config/config.guess`"
	AC_MSG_RESULT([yes])
	# create an objdir if necessary
	if test ! -d "$objdir"; then
		mkdir $objdir
	fi
	AC_MSG_NOTICE([Running configure in $objdir])
	# restart configure in the objdir
	cd $objdir
	../configure --enable-wrapper "$@"
	cd ..
	exit
else
	AC_MSG_RESULT(no)
fi

# Initialise automake
AM_INIT_AUTOMAKE

AC_PREFIX_DEFAULT($HOME/Espresso)

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC

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

# Other tests

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
if test .$enable_debug=.yes; then
	AC_DEFINE_UNQUOTED(DEBUG,"DEBUG",[Turn on debugging?])
fi
# test whether to enable profiling
AC_ARG_ENABLE(profile,AC_HELP_STRING(--enable-profile,[Turn on profiling]),,)
if test .$enable_debug=.yes; then
	AC_DEFINE_UNQUOTED(PROFILING,"PROFILING",[Turn on PROFILING?])
fi

# test whether to use the wrapper or not
AC_ARG_ENABLE(wrapper,
	AC_HELP_STRING(--enable-wrapper,[Create a wrapper for the Espresso binary])
	,,)
if test .$enable_wrapper = .yes; then
   objlibdir=${libdir}/obj-`${srcdir}/config/config.guess`
fi
AM_CONDITIONAL(USE_WRAPPER, test .$enable_wrapper = .yes)
AC_SUBST(objlibdir)

# test whether to enable processor optimisation
AC_ARG_ENABLE(processor-optimization,
	AC_HELP_STRING(--enable-processor-optimization,
	[enable guessing of processor specific optimizations])
	,,enable_processor_optimization=yes)


# TODO: Replace some --enable-FEATURE by --with-PACKAGE
ES_CHECK_MPI
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

AC_DEFINE_UNQUOTED(LAST_CHANGE,"$time_stamp",[When the last change to the code was done])

# Creating the output
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

