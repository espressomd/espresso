dnl Process this file with autoconf to produce a configure script.
dnl Fetch the version and timestamp
m4_define([ESPRESSO_VERSION],
  [esyscmd([awk '/^\* /\
    { for(i=2;i<NF-1;i++) printf("%s ",$i); printf $(NF-1); exit}'\
   RELEASE_NOTES]
  )])
m4_define([ESPRESSO_TIMESTAMP],
  ['esyscmd([awk '/[(][a-z]*[)]/\
    { for(i=2;i<NF;i++) printf("%s ",$i); gsub("[.]", "", $NF); printf $NF; exit}'\
   RELEASE_NOTES]
  )'])

# Initialise autoconf
AC_INIT(ESPResSo, ESPRESSO_VERSION, [espresso@mpip-mainz.mpg.de])
time_stamp=ESPRESSO_TIMESTAMP

AC_PREREQ(2.50)
AC_CONFIG_SRCDIR([$unique_file])
AC_CONFIG_HEADER([archconfig.h])
AC_CONFIG_AUX_DIR(config)

# First test if we are in the source directory
AC_MSG_CHECKING([whether configure is running in the source directory])
if test "$srcdir" = "."; then
	# If we are, do the following
	objdir="obj-`config/config.guess`"
	AC_MSG_RESULT([yes])
	# create an objdir if necessary
	if test ! -d "$objdir"; then
		mkdir $objdir
	fi
	AC_MSG_NOTICE([Starting configure from $objdir])
	# restart configure in the objdir
	cd $objdir
	../configure $@
	cd ..
	exit
else
	AC_MSG_RESULT(no)
fi

AM_INIT_AUTOMAKE

AC_PREFIX_DEFAULT($HOME/Espresso)

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([arpa/inet.h fcntl.h float.h inttypes.h limits.h netdb.h netinet/in.h stddef.h stdlib.h string.h strings.h sys/file.h sys/socket.h sys/time.h unistd.h])

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

AC_ARG_ENABLE(mode,AC_HELP_STRING(--enable-mode,[specify the compilation mode, either production, debug or profiling]),,enable_mode=production)
AC_DEFINE_UNQUOTED(COMPILE_MODE,"$enable_mode",[Compilation mode])

ES_CHECK_MPI
ES_CHECK_EFENCE
ES_CHECK_INLINING
ES_CHECK_LINK
ES_CHECK_TCL
ES_CHECK_TK
ES_CHECK_FFTW
AC_CHECK_PROG(doxygen, doxygen, yes, no)
ES_CHECK_DEPEND

AC_DEFINE_UNQUOTED(LAST_CHANGE,"$time_stamp",[When the last change was done])

# Creating the output
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

