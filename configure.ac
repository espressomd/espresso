#######################################################################
# Used by autogen.sh to produce configure.ac for autoconf             #
#######################################################################

AC_PREREQ(2.50)
AC_INIT(ESPResSo, [v1.9.9h (Seska)], espresso@mpip-mainz.mpg.de)
time_stamp="May 10th, 2006"

if test "$CC" != ""; then	user_defined_CC=yes; fi
if test "$CFLAGS" != ""; then	user_defined_CFLAGS=yes; fi
if test "$LDFLAGS" != ""; then	user_defined_LDFLAGS=yes; fi
if test "$LIBS" != ""; then	user_defined_LIBS=yes; fi

AC_CONFIG_AUX_DIR(config)
AC_CANONICAL_HOST
AC_CANONICAL_TARGET
AC_CONFIG_SRCDIR(main.c)
AC_PREFIX_DEFAULT($HOME/Espresso)

target=$target_cpu-$target_vendor-$target_os

AC_ARG_VAR(LIBS,[use this to specify additional libraries to link against, e.g. -lmx])

AC_ARG_ENABLE(mode,AC_HELP_STRING(--enable-mode,[specify the compilation mode, either production, debug or profiling]),,enable_mode=production)
AC_ARG_ENABLE(processor-optimization,AC_HELP_STRING(--enable-processor-optimization,[enable guessing of processor specific optimizations]),,enable_processor_optimization=yes)

dnl remove the default (optimization) flags
if test .$user_defined_CFLAGS != .yes; then CFLAGS=""; fi
if test .$user_defined_LDFLAGS != .yes; then LDFLAGS=""; fi
if test .$user_defined_LIBS != .yes; then LIBS=""; fi
dnl then load any known configs
ES_KNOWN_CONFIGS

LIBS="-lm $LIBS"

ES_CHECK_MPI
ES_CHECK_EFENCE
ES_CHECK_INLINING
if test .$enable_processor_optimization = .yes; then
	CF_CHECK_OPTFLAGS
fi
CF_SET_ARCH_FLAGS
ES_CHECK_LINK
ES_CHECK_TCL
ES_CHECK_TK
ES_CHECK_FFTW
AC_CHECK_PROG(doxygen, doxygen, yes, no)
ES_CHECK_DEPEND

CFLAGS="$CPPFLAGS -Iobj-$target $CFLAGS"
AC_DEFINE_UNQUOTED(LAST_CHANGE,"$time_stamp")
AC_DEFINE_UNQUOTED(COMPILE_MODE,"$enable_mode")
AC_SUBST(DEPEND)
AC_SUBST(LD)
AC_SUBST(ADD_SOURCES)
AC_SUBST(MPI_INVOCATION)
# check for socket
AC_MSG_CHECKING([whether socket is available without -lsocket])
AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY(socket)],
    [AC_MSG_RESULT(yes)],[
	AC_MSG_RESULT(no)
	AC_CHECK_LIB(socket,socket,,AC_MSG_ERROR([could not find location of socket().]))
	LDFLAGS="$LDFLAGS -lsocket"
    ])
# check for gethostbyname
AC_MSG_CHECKING([whether gethostbyname is available without -lnsl])
AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY(gethostbyname)],
    [AC_MSG_RESULT(yes)],[
	AC_MSG_RESULT(no)
	AC_CHECK_LIB(nsl,gethostbyname,,AC_MSG_ERROR([could not find location of gethostbyname().]))
	LDFLAGS="$LDFLAGS -lnsl"
    ])


AC_CONFIG_HEADERS([obj-$target/archconfig.h:archconfig.h.in])
AC_CONFIG_FILES([Makefile.$target:Makefile.in obj-$target/Espresso:config/Espresso.in])
AC_CONFIG_COMMANDS([starter],[
	sed -e 's/@ARGUMENTS@/\$\*/g' -e s/@NP@/\$NP/g obj-$target/Espresso > obj-$target/Espresso.new
	mv obj-$target/Espresso.new obj-$target/Espresso
	chmod a+x obj-$target/Espresso],[target=$target])
AC_OUTPUT

echo "****************************************************************"
echo "*                   configuration summary                      *"
echo "****************************************************************"
echo ""
echo "platform detected = $target"
echo ""
echo "compiler         = $CC"
echo "linker           = $LD"
echo "c compiler flags = $CFLAGS"
echo "linker flags     = $LDFLAGS"
echo "libraries        = $LIBS"
echo "compilation mode = $enable_mode"
echo "MPI library      = $enable_mpi"
if test .$enable_mpi = .generic; then
echo "    WARNING: a generic MPI version may or may not work, look at the Espresso script in obj-$target"
fi
echo "Tcl library      = $tclversion"
if test $tclversion != tcl8.5 && test $tclversion != tcl8.4; then
echo "       WARNING: Not all supplied scripts work properly with Tcl version before 8.4, consider upgrading"
fi 
if test $tkversion != no; then
echo "Tk library       = $tkversion"
else
echo "Tk disabled (enable with --enable-tk)"
fi
echo "FFTW version     = $enable_fftw"
if test .$doxygen = .no; then
echo "To be able to build the documentation, add doxygen to your PATH."
fi
