dnl -*- mode: autoconf -*-
dnl This file contains macros that test for different features of the
dnl C-compiler and linker.

AC_DEFUN([ES_CHECK_COMPILER],[
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
	AC_ARG_ENABLE(debug,
		AC_HELP_STRING(--enable-debug,
			[Turn on debugging]),
		,enable_debug=no)
	if test .$enable_debug = .yes; then
	   AC_DEFINE(DEBUG,,[Turn on debugging?])
	fi

	# test whether to enable profiling
	AC_ARG_ENABLE(profiling,
		AC_HELP_STRING(--enable-profiling,
			[Turn on profiling]),
		,enable_profiling=no)
	if test .$enable_profiling = .yes; then
	   AC_DEFINE(PROFILING,,[Turn on profiling?])
	fi

	ES_CHECK_INLINING

	# test whether to enable processor optimisation
	AC_ARG_ENABLE(processor-optimization,
		AC_HELP_STRING(--enable-processor-optimization,
		[enable guessing of processor specific optimizations])
		,,enable_processor_optimization=yes)
	if test .$enable_processor_optimization = .yes; then
	        ES_CHECK_OPTFLAGS
	fi

	ES_CHECK_LINK
])


dnl ********************************** inlining ******************************
AC_DEFUN([ES_CHECK_INLINING],[
	if test .$enable_debug = .yes || test .$enable_profiling = .yes; then
		mdinline=static
	else
		AC_MSG_CHECKING([how to inline functions])
		for mdinline in "static inline" "inline static" "inline" "static"; do
		    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
/**** This HAS to be at the beginning of the line for some compilers ****/
#define MDINLINE $mdinline
MDINLINE void test() {}],
			[test();])],[works=yes],[works=no])
		    if test .$works = .yes; then break; fi
		done
		if test .$works = .no; then
		   AC_MSG_ERROR([your compiler does not even support "static"])
		fi
		AC_MSG_RESULT([$mdinline])
	fi
	AC_DEFINE_UNQUOTED(MDINLINE,$mdinline,[How to inline functions])
])

dnl ***************************** compiler type ******************************

AC_DEFUN([ES_CHECK_DEC],[
	if $CC -V 2>&1 | grep "Compaq C" >/dev/null 2>&1; then
		with_compilertype=dec
	fi
	if $CC -V 2>&1 | grep "Digital UNIX Compiler" >/dev/null 2>&1; then
		with_compilertype=dec
	fi
])

AC_DEFUN([ES_CHECK_ICC],[
	if $CC -V 2>&1 | grep "Intel" >/dev/null 2>&1; then
		with_compilertype=icc
	fi
])

AC_DEFUN([ES_CHECK_GCC],[
	if $CC -v 2>&1 | grep "gcc" >/dev/null 2>&1; then
		with_compilertype=gcc
	fi
])

AC_DEFUN([ES_CHECK_XLC],[
	if $CC 2>&1 | grep "xlc" >/dev/null 2>&1; then
		with_compilertype=xlc
	fi
])

AC_DEFUN([ES_CHECK_COMPILERTYPE],[
	AC_ARG_WITH(compilertype,
                AC_HELP_STRING(--with-compilertype=TYPE,
                        [type of the compiler]),
                [ test .$with_compilertype = .no && with_compilertype=unknown ],
                [ with_compilertype=unknown ])
	AC_MSG_CHECKING([the compiler type])
	ES_CHECK_XLC
	if test .$with_compilertype = .unknown; then
		ES_CHECK_DEC
		if test .$with_compilertype = .unknown; then
			ES_CHECK_ICC
			if test .$with_compilertype = .unknown; then
				ES_CHECK_GCC
			fi
		fi
	fi
	AC_MSG_RESULT($with_compilertype)
])

dnl ******************************* optimizations ********************************

AC_DEFUN([ES_TRY_ADD_CFLAG],[
	AC_MSG_CHECKING([whether the compiler accepts $1])
	saved_CFLAGS=$CFLAGS
	CFLAGS="$1 $CFLAGS"
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],[AC_MSG_RESULT(yes);try_add_flag_res=yes],[AC_MSG_RESULT(no); CFLAGS=$saved_CFLAGS; try_add_flag_res=no ])
	if test .$2 != . ; then
		$2=$try_add_flag_res
	fi
])

AC_DEFUN([ES_TRYLINK_ADD_CFLAG],[
	AC_MSG_CHECKING([whether the compiler accepts $1])
	saved_CFLAGS=$CFLAGS
	CFLAGS="$1 $CFLAGS"
	AC_LINK_IFELSE([AC_LANG_PROGRAM([],[])],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no); CFLAGS=$saved_CFLAGS ])
])

AC_DEFUN([ES_CHECK_GCC_FLAGS],[
	case $target_cpu in
	i*86)		type=ia; m=32; arch=$target_cpu ;;
	Core*)		type=ia; m=64; arch=core2; sse=3 ;;
	Pentium_III)	type=ia; m=32; arch=pentium3; sse=1 ;;
	Atom)		type=ia; m=32; arch=pentiumpro; sse=2 ;;
	Celeron)	type=ia; m=32; arch=pentium4; sse=2 ;;
	Pentium_4)	type=ia; m=32; arch=pentium4; sse=2 ;;
	Pentium_M)	type=ia; m=32; arch=pentium-m; sse=2 ;;
	Xeon)		type=ia; m=32; arch=pentium4; sse=2 ;;
	Pentium_4_64)	type=ia; m=64; arch=nocona; sse=2 ;;
	Xeon_64)	type=ia; m=64; arch=nocona; sse=2 ;;
	Athlon)		type=ia; m=32; arch=athlon ;;
	Athlon_MP)	type=ia; m=32; arch=athlon-mp; sse=1 ;;
	Athlon_XP)	type=ia; m=32; arch=athlon-xp; sse=1 ;;
	Opteron)	type=ia; m=64; arch=opteron; sse=2 ;;
	Phenom)		type=ia; m=64; arch=barcelona; sse=3 ;;
	Athlon_64)	type=ia; m=64; arch=athlon64; sse=2;;
	x86_64)		type=ia; m=64; arch=k8 ;;
	Power)		type=pwr; cpu=power ;;
	Power2)		type=pwr; cpu=power2 ;;
	Power*)		type=pwr; cpu=powerpc ;;
	ppc7450)	type=pwr; cpu=7450; altivec=yes ;;
	ppc*)		type=pwr; cpu=powerpc; altivec=yes ;;
	EV67*)		type=alpha; cpu=ev67 ;;
	EV6*)		type=alpha; cpu=ev6 ;;
	EV56*)		type=alpha; cpu=ev56 ;;
	EV5*)		type=alpha; cpu=ev5 ;;
	EV4*)		type=alpha; cpu=ev4 ;;
	sparcv7)	type=sparc; cpu=v7 ;;
	sparcv8)	type=sparc; cpu=v8 ;;
	sparcv9)	type=sparc; cpu=v9 ;;
	sparclite)	type=sparc; cpu=sparclite ;;
	sparclet)	type=sparc; cpu=sparclet ;;
	*)	AC_MSG_WARN([could not recognize your cpu type, relying on generic optimization])
	esac

	# first see if the compiler can autotune

	# try whether the compiler accepts -fastsse (PGI)
	ES_TRY_ADD_CFLAG(-fastsse, accept_fastsse)
	if test .$accept_fastsse = .yes; then
		accept_fast=yes;
	else
		# try just -fast (Apple/PGI)
		ES_TRY_ADD_CFLAG(-fast, accept_fast)
	fi

	# on a 64bit machine, set a flag to test for
	if test .$m = .64; then
	   	ES_64BIT=1;
	fi

	# try CPU specific flags
	if test .$m != .; then
		ES_TRY_ADD_CFLAG(-m$m)
	fi
	if test .$arch != .; then
		ES_TRY_ADD_CFLAG(-march=$arch)
	fi
	if test .$cpu != .; then
		ES_TRY_ADD_CFLAG(-mcpu=$cpu)
	fi

	# if fast works, then do not try CPU specific manual optimization flags
	if test .$accept_fast != .yes ; then

	if test .$altivec != .; then
		ES_TRY_ADD_CFLAG(-maltivec)
	fi
	if test .$sse = .1; then
		ES_TRY_ADD_CFLAG([-msse -mfpmath=sse])
	else
		if test .$sse != .; then
			ES_TRY_ADD_CFLAG([-msse$sse -mfpmath=sse])
		fi
	fi
	if test .$type = .ia; then
		ES_TRY_ADD_CFLAG(-malign-double)
	fi

	fi
	# end of CPU specific flags

	ES_TRY_ADD_CFLAG(-O3)
	ES_TRY_ADD_CFLAG(-Wno-unused-function)
	ES_TRY_ADD_CFLAG(-Wall)
	ES_TRY_ADD_CFLAG(-ffast-math)
	ES_TRY_ADD_CFLAG(-floop-optimize)
	ES_TRY_ADD_CFLAG(-funroll-loops)

	if test .$enable_debug = .yes; then
		ES_TRY_ADD_CFLAG(-g)
	else
		ES_TRY_ADD_CFLAG(-fomit-frame-pointer)
	fi
	if test .$enable_profiling = .yes; then
		ES_TRY_ADD_CFLAG(-pg)
		LDFLAGS="-pg $LDFLAGS"
	fi	  
	if test .$enable_debug = .yes || test .$enable_profiling = .yes; then
		ES_TRY_ADD_CFLAG(-fno-inline)
	else
		ES_TRY_ADD_CFLAG(-finline-limit=1000000)
	fi
])

AC_DEFUN([ES_CHECK_ICC_FLAGS],[
	# check for svml for the vectorizer
	AC_MSG_CHECKING([whether we have libsvml])
	save_LIBS=$LIBS
	svml_found=no
	LIBS="$LIBS -lsvml"
	# we just try to link, because the functions are not well documented
	# and can change
	AC_LINK_IFELSE([AC_LANG_SOURCE(int main() {})],[svml_found=yes],[])
	if test .$svml_found = .yes; then
		AC_MSG_RESULT(yes)
	else
		AC_MSG_RESULT(no)
		LIBS=$save_LIBS
	fi

	case $target_cpu in
	Pentium_III)	xflag=K ;;
	Pentium_4)	xflag=N ;;
	Pentium_M)	xflag=B ;;
	Xeon)		xflag=N ;;
	Athlon)		xflag=K ;;
	Athlon_MP)	xflag=K ;;
	Athlon_XP)	xflag=K ;;
	Opteron)	xflag=W ;;
	Athlon_64)	xflag=W ;;
	Itanium)	;;
	*)	AC_MSG_WARN([could not recognize your cpu type, relying on generic optimization])
	esac
	if test .$xflag != .; then
		ES_TRY_ADD_CFLAG(-x$xflag)
	fi
	ES_TRY_ADD_CFLAG(-O3)
	ES_TRY_ADD_CFLAG(-ip)
	ES_TRY_ADD_CFLAG(-unroll)
	ES_TRY_ADD_CFLAG(-w1)
	dnl stupid warnings from version 8.1 and above
	ES_TRY_ADD_CFLAG([-wd 1572,1173])

	if test .$enable_debug = .yes; then
		ES_TRY_ADD_CFLAG(-g)
	fi
	if test .$enable_profiling = .yes; then
		ES_TRY_ADD_CFLAG(-pg)
		LDFLAGS="-pg $LDFLAGS"
	fi
	if test .$enable_debug = .yes || test .$enable_profiling = .yes; then
		ES_TRY_ADD_CFLAG(-Ob0)
	fi
])

AC_DEFUN([ES_CHECK_XLC_FLAGS],[
	case $target_cpu in
	Power)		cpu=pwr ;;
	Power2)		cpu=pwr2 ;;
	Power3)		cpu=pwr3 ;;
	Power4)		cpu=pwr4 ;;
	Power5)		cpu=pwr5 ;;
	ppc*)		cpu=auto ;;
	*)	AC_MSG_WARN([could not recognize your cpu type, relying on generic optimization])
	esac
	if test .$cpu != .; then
		ES_TRY_ADD_CFLAG(-qarch=$cpu)
		ES_TRY_ADD_CFLAG(-qtune=$cpu)
	fi
	ES_TRY_ADD_CFLAG(-O3)
	ES_TRY_ADD_CFLAG([-Wc,-qipa])
	ES_TRY_ADD_CFLAG(-qhot)
	ES_TRY_ADD_CFLAG(-qfloat=rsqrt:fltint:fold:maf)
	ES_TRY_ADD_CFLAG(-qcpluscmt)
	ES_TRY_ADD_CFLAG(-qmaxmem=-1)

	AC_ARG_ENABLE(xlc-qipa, AC_HELP_STRING(--enable-xlc-qipa,[enable the IPA of the xlc at linking]),,enable_xlc_qipa=yes)
	if test .$enable_xlc_qipa = .yes; then
		LDFLAGS="-qipa $LDFLAGS"
		AC_MSG_WARN([****** WARNING: -qipa is enabled by default ******])
		AC_MSG_WARN([  if you encounter linking problems, such as *.dylib])
		AC_MSG_WARN([  not found, disable it with --disable-xlc-qipa])
	fi

	if test .$enable_debug = .yes; then
		ES_TRY_ADD_CFLAG(-g)
	fi
	if test .$enable_profiling = .yes; then
		ES_TRY_ADD_CFLAG(-pg)
		LDFLAGS="-pg $LDFLAGS"
	fi	  
	if test .$enable_debug! = .yes && test .$enable_profiling=!.yes; then
		ES_TRY_ADD_CFLAG(-qinline=100000)
	fi
])

AC_DEFUN([ES_CHECK_DEC_FLAGS],[
	dnl ok, not good when cross-compiling, but simple...
	ES_TRY_ADD_CFLAG(-tune host)
	ES_TRY_ADD_CFLAG(-fp_reorder)
	ES_TRY_ADD_CFLAG(-std)
	dnl the dec compiler has LOTS of silly warnings
	for flag in uncalled ignorecallval noparmlist uncalled alignment unknownmacro \
        	strctpadding unusedincl needconstext newlocale cxxcomment \
		unusedtop unnecincl nestincl uninit1; do
		ES_TRY_ADD_CFLAG([-msg_disable $flag])
	done
	dnl -w* flags are either for warnings or for the linker
	dnl in the latter case, it compiles nicely, but does not link...
	ES_TRYLINK_ADD_CFLAG(-w0 -check)

	if test .$enable_debug = .yes; then
		ES_TRY_ADD_CFLAG(-g)
	fi
	if test .$enable_profiling = .yes; then
		ES_TRY_ADD_CFLAG(-pg)
		LDFLAGS="-pg $LDFLAGS"
	fi
	if test .$enable_debug! = .yes && test .$enable_profiling! = .yes; then
		ES_TRY_ADD_CFLAG(-O3)
	fi
])

AC_DEFUN([ES_CHECK_OPTFLAGS],[
	ES_CHECK_COMPILERTYPE
	case $with_compilertype in
	gcc)	ES_CHECK_GCC_FLAGS ;;
	icc)	ES_CHECK_ICC_FLAGS ;;
	xlc)	ES_CHECK_XLC_FLAGS ;;
	dec)	ES_CHECK_DEC_FLAGS ;;
	unknown) 
	if test .$enable_debug = .yes; then
		ES_TRY_ADD_CFLAG(-g)
	fi
	if test .$enable_profiling = .yes; then
		ES_TRY_ADD_CFLAG(-pg)
		LDFLAGS="-pg $LDFLAGS"
	fi
	if test .$enable_debug! = .yes && test .$enable_profiling! = .yes; then
		ES_TRY_ADD_CFLAG(-O2)
	fi 
	;;
	esac
])

AC_DEFUN([ES_SET_ARCH_FLAGS],[
	dnl os/binutils dependent stuff
	case $target_os in
	*aix*)	if test .$with_compilertype = .xlc; then
			LDFLAGS="$LDFLAGS -brtl"
		else
			LDFLAGS="$LDFLAGS -Wl,-brtl"
		fi
	esac
])

dnl ******************************* linking ********************************
AC_DEFUN([ES_CHECK_LINK],[
	case $target_os in
	*darwin*) LD="`pwd`/darwinlink.sh $CC" ;;
	*) LD=$CC ;;
	esac
])

AC_DEFUN([ES_NOTE_64BIT],[
if test "$ES_64BIT"; then
   AC_MSG_NOTICE([
********************************************************************************
* Note, that you are using a 64bit machine and ESPResSo tries to build a 64bit *
* binary. Please make sure that you are using a 64bit version of the library!  *
********************************************************************************])
fi
])
