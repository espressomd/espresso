dnl ********************************** inlining ******************************
AC_DEFUN([CF_CHECK_INLINE],[
	AC_MSG_CHECKING([how to inline functions])
	for mdinline in "static inline" "inline static" "inline" "static"; do
		AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
/**** This HAS to be at the beginning of the line for some compilers ****/
#define MDINLINE $mdinline
MDINLINE void test() {}],[test();])],[works=yes],[works=no])
		if test .$works = .yes; then break; fi
	done
	if test .$works = .no; then
		AC_MSG_ERROR([your compiler does not even support "static"])
	fi
	AC_DEFINE_UNQUOTED(MDINLINE,$mdinline)
	AC_MSG_RESULT([$mdinline])
])

dnl ***************************** compiler type ******************************

AC_DEFUN([CF_CHECK_DEC],[
	if $CC -V 2>&1 | grep "Compaq C" >/dev/null 2>&1; then
		compiler_type=dec
	fi
	if $CC -V 2>&1 | grep "Digital UNIX Compiler" >/dev/null 2>&1; then
		compiler_type=dec
	fi
])

AC_DEFUN([CF_CHECK_ICC],[
	if $CC -V 2>&1 | grep "Intel" >/dev/null 2>&1; then
		compiler_type=icc
	fi
])

AC_DEFUN([CF_CHECK_GCC],[
	if $CC -v 2>&1 | grep "gcc" >/dev/null 2>&1; then
		compiler_type=gcc
	fi
])

AC_DEFUN([CF_CHECK_XLC],[
	if $CC 2>&1 | grep "xlc" >/dev/null 2>&1; then
		compiler_type=xlc
	fi
])

AC_DEFUN([CF_CHECK_COMPILERTYPE],[
	compiler_type=unknown
	AC_MSG_CHECKING([the compiler type])
	CF_CHECK_XLC
	if test $compiler_type = unknown; then
		CF_CHECK_DEC
		if test $compiler_type = unknown; then
			CF_CHECK_ICC
			if test $compiler_type = unknown; then
				CF_CHECK_GCC
			fi
		fi
	fi
	AC_MSG_RESULT($compiler_type)
])

dnl ******************************* optimizations ********************************

AC_DEFUN([CF_TRY_ADD_CFLAG],[
	AC_MSG_CHECKING([whether the compiler accepts $1])
	saved_CFLAGS=$CFLAGS
	CFLAGS="$1 $CFLAGS"
	AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],[])],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no); CFLAGS=$saved_CFLAGS ])
])

AC_DEFUN([CF_TRYLINK_ADD_CFLAG],[
	AC_MSG_CHECKING([whether the compiler accepts $1])
	saved_CFLAGS=$CFLAGS
	CFLAGS="$1 $CFLAGS"
	AC_LINK_IFELSE([AC_LANG_PROGRAM([],[])],[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no); CFLAGS=$saved_CFLAGS ])
])

AC_DEFUN([CF_CHECK_GCC_FLAGS],[
	case $target_cpu in
	Pentium_III)	type=ia; m=32; arch=pentium3; sse=yes ;;
	Pentium_4)	type=ia; m=32; arch=pentium4; sse2=yes ;;
	Pentium_M)	type=ia; m=32; arch=pentium-m; sse2=yes ;;
	Xeon)		type=ia; m=32; arch=pentium4; sse2=yes ;;
	Athlon)		type=ia; m=32; arch=athlon ;;
	Athlon_MP)	type=ia; m=32; arch=athlon-mp; sse=yes ;;
	Athlon_XP)	type=ia; m=32; arch=athlon-xp; sse=yes ;;
	Opteron)	type=ia; m=64; arch=opteron; sse2=yes ;;
	Athlon_64)	type=ia; m=64; arch=athlon64; sse2=yes;;
	Power)		type=pwr; cpu=power ;;
	Power2)		type=pwr; cpu=power2 ;;
	Power*)		type=pwr; cpu=powerpc ;;
	ppc*)		case $target_os in
			*darwin*)	type=pwr
					case $target_cpu in
					ppc970)	 cpu=G5 ;;
					ppc7450) cpu=7450 ;;
					ppc*)    cpu=powerpc ;;
					esac
					altivec=yes ;;
			*)		type=pwr; cpu=powerpc;
			esac ;;
	EV67*)		type=alpha; cpu=ev67 ;;
	EV6*)		type=alpha; cpu=ev6 ;;
	EV56*)		type=alpha; cpu=ev56 ;;
	EV5*)		type=alpha; cpu=ev5 ;;
	EV4*)		type=alpha; cpu=ev4 ;;
	*)	AC_MSG_WARN([could not recognize your cpu type, relying on generic optimization])
	esac
	if test .$m != .; then
		CF_TRY_ADD_CFLAG(-m$m)
	fi
	if test .$arch != .; then
		CF_TRY_ADD_CFLAG(-march=$arch)
	fi
	if test .$cpu != .; then
		CF_TRY_ADD_CFLAG(-mcpu=$cpu)
	fi
	if test .$altivec != .; then
		CF_TRY_ADD_CFLAG(-maltivec)
	fi
	if test .$sse != .; then
		CF_TRY_ADD_CFLAG([-msse -mfpmath=sse])
	fi
	if test .$sse2 != .; then
		CF_TRY_ADD_CFLAG([-msse2 -mfpmath=sse])
	fi
	if test .$type = .ia; then
		CF_TRY_ADD_CFLAG(-malign-double)
	fi
	CF_TRY_ADD_CFLAG(-O3)
	CF_TRY_ADD_CFLAG(-Wall)
	CF_TRY_ADD_CFLAG(-ffast-math)
	CF_TRY_ADD_CFLAG(-floop-optimize)
	CF_TRY_ADD_CFLAG(-funroll-loops)

	case $enable_mode in
	production) CF_TRY_ADD_CFLAG(-finline-limit=1000000)
		CF_TRY_ADD_CFLAG(-fomit-frame-pointer)
		;;
	debug)	CF_TRY_ADD_CFLAG(-g)
		CF_TRY_ADD_CFLAG(-fno-inline)
		;;
	profile) CF_TRY_ADD_CFLAG(-pg)
		CF_TRY_ADD_CFLAG(-fno-inline)
		CF_TRY_ADD_CFLAG(-fomit-frame-pointer)
		LDFLAGS="-pg $LDFLAGS"
		;;
	esac
])

AC_DEFUN([CF_CHECK_ICC_FLAGS],[
	LIBS="$LIBS -lsvml"

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
	*)	AC_MSG_WARN([could not recognize your cpu type, relying on generic optimization])
	esac
	if test .$xflag != .; then
		CF_TRY_ADD_CFLAG(-x$xflag)
	fi
	CF_TRY_ADD_CFLAG(-O3)
	CF_TRY_ADD_CFLAG(-ip)
	CF_TRY_ADD_CFLAG(-unroll)
	CF_TRY_ADD_CFLAG(-w1)
	dnl stupid warnings from version 8.1 and above
	CF_TRY_ADD_CFLAG([-wd 1572,1173])

	case $enable_mode in
	production)
		;;
	debug)	CF_TRY_ADD_CFLAG(-g)
		CF_TRY_ADD_CFLAG(-Ob0)
		;;
	profile) CF_TRY_ADD_CFLAG(-pg)
		CF_TRY_ADD_CFLAG(-Ob0)
		LDFLAGS="-pg $LDFLAGS"
		;;
	esac
])

AC_DEFUN([CF_CHECK_XLC_FLAGS],[
	case $target_cpu in
	Power)		cpu=pwr ;;
	Power2)		cpu=pwr2 ;;
	Power3)		cpu=pwr3 ;;
	Power4)		cpu=pwr4 ;;
	ppc*)		cpu=auto ;;
	*)	AC_MSG_WARN([could not recognize your cpu type, relying on generic optimization])
	esac
	if test .$cpu != .; then
		CF_TRY_ADD_CFLAG(-qarch=$cpu)
		CF_TRY_ADD_CFLAG(-qtune=$cpu)
	fi
	CF_TRY_ADD_CFLAG(-O3)
	CF_TRY_ADD_CFLAG(-qipa)
	CF_TRY_ADD_CFLAG(-qhot)
	CF_TRY_ADD_CFLAG(-qfloat=rsqrt:fltint:fold:maf)
	CF_TRY_ADD_CFLAG(-qcpluscmt)
	CF_TRY_ADD_CFLAG(-qmaxmem=-1)

	AC_ARG_ENABLE(xlc-qipa, AC_HELP_STRING(--enable-xlc-qipa,[enable the IPA of the xlc at linking]),,enable_xlc_qipa=yes)
	if test .$enable_xlc_qipa = .yes; then
		LDFLAGS="-qipa $LDFLAGS"
		AC_MSG_WARN([****** WARNING: -qipa is enabled by default ******])
		AC_MSG_WARN([  if you encounter linking problems, such as *.dylib])
		AC_MSG_WARN([  not found, disable it with --disable-xlc-qipa])
	fi

	case $enable_mode in
	production) CF_TRY_ADD_CFLAG(-qinline=100000)
		;;
	debug)	CF_TRY_ADD_CFLAG(-g)
		;;
	profile) CF_TRY_ADD_CFLAG(-pg)
		LDFLAGS="-pg $LDFLAGS"
		;;
	esac
])

AC_DEFUN([CF_CHECK_DEC_FLAGS],[
	dnl ok, not good when cross-compiling, but simple...
	CF_TRY_ADD_CFLAG(-tune host)
	CF_TRY_ADD_CFLAG(-fp_reorder)
	CF_TRY_ADD_CFLAG(-std)
	dnl the dec compiler has LOTS of silly warnings
	for flag in uncalled ignorecallval noparmlist uncalled alignment unknownmacro \
        	strctpadding unusedincl needconstext newlocale cxxcomment \
		unusedtop unnecincl nestincl uninit1; do
		CF_TRY_ADD_CFLAG([-msg_disable $flag])
	done
	dnl -w* flags are either for warnings or for the linker
	dnl in the latter case, it compiles nicely, but does not link...
	CF_TRYLINK_ADD_CFLAG(-w0 -check)

	case $enable_mode in
	production) CF_TRY_ADD_CFLAG(-O3)
		;;
	debug)	CF_TRY_ADD_CFLAG(-g)
		;;
	profile) CF_TRY_ADD_CFLAG(-pg)
		LDFLAGS="-pg $LDFLAGS"
		;;
	esac
])

AC_DEFUN([CF_CHECK_OPTFLAGS],[
	CF_CHECK_COMPILERTYPE
	case $compiler_type in
	gcc)	CF_CHECK_GCC_FLAGS ;;
	icc)	CF_CHECK_ICC_FLAGS ;;
	xlc)	CF_CHECK_XLC_FLAGS ;;
	dec)	CF_CHECK_DEC_FLAGS ;;
	unknown) case $enable_mode in
		production)	CFLAGS="-O2 $CFLAGS" ;;
		debug)		CFLAGS="-g -O2 $CFLAGS" ;;
		debug)		CFLAGS="-pg-O2 $CFLAGS"
		 		LDFLAGS="-pg $LDFLAGS" ;;
		esac
	esac
])

AC_DEFUN([CF_SET_ARCH_FLAGS],[
	dnl os/binutils dependent stuff
	case $target_os in
	*aix*)	if test .$compiler_type = .xlc; then
			LDFLAGS="$LDFLAGS -brtl"
		else
			LDFLAGS="$LDFLAGS -Wl,-brtl"
		fi
	esac
])
