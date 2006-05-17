m4_include([config/knownconfigs.m4])
m4_include([config/mpi.m4])
m4_include([config/compilerflags.m4])

dnl search for a library in additional paths
AC_DEFUN([ES_ADDPATH_CHECK_LIB],[
	AC_MSG_CHECKING([for lib$1])
	save_LDFLAGS=$LDFLAGS
	save_LIBS=$LIBS
	adp_found=no
	dnl let's see whether it's in the default paths
	LIBS="-l$1 $save_LIBS"
	AC_LINK_IFELSE([AC_LANG_CALL([],[$2])],[adp_found=yes],[])

	if test .$adp_found = .no; then
		for path in $5 /sw/lib /usr/lib64 /usr/local/lib64 /opt/lib64 /usr/lib /usr/local/lib /opt/lib; do
			LDFLAGS="$save_LDFLAGS -L$path"
			AC_LINK_IFELSE([AC_LANG_CALL([],[$2])],[adp_found=yes],[])
			if test .$adp_found = .yes; then break; fi
		done
	fi
	if test .$adp_found = .yes; then
		AC_MSG_RESULT(yes)
	else
		AC_MSG_RESULT(no)
		LDFLAGS=$save_LDFLAGS
		LIBS=$save_LIBS
	fi
	AS_IF([test .$adp_found = .yes], [$3],[$4])
])

dnl search for a header file in additional paths
AC_DEFUN([ES_ADDPATH_CHECK_HEADER],[
	AC_MSG_CHECKING([for $1])
	save_CFLAGS=$CFLAGS
	adp_found=no
	dnl let's see whether it's in the default paths
	AC_COMPILE_IFELSE([
		#include <$1>
	],[adp_found=yes],[])

	if test .$adp_found = .no; then
		for path in $4 /sw/include /usr/include /usr/local/include /opt/include; do
			CFLAGS="$save_CFLAGS -I$path"
			AC_COMPILE_IFELSE([
				#include <$1>
			],[adp_found=yes],[])
			if test .$adp_found = .yes; then break; fi
		done
	fi
	if test .$adp_found = .yes; then
		AC_MSG_RESULT(yes)
	else
		AC_MSG_RESULT(no)
		CFLAGS=$save_CFLAGS
	fi
	AS_IF([test .$adp_found = .yes], [$2],[$3])
])

AC_DEFUN([ES_CHECK_LINK],[
	case $target_os in
	*darwin*) LD="`pwd`/darwinlink.sh $CC" ;;
	*) LD=$CC ;;
	esac
])

AC_DEFUN([ES_CHECK_EFENCE],[
	AC_ARG_ENABLE(efence,AC_HELP_STRING(--enable-efence,[use ElectricFence memory debugging for the debug binary]),,enable_efence=no)
	if test $enable_efence = yes; then
		AC_CHECK_LIB(efence,malloc,,AC_MSG_ERROR([could not link against the efence library]))
		CFLAGS="$CFLAGS -DEFENCE"
		LDFLAGS="$LDFLAGS -lefence"
	fi
])

AC_DEFUN([ES_CHECK_MPI],[
	AC_ARG_ENABLE(mpi,AC_HELP_STRING(--enable-mpi=<type>,[compile with MPI (parallelization) support
	if you specify <type>=fake,generic,lam,mpich,poe or dmpi,
	the corresponding MPI implementation will used,
	otherwise the native environment of the platform is used. If
	none is found, the fake implementation for only one processor
	is used]),,enable_mpi=yes)
	if test $enable_mpi = fake; then
		MPI_SETUP_FAKE
	elif test $enable_mpi = no; then
		enable_mpi=fake
		MPI_SETUP_FAKE
	else
		if test $enable_mpi = yes; then
			if test .$use_mpi != .; then
				enable_mpi=$use_mpi
			fi
		fi
		if test $enable_mpi = yes; then
			MPI_GUESS_ENV
		else
			MPI_SETUP($enable_mpi)
		fi
	fi
])

AC_DEFUN([ES_CHECK_DEPEND],[
	AC_CHECK_PROGS(DEPEND, [makedepend mkdep])
	if test .$DEPEND = .; then
		if test .$compiler_type = .gcc; then
			DEPEND=$CC
		fi
	fi
])

AC_DEFUN([ES_CHECK_INLINING],[
	if test $enable_mode = production; then
		CF_CHECK_INLINE
	else
		CFLAGS="$CFLAGS -DMDINLINE=\"static\""
	fi
])

AC_DEFUN([ES_CHECK_TCL],[
	AC_DEFINE(USE_NON_CONST)

	AC_ARG_ENABLE(tcl,AC_HELP_STRING([--enable-tcl],[specify the tcl library to use (e.g. tcl8.4)]),
		[tclversion=$enable_tcl], [tclversion=yes])
	if test .$tclversion = .yes; then
		tclversion=""
		for version in $TCL_VERSION tcl8.5 tcl8.4 tcl8.3 tcl8.2 tcl; do
			ES_ADDPATH_CHECK_LIB($version, Tcl_Init, [tclversion=$version], [])
			if test .$tclversion != .; then break; fi
		done
	else
		ES_ADDPATH_CHECK_LIB($tclversion, Tcl_Init, [], [tclversion=""])
	fi
	if test .$tclversion = .; then
		case $target_os in
		*darwin*) AC_MSG_NOTICE([If you have Tcl installed, make sure that in one of the library paths, e.g. /usr/local/lib,
	there is a link from lib<tclversion>.dylib to the Tcl library, which usually is
	/Library/Frameworks/Tcl.framework/Tcl.]) ;;
		*)	AC_MSG_ERROR([Tcl library $tclversion not found]) ;;
		esac
	fi
	case $target_os in
	*darwin*) extrapaths=/Library/Frameworks/Tcl.framework/Headers ;;
	*) ;;
	esac
	ES_ADDPATH_CHECK_HEADER(tcl.h, [], [AC_MSG_ERROR(TCL headers not found)],$extrapaths)
])

AC_DEFUN([ES_CHECK_TK],[
	AC_ARG_ENABLE(tk,AC_HELP_STRING([--enable-tk],[tk version for the GUI to use, yes for auto check, no to disable]),
		[tkversion=$enable_tk], [tkversion=no])
	if test .$tkversion != .no; then
		dnl test for X11
		AC_PATH_XTRA
		saved_CFLAGS=$CFLAGS
		saved_LDFLAGS=$LDFLAGS
		saved_LIBS=$LIBS
		CFLAGS="$CFLAGS $X_CFLAGS"
		LDFLAGS="$LDFLAGS $X_LIBS"
		LIBS="$LIBS $X_PRE_LIBS -lX11 $X_EXTRA_LIBS"
		AC_LINK_IFELSE([AC_LANG_CALL([],[XOpenDisplay])],[x11_works=yes],[x11_works=no])
		if test $x11_works = no ;then
			AC_MSG_WARN([could not link against X11, hoping Tk works without])
			CFLAGS=$saved_CFLAGS
			LDFLAGS=$saved_LDFLAGS
			LIBS=$saved_LIBS
		fi
		if test .$tkversion = .yes; then
			tkversion=""
			for version in $TK_VERSION tk8.5 tk8.4 tk8.3 tk8.2 tk; do
				ES_ADDPATH_CHECK_LIB($version, Tk_Init, [tkversion=$version], [])
				if test .$tkversion != .; then
					break;
				fi
			done
		else
			ES_ADDPATH_CHECK_LIB($tkversion, Tk_Init, [], [tkversion=""])
		fi
		if test .$tkversion = .; then
			case $target_os in
			*darwin*) AC_MSG_ERROR([If you have Tk installed, make sure that in one of the library paths, e.g. /usr/local/lib,
	there is a link from lib<tkversion>.dylib to the Tk library, which usually is
	/Library/Frameworks/Tk.framework/Tk.]) ;;
			*)	AC_MSG_ERROR([Tk library $tkversion not found]) ;;
			esac
		fi
		if test $tkversion = tk; then
			if test .$tclversion != .tcl; then
				AC_MSG_WARN([You are using a generic Tk version, but a defined Tcl version. This may cause problems.
	Try --tcl-version=tcl to also use a generic Tcl version, which may fit better.])
			fi
		fi
		case $target_os in
		*darwin*) extrapaths=/Library/Frameworks/Tk.framework/Headers ;;
		*) ;;
		esac
		ES_ADDPATH_CHECK_HEADER(tk.h, [], [AC_MSG_ERROR(Tk headers not found)],$extrapaths)
		AC_DEFINE(TK)
	fi
])

AC_DEFUN([ES_CHECK_FFTW3],[
 	ES_ADDPATH_CHECK_LIB(fftw3, fftw_plan_many_dft, [fftw3_found=yes], [fftw3_found=no])
	if test .$fftw3_found = .yes; then
		ES_ADDPATH_CHECK_HEADER(fftw3.h, [], [fftw3_found=no])
	fi
])

AC_DEFUN([ES_CHECK_FFTW2],[
	dnl we just assume that fftw and rfftw are in the same directory
	dnl if this is not the case for you, consider cleaning up your system...
	dnl first we check for the d* (SuSE) versions, then the normal ones
	dnl At the end we have to include rfftw before fftw on some systems, but testing
	dnl is only possible for fftw
	saved_LIBS=$LIBS
 	ES_ADDPATH_CHECK_LIB(dfftw, fftw_create_plan_specific, [fftw2_found=yes], [fftw2_found=no])
	if test .$fftw2_found = .yes; then
		LIBS="$saved_LIBS -ldrfftw -ldfftw"
	else
	 	ES_ADDPATH_CHECK_LIB(fftw, fftw_create_plan_specific, [fftw2_found=yes], [fftw2_found=no])
		if test .$fftw2_found = .yes; then
			LIBS="$LIBS -lrfftw -lfftw"
		fi
	fi
	if test .$fftw2_found = .yes; then
		ES_ADDPATH_CHECK_HEADER(fftw.h, [], [fftw2_found=no])
	fi
])

AC_DEFUN([ES_CHECK_FFTW],[
	AC_ARG_ENABLE(fftw,AC_HELP_STRING([--enable-fftw=<version>],[specify the version of FFTW to use (2 or 3)]),
		[], [enable_fftw=guess])
	if test $enable_mpi = guess; then
		if test .$use_fftw != .; then
			enable_mpi=$use_fftw
		fi
	fi
	if test $enable_fftw = 3; then
		AC_DEFINE(USEFFTW3)
		ES_CHECK_FFTW3
		if test .$fftw3_found != .yes; then
			AC_MSG_ERROR([could not link against FFTW3, please specify its header and library locations in CFLAGS and LDFLAGS])
		fi
	elif test $enable_fftw = 2; then
		ES_CHECK_FFTW2
		if test .$fftw2_found != .yes; then
			AC_MSG_ERROR([could not link against FFTW2, please specify its header and library locations in CFLAGS and LDFLAGS])
		fi
	else
		ES_CHECK_FFTW3
		if test .$fftw3_found = .yes; then
			enable_fftw=3
			AC_DEFINE(USEFFTW3)
		else
			ES_CHECK_FFTW2
			if test .$fftw2_found != .yes; then
				AC_MSG_ERROR([no FFTW found])
			fi
			enable_fftw=2
		fi
	fi
])
