dnl -*- mode: autoconf -*-

AC_DEFUN([ES_INIT_CHOOSER],[
	# test whether to use the chooser or not
	AC_ARG_ENABLE(chooser,
		AC_HELP_STRING([--enable-chooser],
			[create a chooser for the Espresso binary [[no]]]
			),
		,[enable_chooser=no])
	if test .$enable_chooser = .yes; then
	   execpkglibdir='${pkglibdir}'/obj-`$srcdir/config/config.guess`
	   scriptsdir='${pkglibdir}/scripts'
	else
	   scriptsdir='${pkgdatadir}/scripts'
	fi
	AM_CONDITIONAL(USE_CHOOSER, test .$enable_chooser = .yes)
	AC_SUBST(scriptsdir)
	AC_SUBST(execpkglibdir)
])

AC_DEFUN([ES_INIT_MYCONFIG],[
	# Handling the myconfig-header
	AC_ARG_WITH(myconfig, 
		AC_HELP_STRING(--with-myconfig=FILE,
			[name of the local config file [[myconfig.h]]]),
		[ test .$with_myconfig = .no && with_myconfig=myconfig.h],
		[ with_myconfig=myconfig.h ])
	MYCONFIG_H=$with_myconfig
	AC_SUBST(MYCONFIG_H)
	AC_DEFINE_UNQUOTED(MYCONFIG_H,"$with_myconfig",[The name of the local configure header])
	AH_BOTTOM([
#ifdef MYCONFIG_H
#include MYCONFIG_H
#endif
	])
])

AC_DEFUN([ES_CHECK_EFENCE],[
	# check for efence
	AC_ARG_WITH(efence,
		AC_HELP_STRING([--with-efence],[use ElectricFence memory debugging for the debug binary]),
		,with_efence=no)
	AC_ARG_ENABLE(efence,,
		      [
		      with_efence=$enable_efence
		      AC_MSG_WARN([
********************************************************************************
* The option --enable-efence is deprecated and will be removed in a future     *
* version of Espresso! Please use --with-efence instead!                       *
********************************************************************************\
		       ])])
	if test $with_efence = yes; then
		AC_CHECK_LIB(efence,malloc,,AC_MSG_ERROR([could not link against the efence library]))
		CPPFLAGS="$CPPFLAGS -DEFENCE"
		LDFLAGS="$LDFLAGS -lefence"
	fi
])

AC_DEFUN([ES_CHECK_TCL],[
	AC_BEFORE([$0],[ES_CHECK_TK])

	# check for tcl
        AC_DEFINE([USE_NON_CONST],1,[prevent TCL from defining const ptrs])

	AC_ARG_WITH(tcl,AC_HELP_STRING([--with-tcl=VERSION],[specify the tcl library to use (e.g. tcl8.4)]),
		[tclversion=$with_tcl], [tclversion=yes])

	AC_ARG_ENABLE(tcl,,
		[
		tclversion=$enable_tcl
		AC_MSG_WARN([
********************************************************************************
* The option --enable-tcl is deprecated and will be removed in a future        *
* version of Espresso! Please use --with-tcl instead!                          *
********************************************************************************\
		])])

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
	AC_ARG_WITH(tk,
		AC_HELP_STRING([--with-tk=VERSION],
			[whether to use Tk and which version to use for the GUI]),
		[tkversion=$with_tk], [tkversion=no])
	AC_ARG_ENABLE(tk,,
		[
		tkversion=$enable_tk
		AC_MSG_WARN([
********************************************************************************
* The option --enable-tk is deprecated and will be removed in a future         *
* version of Espresso! Please use --with-tk instead!                           *
********************************************************************************\
		])])
	if test .$tkversion != .no; then
		dnl test for X11
		AC_PATH_XTRA
		saved_CPPFLAGS=$CPPFLAGS
		saved_LDFLAGS=$LDFLAGS
		saved_LIBS=$LIBS
		CPPFLAGS="$CPPFLAGS $X_CFLAGS"
		LDFLAGS="$LDFLAGS $X_LIBS"
		LIBS="$LIBS $X_PRE_LIBS -lX11 $X_EXTRA_LIBS"
		AC_LINK_IFELSE([AC_LANG_CALL([],[XOpenDisplay])],[x11_works=yes],[x11_works=no])
		if test $x11_works = no ;then
			AC_MSG_WARN([could not link against X11, hoping Tk works without])
			CPPFLAGS=$saved_CPPFLAGS
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
		AC_DEFINE(TK,,[Whether to use Tk])
	fi
])


dnl ******************** helper macros ********************
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
	save_CPPFLAGS=$CPPFLAGS
	adp_found=no
	dnl let's see whether it's in the default paths
	AC_COMPILE_IFELSE([
		#include <$1>
	],[adp_found=yes],[])

	if test .$adp_found = .no; then
		for path in $4 /sw/include /usr/include /usr/local/include /opt/include; do
			CPPFLAGS="$save_CPPFLAGS -I$path"
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
		CPPFLAGS=$save_CPPFLAGS
	fi
	AS_IF([test .$adp_found = .yes], [$2],[$3])
])

