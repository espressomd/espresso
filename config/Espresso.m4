dnl -*- mode: autoconf -*-
dnl Copyright (C) 2010 The ESPResSo project
dnl Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
dnl  
dnl This file is part of ESPResSo.
dnl  
dnl ESPResSo is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl  
dnl ESPResSo is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl  
dnl You should have received a copy of the GNU General Public License
dnl along with this program.  If not, see <http://www.gnu.org/licenses/>. 
AC_DEFUN([ES_INIT_CHOOSER],[
	# test whether to use the chooser or not
	AC_ARG_ENABLE(chooser,
		AC_HELP_STRING([--enable-chooser],
			[create a chooser for the Espresso binary [[no]]]
			),
		,[enable_chooser=no])
	AC_MSG_CHECKING([whether to enable the chooser])
	if test .$enable_chooser = .yes; then
	   AC_MSG_RESULT(yes)
	   execpkglibdir='${pkglibdir}'/obj-`$srcdir/config/config.guess`
	   scriptsdir='${pkglibdir}/scripts'
	else
	   AC_MSG_RESULT(no)
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
			[default name of the local config file [[myconfig.h]]]),
		[ test .$with_myconfig = .no && with_myconfig=myconfig.h],
		[ with_myconfig=myconfig.h ])
	AC_MSG_CHECKING([the name of the configuration header])
	myconfig=$with_myconfig
	AC_MSG_RESULT($myconfig)
	AC_SUBST(myconfig)
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
	if test .$with_efence = .yes; then
		AC_CHECK_LIB(efence,malloc,,AC_MSG_FAILURE([could not link against the efence library]))
		CPPFLAGS="$CPPFLAGS -DEFENCE"
		LDFLAGS="$LDFLAGS -lefence"
	fi
])

AC_DEFUN([ES_CHECK_TCL],[
	AC_BEFORE([$0],[ES_CHECK_TK])

	# check for tcl
        AC_DEFINE([USE_NON_CONST],1,[prevent TCL from defining const ptrs])

	AC_ARG_WITH(tcl,AC_HELP_STRING([--with-tcl=VERSION],[specify the tcl library to use (e.g. tcl8.4)]),
	            [], [with_tcl=yes])

	AC_ARG_ENABLE(tcl,,
		[
		with_tcl=$enable_tcl
		AC_MSG_WARN([
********************************************************************************
* The option --enable-tcl is deprecated and will be removed in a future        *
* version of Espresso! Please use --with-tcl instead!                          *
********************************************************************************\
		])])

	dnl with_tcl=yes  try to find a working TCL version, bail out if none is found
	dnl with_tcl=no   bail out
	dnl otherwise     use the specified version

	if test .$with_tcl = .no; then
	   AC_MSG_ERROR([Tcl is required by ESPResSo!])
	elif test .$with_tcl = .yes; then
		for version in $TCL_VERSION tcl8.6 tcl8.5 tcl8.4 tcl8.3 tcl8.2 tcl; do
			ES_ADDPATH_CHECK_LIB($version, Tcl_Init, [use_tcl=$version])
			if test .$use_tcl != .; then break; fi
		done
	else
		ES_ADDPATH_CHECK_LIB($with_tcl, Tcl_Init, [use_tcl=$with_tcl])
	fi

	# check the results
	if test .$use_tcl = .; then
	   # Tcl was not found
	   case $target_os in
		*darwin*) AC_MSG_NOTICE(
[If you have Tcl installed, make sure that in one of the library paths, e.g. /usr/local/lib,
there is a link from lib<tclversion>.dylib to the Tcl library, which usually is
/Library/Frameworks/Tcl.framework/Tcl.]) ;;
		*) ES_NOTE_64BIT
		   AC_MSG_FAILURE([
********************************************************************************
* Could not link against the (static) Tcl library (libtcl*.a).                 *
* Please add the library path to LDFLAGS (e.g. configure LDFLAGS=-L/usr/lib)!  *
********************************************************************************
]) ;;
		esac
	fi
	case $target_os in
	*darwin*) extrapaths=/Library/Frameworks/Tcl.framework/Headers ;;
	*linux*)  # path used by *buntu
		  extrapaths=/usr/include/$version ;;
	esac
	ES_ADDPATH_CHECK_HEADER(tcl.h, [], 
		[AC_MSG_FAILURE([
********************************************************************************
* Could not find the Tcl header files (tcl.h).                                 *
* Please add the include path to CPPFLAGS                                      *
* (e.g. configure CPPFLAGS=-I/usr/include)!                                    *
********************************************************************************
])],
		$extrapaths)

	if test .$use_tcl = .; then
	   use_tcl=none
	fi
])


AC_DEFUN([ES_CHECK_TK],[
	AC_ARG_WITH(tk,
		AC_HELP_STRING([--with-tk=VERSION],
			[whether to use Tk, and which version to use]),
		[], [with_tk=no])
	AC_ARG_ENABLE(tk,,
		[
		with_tk=$enable_tk
		AC_MSG_WARN([
********************************************************************************
* The option --enable-tk is deprecated and will be removed in a future         *
* version of Espresso! Please use --with-tk instead!                           *
********************************************************************************\
		])])
	dnl with_tk=no   don't use Tk
	dnl with_tk=yes  try to find a working Tk version, bail out if none is found
	dnl otherwise    use the specified version
	if test .$with_tk != .no; then
		# test for X11
		AC_PATH_XTRA
		saved_CPPFLAGS=$CPPFLAGS
		saved_LDFLAGS=$LDFLAGS
		saved_LIBS=$LIBS
		CPPFLAGS="$CPPFLAGS $X_CFLAGS"
		LDFLAGS="$LDFLAGS $X_LIBS"
		LIBS="$LIBS $X_PRE_LIBS -lX11 $X_EXTRA_LIBS"
		AC_LINK_IFELSE([AC_LANG_CALL([],[XOpenDisplay])],[x11_works=yes],[x11_works=no])
		if test $x11_works = no; then
			AC_MSG_WARN([could not link against X11, hoping Tk works without])
			CPPFLAGS=$saved_CPPFLAGS
			LDFLAGS=$saved_LDFLAGS
			LIBS=$saved_LIBS
		fi
		# now test whether Tk can be found
		if test .$with_tk = .yes; then
			for version in $TK_VERSION tk8.5 tk8.4 tk8.3 tk8.2 tk; do
				ES_ADDPATH_CHECK_LIB($version, Tk_Init, [use_tk=$version], [])
				if test .$use_tk != .; then break; fi
			done
		else
			ES_ADDPATH_CHECK_LIB($with_tk, Tk_Init, [use_tk=$with_tk], [])
		fi
		if test .$use_tk = .; then
			case $target_os in
			(*darwin*) AC_MSG_ERROR([If you have Tk installed, make sure that in one of the library paths, e.g. /usr/local/lib,
	there is a link from lib<tkversion>.dylib to the Tk library, which usually is
	/Library/Frameworks/Tk.framework/Tk.]) ;;
			(*) ES_NOTE_64BIT
			    AC_MSG_FAILURE([Tk library $with_tk not found]) ;;
			esac
		fi
		if test .$use_tk = .tk; then
			if test .$use_tcl != .tcl; then
				AC_MSG_WARN([You are using a generic Tk version, but a defined Tcl version. This may cause problems.
	Try --with-tcl=tcl to also use a generic Tcl version, which may fit better.])
			fi
		fi
		case $target_os in
		*darwin*) extrapaths=/Library/Frameworks/Tk.framework/Headers ;;
		*linux*)  # path used by *buntu
			  extrapaths="/usr/include/$version /usr/include/$use_tcl" ;;
		(*) ;;
		esac
		ES_ADDPATH_CHECK_HEADER(tk.h, [], 
			[AC_MSG_ERROR([Tk headers not found. Please add the include path to CPPFLAGS (e.g. configure CPPFLAGS=-I/usr/include/tcl8.4).])
			]
			,$extrapaths)
		AC_DEFINE_UNQUOTED(TK,$use_tk,[Whether to use Tk])
	else
		use_tk=none
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
	AC_COMPILE_IFELSE([AC_LANG_SOURCE([
		#include <$1>
	])],[adp_found=yes],[])

	if test .$adp_found = .no; then
		for path in $4 /sw/include /usr/include /usr/local/include /opt/include; do
			CPPFLAGS="$save_CPPFLAGS -I$path"
			AC_COMPILE_IFELSE([AC_LANG_SOURCE([
				#include <$1>
			])],[adp_found=yes],[])
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


