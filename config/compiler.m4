dnl -*- mode: autoconf -*-
dnl This file contains macros that test for different features of the
dnl C-compiler and linker.
dnl
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

AC_DEFUN([ES_CHECK_COMPILER],[
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

	ES_TRY_ADD_CFLAG(-O5)
	ES_CHECK_INLINING
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


dnl ******************************* linking ********************************
AC_DEFUN([ES_CHECK_LINK],[
	case $target_os in
	*darwin*) LD="`pwd`/darwinlink.sh $CC" ;;
	*) LD=$CC ;;
	esac
])
