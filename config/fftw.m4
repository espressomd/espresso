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
AC_DEFUN([ES_CHECK_FFTW],[
	AC_ARG_WITH([fftw],
		AC_HELP_STRING([--with-fftw=VERSION],
			[specify the version of FFTW to use (2 or 3)]),
		, with_fftw=guess)
	dnl with_fftw=no    don't use FFTW
	dnl with_fftw=yes   try to find a working FFTW, bail out if none is found
	dnl with_fftw=  (not set) try to find a working FFTW, continue if none is found
	dnl otherwise       use the specified version

        LIBS=" $LIBS -lm "

	if test x$with_fftw = xguess || test x$with_fftw = xyes; then
     		# determine fftw automatically
		ES_CHECK_FFTW3
		if test x$fftw3_found = xyes; then
		   use_fftw=3
		else
		  ES_CHECK_FFTW2
		  if test x$fftw2_found = xyes; then
		     use_fftw=2
		  elif test x$with_fftw = xguess; then
		      # no fftw was found, and none was requested
		      use_fftw=none
		  else
                      # no fftw was found, but the user wanted one (--with-fftw)
		      ES_NOTE_64BIT
		      AC_MSG_FAILURE([
********************************************************************************
* Could not find FFTW!                                                         *
********************************************************************************
])
		  fi
		fi
	elif test x$with_fftw = x3; then
                use_fftw=3
		ES_CHECK_FFTW3
		if test x$fftw3_found != xyes; then
		   ES_NOTE_64BIT
		   AC_MSG_FAILURE([
********************************************************************************
* Could not find FFTW3!                                                        *
* Please add the library path to LDFLAGS (e.g. configure LDFLAGS=-L/usr/lib)   *
* Please add the include path to CPPFLAGS                                      *
* (e.g. configure CPPFLAGS=-I/usr/include).                                    *
********************************************************************************
])
		fi
	elif test x$with_fftw = x2; then
		use_fftw=2
		ES_CHECK_FFTW2
		if test x$fftw2_found != xyes; then
		   ES_NOTE_64BIT
		   AC_MSG_FAILURE([
********************************************************************************
* Could not find FFTW2!                                                        *
* Please add the library path to LDFLAGS (e.g. configure LDFLAGS=-L/usr/lib)   *
* Please add the include path to CPPFLAGS                                      *
* (e.g. configure CPPFLAGS=-I/usr/include).                                    *
********************************************************************************
])
		fi
	elif test x$with_fftw = xno; then
	     use_fftw=none
	else
	  AC_MSG_ERROR([specified bad FFTW version ($with_fftw)])
	fi

	# after this, use_fftw should be set correctly
	if test x$use_fftw != xnone; then
		# save the result
		AC_DEFINE_UNQUOTED(FFTW, $use_fftw, [Whether to use the FFTW library, and which version to use])
	fi
])

AC_DEFUN([ES_CHECK_FFTW3],[
 	ES_ADDPATH_CHECK_LIB(fftw3, fftw_plan_many_dft, [fftw3_found=yes], [fftw3_found=no])
	if test x$fftw3_found = xyes; then
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
	if test x$fftw2_found = xyes; then
		LIBS="$saved_LIBS -ldrfftw -ldfftw"
	else
	 	ES_ADDPATH_CHECK_LIB(fftw, fftw_create_plan_specific, [fftw2_found=yes], [fftw2_found=no])
		if test x$fftw2_found = xyes; then
			LIBS="$LIBS -lrfftw -lfftw"
		fi
	fi
	if test x$fftw2_found = xyes; then
		ES_ADDPATH_CHECK_HEADER(fftw.h, [], [fftw2_found=no])
	fi
])

