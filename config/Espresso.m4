dnl -*- mode: autoconf -*-
dnl Copyright (C) 2010,2011,2012,2013 The ESPResSo project
dnl Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
dnl   Max-Planck-Institute for Polymer Research, Theory Group
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


