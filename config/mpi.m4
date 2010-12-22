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
dnl

dnl Recognize the MPI compiler
AC_DEFUN([ES_CHECK_MPI],[
  dnl The following is needed to avoid a "macro is expanded before it
  dnl is required" warning. The function _ES_CHECK_MPI will check for
  dnl the compiler
  AC_REQUIRE([_ES_CHECK_MPI])

  AS_IF([test x"$with_mpi" != xno], [
    # test whether MPI_Init is available now
    AC_CHECK_FUNC(MPI_Init,, [
      # if not, try to find it in a library
      AC_SEARCH_LIBS(MPI_Init, [mpi mpich],, [
        # if not, give up or use fake
        if test xyes = x"$with_mpi"; then
          AC_MSG_FAILURE([MPI compiler requested, but couldn't compile with MPI.])
        else
          AC_MSG_WARN([No MPI compiler found, will use fake implementation!])
          use_mpi_fake="yes"
        fi
      ])
    ])

    # test wether the MPI headers are there
    AS_IF([test x"$use_mpi_fake" != xyes], [
      AC_MSG_CHECKING([for mpi.h])
      AC_TRY_COMPILE([#include <mpi.h>],,[
        AC_MSG_RESULT(yes)
      ], [
        AC_MSG_RESULT(no)
	AC_MSG_ERROR([Header mpi.h was not found!])
      ])
    ])
  ])

  # if requested, use the fake implementation
  AM_CONDITIONAL(MPI_FAKE, [test x"$use_mpi_fake" = xyes])

  # determine ESPRESSO_MPIEXEC
  AS_IF([test x"$use_mpi_fake" = xyes], [
    ESPRESSO_MPIEXEC=""
  ],[
    # mpiexec executable
    AC_ARG_VAR([MPIEXEC], [MPI command mpiexec])
    AS_IF([test x"$MPIEXEC" = x], [
      AC_PATH_PROG([MPIEXEC], [mpiexec], [no])
    ])

    AC_MSG_CHECKING([for the mympiexec user script])
    AC_ARG_WITH([mympiexec],
      AS_HELP_STRING([--with-mympiexec@<:@=SCRIPT@:>@],
        [specify the mpiexec-like program or script that should be
         used to run ESPResSo in parallel. If the script doesn't
         exist, it will try to use mpiexec. 
	 @<:@SCRIPT=./mympiexec.sh@:>@]),
      [if test x"$with_mympiexec" = xno; then
         MYMPIEXEC=""
       else
         MYMPIEXEC="$with_mympiexec"
         dir=`AS_DIRNAME([$MYMPIEXEC])`
         if test x"$dir" = x.; then
           MYMPIEXEC="`pwd`/$MYMPIEXEC"
         fi
       fi],
      [ MYMPIEXEC="`pwd`/mympiexec.sh" ])
    AC_MSG_RESULT($MYMPIEXEC)
    AC_SUBST(MYMPIEXEC)

    ESPRESSO_MPIEXEC="`pwd`/tools/es_mpiexec"
 ])
])

dnl We need to split the main function (ES_CHECK_MPI) into two parts
dnl to avoid an "macro is expanded before it is required" warning.
AC_DEFUN([_ES_CHECK_MPI],[
  AC_MSG_CHECKING([whether to compile using MPI])
  # Check for --with-mpi
  AC_ARG_WITH(mpi, [AC_HELP_STRING([--with-mpi],
      [compile with MPI (parallelization) support. If none is found, the
       fake implementation for only one processor is used. Default: guess])
  ], [], [
    with_mpi=auto
  ])
  AC_MSG_RESULT($with_mpi)

  AC_ARG_VAR(MPICC,[MPI C compiler command])
  # if MPI is wanted, look for MPI compiler
  if test x"$with_mpi" != xno; then
    if test -z "$CC" && test -n "$MPICC"; then
      CC="$MPICC"
    else
      AC_CHECK_TOOLS([CC], [mpicc hcc mpxlc_r mpxlc mpcc cmpicc])
    fi
  else
    use_mpi_fake="yes"
  fi
  AC_PROG_CC
])

