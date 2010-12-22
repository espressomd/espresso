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
  # Check for --with-mpi
  AC_ARG_WITH(mpi, [AC_HELP_STRING([--with-mpi],
      [compile with MPI (parallelization) support. If none is found, the
       fake implementation for only one processor is used. Default: auto])
  ],,[with_mpi=auto])

  AX_MPI_ONLY([test x"$with_mpi" != xno],,[
    if test xyes = x"$with_mpi"; then
      AC_MSG_FAILURE([MPI compiler requested, but couldn't compile with MPI.])
    else
      AC_MSG_WARN([No MPI compiler found, will use fake implementation!])
      use_mpi_fake="yes"
    fi
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

