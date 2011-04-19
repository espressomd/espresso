dnl -*- mode: autoconf -*-
# ===========================================================================
#       http://www.gnu.org/software/autoconf-archive/ax_mpi_only.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_MPI_ONLY([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND[, MPI-WANTED-TEST]]])
#
# DESCRIPTION
#
#   This macro tries to find out how to compile programs that use MPI
#   (Message Passing Interface), a standard API for parallel process
#   communication (see http://www-unix.mcs.anl.gov/mpi/)
#
# LICENSE
#
#   Copyright (c) 2010 Olaf Lenz <olenz@icp.uni-stuttgart.de>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 1

AC_DEFUN([AX_MPI_ONLY], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

# Check for compiler
# Needs to be split off into an extra macro to ensure right expansion
# order.
AC_REQUIRE([_AX_MPI_ONLY],[_AX_MPI_ONLY([$1])])

AS_IF([test x"$_ax_mpi_only_mpi_wanted" = xno], 
  [ ax_mpi_only_mpi_found=no ],
  [
    # test whether MPI_Init is available in one oft he libraries
    AC_SEARCH_LIBS(MPI_Init, [mpi mpich],
      [ ax_mpi_only_mpi_found=yes ],
      [ ax_mpi_only_mpi_found=no ])

    # Check for header
    AS_IF([test x"$ax_only_mpi_found" = xyes], [
      AC_MSG_CHECKING([for mpi.h])
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([#include <mpi.h>])],
        [ AC_MSG_RESULT(yes)], 
        [ AC_MSG_RESULT(no)
	  ax_mpi_only_mpi_found=no
      ])
    ])
])

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
AS_IF([test x"$ax_mpi_only_mpi_found" = xyes], [
        ifelse([$2],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$2])
        :
],[
        $3
        :
])

])dnl AX_MPI_ONLY

dnl _AX_MPI_ONLY is an internal macro required by AX_MPI_ONLY.
dnl To ensure the right expansion order, the main function AX_MPI_ONLY
dnl has to be split into two parts.
AC_DEFUN([_AX_MPI_ONLY], [
  AC_ARG_VAR(MPICC,[MPI C compiler command])
  ifelse([$1],,[_ax_mpi_only_mpi_wanted=yes],[
    AC_MSG_CHECKING([whether to compile using MPI])
    if $1; then
      _ax_mpi_only_mpi_wanted=yes
    else
      _ax_mpi_only_mpi_wanted=no
    fi
    AC_MSG_RESULT($_ax_mpi_only_mpi_wanted)
  ])
  if test x"$_ax_mpi_only_mpi_wanted" = xyes; then
    if test -z "$CC" && test -n "$MPICC"; then
      CC="$MPICC"
    else
      AC_CHECK_TOOLS([CC], [mpicc hcc mpxlc_r mpxlc mpcc cmpicc])
    fi
  fi
  AC_PROG_CC

])dnl _AX_MPI_ONLY
