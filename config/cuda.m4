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
AC_DEFUN([ES_CHECK_CUDA],[
	# search for the nvcc CUDA compiler
	AC_ARG_WITH([cuda],
		AC_HELP_STRING([--with-cuda=CUDAINSTALLDIR],
			[specify where CUDA is installed.
			The cuda compiler can also be specified by setting the NVCC
			environment variable. The CUDA library and header can be
			manually specified by using CPPFLAGS, LDFLAGS and LIBS.]),
		, with_cuda=no)

	if test x$with_cuda != xno; then
	   # if installation dir is given, set the paths
	   if test x$with_cuda != xyes; then
	       if test x$NVCC = x; then
		   NVCC=$with_cuda/bin/nvcc
	       fi
	       if test -d $with_cuda/lib64; then
	           LDFLAGS="$LDFLAGS -L$with_cuda/lib64"
	       else
		   LDFLAGS="$LDFLAGS -L$with_cuda/lib"
	       fi
	       NVCCFLAGS="$NVCCFLAGS -I$with_cuda/include"
           fi

	   # NVCC
	   AC_ARG_VAR(NVCC,[NVIDIA CUDA compiler command])
	   AC_ARG_VAR(NVCCFLAGS,[special compiler flags for the NVIDIA CUDA compiler])
	   AC_PATH_PROG(NVCC, nvcc, no)
	   if test x$NVCC = xno; then
               AC_MSG_FAILURE([CUDA compiler nvcc was not found, specify location using the NVCC variable])
	   fi

	   # runtime library
	   AC_CHECK_LIB(cudart, cudaGetDevice, [LIBS="$LIBS -lcudart"], [AC_MSG_FAILURE([could not find cuda runtime library (cudart), specify location using LDFLAGS])])

	   # NVCC compile check
	   AC_MSG_CHECKING([whether CUDA works])
	   # if no other compute capability is defined by the user, we need at least 1.1
	   case "$NVCCFLAGS" in
	       *-arch=*) ;;
	       *) NVCCFLAGS="$NVCCFLAGS --ptxas-options=-v -arch=compute_20 -code=sm_20"
	   esac
	   # use nvcc
	   save_CC=$CC
	   save_CFLAGS=$CFLAGS
	   CC=$NVCC
	   CFLAGS="$NVCCFLAGS -x cu"
	   AC_COMPILE_IFELSE([AC_LANG_PROGRAM([
		#include <cuda.h>
		__global__ void test()
		{}
	   ],[
		dim3 block(1), grid(1);
		test<<<block,grid>>>();
	   ])], [AC_MSG_RESULT(yes)], [AC_MSG_FAILURE([cannot compile CUDA code. Look at config.log for more details.])])
	   CC=$save_CC
	   CFLAGS=$save_CFLAGS
	   AC_DEFINE(CUDA,[],[Whether CUDA is available])
	fi

	AM_CONDITIONAL(CUDA, [test x$with_cuda != xno])
])
