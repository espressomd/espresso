dnl -*- mode: autoconf -*-

AC_DEFUN([ES_CHECK_FFTW],[
	AC_ARG_WITH(fftw,
		AC_HELP_STRING([--with-fftw=VERSION],
		[specify the version of FFTW to use (2 or 3)]),
		[], [with_fftw=guess])
	AC_ARG_ENABLE(fftw,,
		[
		with_fftw=$enable_fftw
		AC_MSG_WARN([
********************************************************************************
* The option --enable-fftw is deprecated and will be removed in a future       *
* version of Espresso! Please use --with-fftw instead!                         *
********************************************************************************\
		])])

	if test $with_fftw = guess; then
		if test .$known_fftw != .; then
			with_fftw=$known_fftw
		fi
	fi
	if test $with_fftw = 3; then
		AC_DEFINE(USEFFTW3,,[Whether to use the FFTW3 library])
		ES_CHECK_FFTW3
		if test .$fftw3_found != .yes; then
			AC_MSG_ERROR([could not link against FFTW3, please specify its header and library locations in CPPFLAGS and LDFLAGS])
		fi
	elif test $with_fftw = 2; then
		ES_CHECK_FFTW2
		if test .$fftw2_found != .yes; then
			AC_MSG_ERROR([could not link against FFTW2, please specify its header and library locations in CPPFLAGS and LDFLAGS])
		fi
	else
		ES_CHECK_FFTW3
		if test .$fftw3_found = .yes; then
			with_fftw=3
			AC_DEFINE(USEFFTW3)
		else
			ES_CHECK_FFTW2
			if test .$fftw2_found != .yes; then
				AC_MSG_ERROR([no FFTW found])
			fi
			with_fftw=2
		fi
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

