AC_DEFUN([MPI_SETUP],[
	if test .$1 = .; then MPI_GUESS_ENV
	else
		case $1 in
		poe)	MPI_FIND_MPICC
			MPI_FIND_POE
			if test $poe_found != yes; then AC_MSG_ERROR([IBM POE MPI not working]); fi
			;;
		dmpi)	MPI_FIND_DMPI
			if test $dmpi_found != yes; then AC_MSG_ERROR([Tru64/OSF MPI not working]); fi
			;;
		generic) MPI_FIND_MPICC
			MPI_FIND_GENERIC
			if test $generic_found != yes; then AC_MSG_ERROR([generic MPI not working]); fi
			;;
		lam)	MPI_FIND_MPICC
			MPI_FIND_LAM
			if test $lam_found != yes; then	AC_MSG_ERROR([LAM/MPI not working]); fi
			;;
		mpich)	MPI_FIND_MPICH
			if test $mpich_found != yes; then AC_MSG_ERROR([MPICH not working]); fi
			;;
		esac
	fi
])

AC_DEFUN([MPI_GUESS_ENV],[
	enable_mpi=guess
	dnl look for special implementations
	case $target_os in
	*aix*)	MPI_FIND_MPICC
		mpicc_done=y
		MPI_FIND_POE
		if test .$poe_found = .yes; then enable_mpi=poe; fi
		;;
	*osf*)	MPI_FIND_DMPI
		if test .$dmpi_found = .yes; then enable_mpi=dmpi; fi
		;;
	esac

	dnl always check for generic, LAM/MPI and MPICH, as they are multiplatform

	if test .$enable_mpi = .guess; then
		if test .$mpicc_done = .; then
			MPI_FIND_MPICC
		fi
		MPI_FIND_LAM
		if test .$lam_found = .yes; then enable_mpi=lam; fi
	fi

	if test .$enable_mpi = .guess; then
		MPI_FIND_MPICH
		if test .$mpich_found = .yes; then enable_mpi=mpich; fi
	fi

	dnl generic last, since here we just assume how to run mpirun

	if test .$enable_mpi = .guess; then
		if test .$mpicc_done = .; then
			MPI_FIND_MPICC
		fi
		MPI_FIND_GENERIC
		if test .$generic_found = .yes; then enable_mpi=generic; fi
	fi

	dnl out of guesses
	if test .$enable_mpi = .guess; then
		AC_MSG_WARN([
********************************************************************
* could neither detect LAM nor MPICH nor a native MPI environment, *
* using the FAKE implementation for one processor only.            *
* If you have an MPI environment, please specify its type and      *
* the compiler or includes/libraries of your MPI implementation    *
* manually, or, even better, add your MPI environment to           * 
* config/mpi.m4                                                    *
********************************************************************
])
		enable_mpi=fake
		MPI_SETUP_FAKE
	fi
])

AC_DEFUN([MPI_SETUP_FAKE],[
	CFLAGS="$CFLAGS -I."
	AC_DEFINE(MPI,"fake",[Which MPI implementation to use])
	ADD_SOURCES="$ADD_SOURCES mpi"
	MPI_INVOCATION="\$ESPRESSO_SOURCE/obj-$target/Espresso_bin @ARGUMENTS@"
])

AC_DEFUN([MPI_FIND_MPICC],[
	# only test the compiler if not overridden by the user
	if test .$user_defined_CC != .yes; then
		AC_CHECK_PROGS(MPICC,[mpcc mpxlc mpicc mpicci mpiccg hcc cc icc gcc])
		CC=$MPICC
	fi
	AC_MSG_CHECKING([whether the $CC command works out of the box for MPI])
	AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY(MPI_Init)],
    		[AC_MSG_RESULT(yes); mpicc_works=yes],[AC_MSG_RESULT(no); mpicc_works=no])
])

AC_DEFUN([MPI_FIND_GENERIC],[
	AC_MSG_NOTICE([trying to find a generic MPI])
	generic_found=yes
	if test .$mpicc_works != .yes; then
		generic_found=no
	else
		AC_PATH_PROG(generic_bin, mpirun, mpirun, [])
		if test .$generic_bin = .; then
			generic_found=no
		fi
	fi
	if test $generic_found = yes; then
		AC_DEFINE(MPI,"generic")
		MPI_INVOCATION="$generic_bin -np @NP@ \$ESPRESSO_SOURCE/obj-$target/Espresso_bin @ARGUMENTS@"
		AC_MSG_NOTICE([found a generic MPI environment])
	fi
])

AC_DEFUN([MPI_FIND_LAM],[
	AC_MSG_NOTICE([trying to find a LAM environment])
	lam_found=yes
	if test $mpicc_works = yes; then
		lam_cflags=$CFLAGS
		lam_ldflags=$LDFLAGS
		lam_libs=$LIBS
	else
		save_cflags=$CFLAGS
		save_ldflags=$LDFLAGS
		save_libs=$LIBS

		AC_MSG_CHECKING([whether the $CC command works without additional flags])
		lam_cflags=$CFLAGS
		lam_ldflags=$LDFLAGS
		lam_libs="$LIBS -lmpi -llam -lpthread"
		LIBS=$lam_libs
		AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY(MPI_Init)],
    			[AC_MSG_RESULT(yes); mpicc_works_wl=yes],[AC_MSG_RESULT(no); mpicc_works_wl=no])
		if test .$mpicc_works_wl = .no; then
			dnl guess the library paths or the include paths.
			for hdrprf in /opt/lam /usr/lam /usr/local/lam /usr /usr/local; do
				if test -f $hdrprf/include/lam.h; then lam_hdr_prf=$hdrprf; lam_hdr=$hdrprf/include; break; fi
			done
			if test $lam_hdr. = .; then
				AC_MSG_NOTICE([did not find lam.h, if you want to use LAM, please specify its location via CFLAGS])
				lam_found=no
			fi
			for libprf in /opt/lam /usr/lam /usr/local/lam /usr /usr/local; do
				for lib in $libprf/lib64 $libprf/lib; do
					if test -f $lib/liblam.a; then lam_lib_prf=$libprf; lam_lib=$lib; break; fi
				done
				if test .$lam_lib != .; then break; fi
			done
			if test $lam_lib. = .; then
				AC_MSG_NOTICE([did not find lam library, if you want to use LAM, please specify its location and name via LDFLAGS])
				lam_found=no
			fi
			if test $lam_found = yes; then
				lam_cflags="$CFLAGS -I$lam_hdr"
				lam_ldflags="$LDFLAGS -L$lam_lib"
				CFLAGS=$lam_cflags
				LDFLAGS=$lam_ldflags
				AC_MSG_CHECKING([whether the $CC command works with guessed headers and libraries])
	 			AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY(MPI_Init)],
		    			AC_MSG_RESULT(yes),[AC_MSG_RESULT(no); lam_found=no])
			fi
		fi
		CFLAGS=$save_cflags
		LDFLAGS=$save_ldflags
		LIBS=$save_libs
	fi
	dnl last but not least the invocation, checking for mpirun
	AC_PATH_PROG(lam_bin, mpirun, mpirun, $PATH:$lam_lib_prf/bin:$lam_hdr_prf/bin)
	dnl check the program that it actually belongs to LAM
	$lam_bin -h > conftmp.log 2>&1
	if ! grep LAM conftmp.log >/dev/null 2>&1; then
		AC_MSG_NOTICE([mpirun binary for LAM not found, if you want to use LAM, please add it to your PATH])
		lam_found=no
	fi
	rm -f conftmp.log
	if test $lam_found = yes; then
		AC_DEFINE(MPI,"lam")
		CFLAGS=$lam_cflags
		LDFLAGS=$lam_ldflags
		LIBS=$lam_libs
		MPI_INVOCATION="$lam_bin -np @NP@ -nsigs \$ESPRESSO_SOURCE/obj-$target/Espresso_bin @ARGUMENTS@"
		AC_MSG_NOTICE([found a working LAM environment])
	fi
])

AC_DEFUN([MPI_FIND_MPICH],[
	AC_MSG_NOTICE([trying to find a MPICH environment])
	mpich_found=yes
	save_cflags=$CFLAGS
	save_ldflags=$LDFLAGS
	save_libs=$LIBS

	AC_MSG_CHECKING([whether the $CC command works without flags])
	mpich_cflags=$CFLAGS
	mpich_ldflags=$LDFLAGS
	mpich_libs="$LIBS -lmpich"
	LIBS=$mpich_libs
	AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY(MPI_Init)],
    		[AC_MSG_RESULT(yes); cc_works_wl=yes],[AC_MSG_RESULT(no); cc_works_wl=no])
	if test .$cc_works_wl = .no; then
		dnl guess the library paths or the include paths.
		for hdrpprf in /opt/mpich /usr/mpich /usr/local/mpich /usr /usr/local; do
			dnl test directly ...
			if test -f $hdrpprf/include/mpi.h; then mpich_hdr_prf=$hdrpprf; mpich_hdr=$hdrpprf/include; break; fi
			dnl ... but it may as well be a subdir
			for hdrprf in $hdrpprf/*; do
				if test -d $hdrprf; then
					if test -f $hdrprf/include/mpi.h; then mpich_hdr_prf=$hdrprf; mpich_hdr=$hdrprf/include; break; fi
				fi
				if test .$mpich_hdr != .; then break; fi
			done
			if test .$mpich_hdr != .; then break; fi
		done
		if test $mpich_hdr. = .; then
			AC_MSG_NOTICE([did not find mpi.h, please specify its location via CFLAGS])
			mpich_found=no
		fi
		for libpprf in /opt/mpich /usr/mpich /usr/local/mpich /usr /usr/local; do
			dnl test directly ...
			for lib in $libpprf/lib64 $libpprf/lib; do
				if test -f $lib/libmpich.a; then mpich_lib_prf=$libpprf; mpich_lib=$lib; break; fi
			done
			dnl ... but it may as well be a subdir
			for libprf in $libpprf/*; do
				if test -d $libprf; then
					for lib in $libprf/lib64 $libprf/lib; do
						if test -f $lib/libmpich.a; then mpich_lib_prf=$libprf; mpich_lib=$lib; break; fi
					done
				fi
				if test .$mpich_lib != .; then break; fi
			done
			if test .$mpich_lib != .; then break; fi			
		done
		if test $mpich_lib. = .; then
			AC_MSG_NOTICE([did not find mpich library, please specify its location and name via LDFLAGS])
			mpich_found=no
		fi
		if test $mpich_found = yes; then
			mpich_cflags="$CFLAGS -I$mpich_hdr"
			mpich_ldflags="$LDFLAGS -L$mpich_lib"
			CFLAGS=$mpich_cflags
			LDFLAGS=$mpich_ldflags
			AC_MSG_CHECKING([whether the $CC command works with guessed headers and libraries])
	 		AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY(MPI_Init)],
		    		AC_MSG_RESULT(yes),[AC_MSG_RESULT(no); mpich_found=no])
			CFLAGS=$save_cflags
			LDFLAGS=$save_ldflags
		fi
	fi
	CFLAGS=$save_cflags
	LDFLAGS=$save_ldflags
	LIBS=$save_libs

	dnl last but not least the invocation, checking for mpirun
	AC_PATH_PROG(mpich_bin, mpirun, mpirun, $mpich_lib_prf/bin:$mpich_hdr_prf/bin:$PATH)
	dnl check the program that it actually belongs to MPICH
	$mpich_bin -h > conftmp.log 2>&1
	if ! grep mpich conftmp.log >/dev/null 2>&1; then
		AC_MSG_NOTICE([mpirun not found or not from MPICH, please make it first in your PATH])
		mpich_found=no
	fi
	rm -f conftmp.log
	if test $mpich_found = yes; then
		AC_DEFINE(MPI,"mpich")
		CFLAGS="$mpich_cflags"
		LDFLAGS=$mpich_ldflags
		LIBS=$mpich_libs
		MPI_INVOCATION="$mpich_bin -np @NP@ \$ESPRESSO_SOURCE/obj-$target/Espresso_bin @ARGUMENTS@"
		AC_MSG_NOTICE([found a working MPICH environment])
	fi
])

AC_DEFUN([MPI_FIND_POE],[
	poe_found=yes
	if test .$mpicc_works != .yes; then
		AC_MSG_NOTICE([the poe environment only works via mpcc, please specify your c compiler])
		poe_found=no
	fi
	if test .$poe_found = .yes; then
		AC_PATH_PROG(poe_bin, poe, poe, $PATH)
		AC_DEFINE(MPI,"poe")
		MPI_INVOCATION="$poe_bin \$ESPRESSO_SOURCE/obj-$target/Espresso_bin @ARGUMENTS@ -procs @NP@"
	fi
])

AC_DEFUN([MPI_FIND_DMPI],[
	saved_libs=$LIBS
	dmpi_found=yes
	AC_MSG_CHECKING([whether the $CC command works for MPI])
	dmpi_libs="$LIBS -lmpi"
	LIBS=$dmpi_libs
	AC_LINK_IFELSE([AC_LANG_FUNC_LINK_TRY(MPI_Init)],
    			[AC_MSG_RESULT(yes)],[AC_MSG_RESULT(no); dmpi_found=no])
	if test .$dmpi_found = .yes; then
		AC_PATH_PROG(dmpi_bin, dmpirun, dmpirun, $PATH)
		AC_DEFINE(MPI, "dmpi")
		MPI_INVOCATION="$dmpi_bin -np @NP@ \$ESPRESSO_SOURCE/obj-$target/Espresso_bin @ARGUMENTS@"
	fi
])
