AC_DEFUN([KC_BLADE],[
	AC_MSG_NOTICE([using preconfigured settings for the Blade Center Garching])
        use_mpi=mpich
        CC=mpicci
        CFLAGS="-I/afs/ipp/@sys/soft/fftw/f95i/include"
        LDFLAGS="-L/afs/ipp/@sys/soft/fftw/f95i/lib -L/afs/ipp/@sys/soft/gnu/lib"
        use_fftw=2
])

AC_DEFUN([KC_DINO],[
	AC_MSG_NOTICE([using preconfigured settings for the Alphas at MPIP])
        use_mpi=dmpi
	CC=cc
        CFLAGS="-I/usr/local/include -I/sw/axp_osf_40/tcl8.2/include"
        LDFLAGS="-L/sw/axp_osf_40/tcl8.2/shlib"
	TCL_VERSION=tcl
	use_fftw=2
])

AC_DEFUN([ES_KNOWN_CONFIGS],[
	AC_ARG_ENABLE(config,AC_HELP_STRING(--enable-config=cfg,[use a predefined configuration for known brain damaged system setups, currently supported: blade (Garching Blade Center), psi (Garching Regattas), thalpha (DEC Alphas at MPIP)]),,)
	if test .$enable_config = .blade; then
		KC_BLADE
	elif test .$enable_config = .thalpha; then
		KC_DINO
	fi
])
