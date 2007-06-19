dnl -*- mode: autoconf -*-
dnl This file defines known configurations.
dnl Sets the variables CC, CPPFLAGS, LDFLAGS, known_mpi and known_fftw
AC_DEFUN([ES_KNOWN_CONFIGS],[
	AC_BEFORE([$0],[ES_CHECK_MPI])
	AC_ARG_ENABLE(config,
		AC_HELP_STRING(--enable-config=cfg,
			[use a predefined configuration for known
			brain damaged system setups, currently
			supported: blade (Garching Blade Center), psi
			(Garching Regattas), thalpha (DEC Alphas at
			MPIP)]),,) 
	if test .$enable_config = .blade; then
		ES_KNOWN_BLADE
	elif test .$enable_config = .thalpha; then
		ES_KNOWN_DINO
	fi
])

AC_DEFUN([ES_KNOWN_BLADE],[
	AC_MSG_NOTICE([using preconfigured settings for the Blade Center Garching])
        known_mpi=mpich
        CC=mpicci
        CPPFLAGS="-I/afs/ipp/@sys/soft/fftw/f95i/include"
        LDFLAGS="-L/afs/ipp/@sys/soft/fftw/f95i/lib -L/afs/ipp/@sys/soft/gnu/lib"
        known_fftw=2
])

AC_DEFUN([ES_KNOWN_DINO],[
	AC_MSG_NOTICE([using preconfigured settings for the Alphas at MPIP])
        known_mpi=dmpi
	CC=cc
	known_fftw=2
])

