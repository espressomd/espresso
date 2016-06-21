if MPI_FAKE
AM_CPPFLAGS = -I$(srcdir)/../mpifake
else
AM_CPPFLAGS=
endif

AM_CPPFLAGS += -I$(srcdir)/.. -I$(builddir)/..

noinst_LTLIBRARIES = libEspressoConstraints.la

libEspressoConstraints_la_SOURCES = \
	ChargedPlate.cpp ChargedPlate.hpp \
	ChargedRod.cpp ChargedRod.hpp \
	Constraint.cpp Constraint.hpp \
	ConstraintList.cpp ConstraintList.hpp \
	ExternalMagneticField.cpp ExternalMagneticField.hpp \
	GeometryConstraint.cpp GeometryConstraint.hpp
