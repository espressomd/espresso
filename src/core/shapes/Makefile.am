AM_CPPFLAGS = -I$(top_srcdir)/src/core
if MPI_FAKE
AM_CPPFLAGS += -I$(top_srcdir)/src/core/mpifake
endif

noinst_LTLIBRARIES = libEspressoShapes.la

libEspressoShapes_la_SOURCES = \
	Shape.hpp \
	ShapeList.cpp ShapeList.hpp \
	Wall.hpp Wall.cpp

