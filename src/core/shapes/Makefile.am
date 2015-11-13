if MPI_FAKE
AM_CPPFLAGS = -I$(srcdir)/../mpifake
else
AM_CPPFLAGS=
endif

AM_CPPFLAGS += -I$(srcdir)/.. -I$(builddir)/..

noinst_LTLIBRARIES = libEspressoShapes.la

libEspressoShapes_la_SOURCES = \
	Shape.hpp Shape.cpp \
	Cylinder.hpp Cylinder.cpp \
	HollowCone.hpp HollowCone.cpp \
	Maze.hpp Maze.cpp \
	Plane.hpp Plane.cpp \
	Pore.hpp Pore.cpp \
	Rhomboid.hpp Rhomboid.cpp \
	Slitpore.hpp Slitpore.cpp \
	Sphere.hpp Sphere.cpp \
	SpheroCylinder.hpp SpheroCylinder.cpp \
	Stomatocyte.hpp Stomatocyte.cpp \
	Wall.hpp Wall.cpp
