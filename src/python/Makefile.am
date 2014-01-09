# Copyright (C) 2013 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

AM_DEFAULT_SOURCE_EXT = .pyx

# Build rules for python libraries
AM_CPPFLAGS = $(PYTHON_CPPFLAGS) -I$(top_srcdir)/src
AM_LDFLAGS = $(PYTHON_LDFLAGS) $(PYTHON_EXTRA_LDFLAGS) -module -avoid-version -no-undefined
LIBS =  ../libEspresso.la $(PYTHON_EXTRA_LIBS) @LIBS@
if MPI_FAKE
AM_CPPFLAGS += -I$(top_srcdir)/src/mpifake
endif

# build rule for myconfig.pxi
myconfig.pxi: ../myconfig-final.hpp
	$(AM_V_GEN)$(PYTHON) $(srcdir)/make_myconfig.py ../myconfig-final.hpp > /dev/null
CLEANFILES = myconfig.pxi
BUILT_SOURCES = myconfig.pxi

# Build rules for cython
SUFFIXES = .pyx
cython_verbose = $(cython_verbose_@AM_V@)
cython_verbose_ = $(cython_verbose_@AM_DEFAULT_V@)
cython_verbose_0 = @echo "  CYTHON   $@";
.pyx.cpp:
	$(cython_verbose)cython -I $(builddir) $< -o $@

.PRECIOUS: .pyx.cpp

#espresso.cpp: espresso.pxd myconfig.pxi
#particle_data.cpp: particle_data.pxd myconfig.pxi
#analyze.cpp: analyze.pxd myconfig.pxi
# cellsystem.cpp: cellsystem.pxd myconfig.pxi
# changeVolume.cpp: changeVolume.pxd myconfig.pxi
# code_info.cpp: code_info.pxd myconfig.pxi
# debye_hueckel.cpp: debye_hueckel.pxd myconfig.pxi
# global_variables.cpp: global_variables.pxd myconfig.pxi
# integrate.cpp: integrate.pxd myconfig.pxi
# interaction_data.cpp: interaction_data.pxd myconfig.pxi
# invalidateSystem.cpp: invalidateSystem.pxd myconfig.pxi
# thermostat.cpp: thermostat.pxd myconfig.pxi
# utils.cpp: utils.pxd myconfig.pxi

CLEANFILES = espresso.cpp \
	particle_data.cpp \
	analyze.cpp \
	cellsystem.cpp \
	changeVolume.cpp \
	code_info.cpp \
	debye_hueckel.cpp \
	global_variables.cpp \
	integrate.cpp \
	interaction_data.cpp \
	invalidateSystem.cpp \
	thermostat.cpp \
	utils.cpp

pyexec_LTLIBRARIES = \
	espresso.la \
	particle_data.la \
	analyze.la \
	cellsystem.la \
	changeVolume.la \
	code_info.la \
	debye_hueckel.la \
	global_variables.la \
	integrate.la \
	interaction_data.la \
	invalidateSystem.la \
	thermostat.la \
	utils.la

#	
#	lb.la
#	cuda_init.la
