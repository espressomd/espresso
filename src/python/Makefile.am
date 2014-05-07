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

AM_DEFAULT_SOURCE_EXT = .cpp

# Build rules for python libraries
AM_CPPFLAGS = $(PYTHON_CPPFLAGS) -I$(top_srcdir)/src 
AM_LDFLAGS = $(PYTHON_LDFLAGS) -module -avoid-version -no-undefined
LIBS =  ../libEspresso.la @LIBS@
if MPI_FAKE
AM_CPPFLAGS += -I$(top_srcdir)/src/mpifake
endif

# Build rules for cython
SUFFIXES = .pyx

PYXOBJECTS =  \
	espresso \
	particle_data \
	analyze \
	cellsystem \
	changeVolume \
	code_info \
	cuda_init \
	debye_hueckel \
	global_variables \
	integrate \
	interaction_data \
	invalidateSystem \
	thermostat \
	utils

# lb has not been updated yet

CPPOBJECTS = $(PYXOBJECTS:=.cpp)

cython_verbose = $(cython_verbose_@AM_V@)
cython_verbose_ = $(cython_verbose_@AM_DEFAULT_V@)
cython_verbose_0 = @echo "  CYTHON   $@";

$(CPPOBJECTS): %.cpp: %.pyx %.pxd myconfig.pxi
	$(cython_verbose)cython -I $(builddir) $< -o $@

CLEANFILES = $(CPPOBJECTS)

LAOBJECTS = $(PYXOBJECTS:=.la)
pyexec_LTLIBRARIES = $(LAOBJECTS)

# build rule for myconfig.pxi
noinst_PROGRAMS = gen_pxiconfig
gen_pxiconfig_SOURCES = gen_pxiconfig.cpp

gen_pxiconfig.cpp: $(top_srcdir)/config/gen_pxiconfig.py \
		   $(top_srcdir)/src/features.def ../myconfig-final.hpp
	$(AM_V_GEN)$(PYTHON) $(top_srcdir)/config/gen_pxiconfig.py \
	  $(top_srcdir)/src/features.def gen_pxiconfig.cpp

myconfig.pxi: gen_pxiconfig
	$(AM_V_GEN)./gen_pxiconfig > myconfig.pxi

CLEANFILES = gen_pxiconfig.cpp myconfig.pxi
BUILT_SOURCES = gen_pxiconfig.cpp myconfig.pxi

