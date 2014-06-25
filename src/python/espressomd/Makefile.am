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
AM_CPPFLAGS = $(PYTHON_CPPFLAGS) \
	-I$(top_srcdir)/src -I$(top_srcdir)/src/core -I$(top_builddir)/src/core
AM_LDFLAGS = $(PYTHON_LDFLAGS) -module -avoid-version -no-undefined
LIBS =  $(top_builddir)/src/core/libEspresso.la @LIBS@

if MPI_FAKE
AM_CPPFLAGS += -I$(top_srcdir)/src/mpifake
endif

# Build rules for cython
SUFFIXES = .pyx

EXTRA_DIST = gen_pxiconfig.py
PYFILES = __init__.py

PYXFILES =  \
	System.pyx \
	particle_data.pyx \
	analyze.pyx \
	cellsystem.pyx \
	changeVolume.pyx \
	code_info.pyx \
	cuda_init.pyx \
	debye_hueckel.pyx \
	global_variables.pyx \
	integrate.pyx \
	interactions.pyx \
	invalidateSystem.pyx \
	thermostat.pyx \
	utils.pyx

PXDFILES = \
	System.pxd \
	particle_data.pxd \
	analyze.pxd \
	c_analyze.pxd \
	cellsystem.pxd \
	changeVolume.pxd \
	code_info.pxd \
	cuda_init.pxd \
	debye_hueckel.pxd \
	global_variables.pxd \
	integrate.pxd \
	interactions.pxd \
	invalidateSystem.pxd \
	thermostat.pxd \
	utils.pxd
# lb has not been updated yet

# Python files have to be copied to the build dir
python-all:
	$(AM_V_GEN)for f in $(PYFILES); do \
          if test ! -f $$f || test $(srcdir)/$$f -nt $$f; then \
            $(INSTALL) -D $(srcdir)/$$f $$f; \
          fi ; \
        done
all-local: python-all

PYCFILES = $(PYFILES:.py=.pyc)
CLEANFILES = $(PYFILES) $(PYCFILES)

EXTRA_DIST += $(PYFILES) $(PYXFILES) $(PXDFILES)

# Generate the ".cpp" files 
CPPFILES = $(PYXFILES:.pyx=.cpp)

cython_verbose = $(cython_verbose_@AM_V@)
cython_verbose_ = $(cython_verbose_@AM_DEFAULT_V@)
cython_verbose_0 = @echo "  CYTHON   $@";

$(CPPFILES): %.cpp: %.pyx %.pxd myconfig.pxi
	$(cython_verbose)cython -I $(builddir) $< -o $@

CLEANFILES += $(CPPFILES)

# Dummy target to generate libtool rules
noinst_PROGRAMS = dummy
dummy_SOURCES = dummy.cpp

# Generate the ".la" files. We do not really need these files, but we
# want to use libtool, so we have to give it a corresponding target.
# Also, like this, the .so-files are installed to the right place.
LAFILES = $(PYXFILES:.pyx=.la)
$(LAFILES): %.la: %.lo
	$(AM_V_CXXLD)$(CXXLINK) -rpath $(pyexecdir) $< $(LIBS)
CLEANFILES += $(LAFILES)

# create links to the ./libs/*.so in the dir
SOFILES = $(PYXFILES:.pyx=.so)
all-local: $(SOFILES)
CLEANFILES += $(SOFILES)

$(SOFILES): %.so: %.la
	$(AM_V_GEN)test -L $@ || $(LN_S) @objdir@/$@ $@

# build rule for myconfig.pxi
# create gen_pxiconfig.cpp from features.def
gen_pxiconfig.cpp: $(srcdir)/gen_pxiconfig.py $(top_srcdir)/src/features.def 
	$(AM_V_GEN)$(PYTHON) $(srcdir)/gen_pxiconfig.py \
	  $(top_srcdir)/src/features.def gen_pxiconfig.cpp

gen_pxiconfig_CPPFLAGS=-I$(top_srcdir)/src/core -I$(top_builddir)/src/core -I$(top_builddir)/src
gen_pxiconfig: gen_pxiconfig.cpp $(top_builddir)/src/core/config.hpp $(top_builddir)/src/core/myconfig-final.hpp
	$(AM_V_CXXLD)$(CXXLINK) $(gen_pxiconfig_CPPFLAGS) $<
	
CLEANFILES += gen_pxiconfig

# create myconfig.pxi from gen_pxiconfig
myconfig.pxi: gen_pxiconfig
	$(AM_V_GEN)./gen_pxiconfig > myconfig.pxi

CLEANFILES += gen_pxiconfig.cpp myconfig.pxi
BUILT_SOURCES = gen_pxiconfig.cpp myconfig.pxi

