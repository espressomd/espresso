# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
	-I$(top_srcdir)/src -I$(top_srcdir)/src/core -I$(top_builddir)/src/core \
	-I$(top_srcdir)/src/core/actor -I$(top_builddir)/scr/core/actor
AM_LDFLAGS = $(PYTHON_LDFLAGS) -module -avoid-version -no-undefined

if MPI_FAKE
AM_CPPFLAGS += -I$(top_srcdir)/src/core/mpifake
endif

#CYTHON compile options e.g. debug
AM_CYFLAGS = $(CYFLAGS)

# Build rules for cython
SUFFIXES = .pyx

EXTRA_DIST = gen_pxiconfig.py
PYFILES = __init__.py \
	  system.py \
	  highlander.py

PYXFILES =  \
	_init.pyx \
	_system.pyx \
	actors.pyx \
	analyze.pyx \
	cellsystem.pyx \
	code_info.pyx \
	cuda_init.pyx \
	debye_hueckel.pyx \
	electrokinetics.pyx \
	electrostatic_extensions.pyx \
	electrostatics.pyx \
	integrate.pyx \
	interactions.pyx \
	lb.pyx \
	magnetostatics.pyx \
	magnetostatic_extensions.pyx \
	particle_data.pyx \
	reaction.pyx \
	thermostat.pyx \
	utils.pyx \
	h5md.pyx \
	grand_canonical.pyx

PXDFILES = \
	_espresso/globals.pxd \
	analyze.pxd \
	c_analyze.pxd \
	cellsystem.pxd \
	code_info.pxd \
	cuda_init.pxd \
	debye_hueckel.pxd \
	electrostatic_extensions.pxd \
	electrostatics.pxd \
	_espresso.pxd \
	_init.pxd \
	integrate.pxd \
	interactions.pxd \
	lb.pxd \
	magnetostatics.pxd \
	magnetostatic_extensions.pxd \
	p3m_common.pxd \
	particle_data.pxd \
	reaction.pxd \
	_system.pxd \
	_tcl.pxd \
	thermostat.pxd \
	utils.pxd \
	actors.pxd \
	h5md.pxd \
	grand_canonical.pxd

LIBS +=	$(top_builddir)/src/core/libEspresso.la

if TCL
LIBS +=	$(top_builddir)/src/tcl/libEspressoTcl4Py.la
PYXFILES += _tcl.pyx
endif

# Python files have to be copied to the build dir
python-all:
	$(AM_V_GEN)for f in $(PYFILES); do \
          if test ! -f $$f || test $(srcdir)/$$f -nt $$f; then \
            $(INSTALL) -D $(srcdir)/$$f $$f; \
          fi ; \
        done
all-local: python-all

PYCFILES = $(PYFILES:.py=.pyc)
CLEANFILES = $(PYCFILES)
clean-local: 
	for file in $(PYFILES); do \
          if test -e $$file && test ! $$file -ef $(abs_srcdir)/$$file; then \
            rm -f $$file; \
          fi; \
	done

EXTRA_DIST += $(PYFILES) $(PYXFILES) $(PXDFILES)

# Generate the ".cpp" files 
CPPFILES = $(PYXFILES:.pyx=.cpp)

cython_verbose = $(cython_verbose_@AM_V@)
cython_verbose_ = $(cython_verbose_@AM_DEFAULT_V@)
cython_verbose_0 = @echo "  CYTHON   $@";

$(CPPFILES): %.cpp: %.pyx myconfig.pxi
	$(cython_verbose)cython $(CYFLAGS) --cplus \
	  -I $(builddir) -I $(srcdir)/_espresso $< -o $@

CLEANFILES += $(CPPFILES)

# Dummy target to generate libtool rules
noinst_PROGRAMS = dummy
dummy_SOURCES = dummy.cpp

# Generate the ".la" files. We do not really need these files, but we
# want to use libtool, so we have to give it a corresponding target.
# Also, like this, the .so-files are installed to the right place.
LAFILES = $(PYXFILES:.pyx=.la)
$(LAFILES): %.la: %.lo
	$(AM_V_CXXLD)$(CXXLINK) -rpath $(abs_builddir) -L$(top_builddir)/src/core/.libs $(LIBS) -lEspresso -rpath $(pyexecdir) $<
CLEANFILES += $(LAFILES)

# create links to the ./libs/*.so in the dir
SOFILES = $(PYXFILES:.pyx=.so)
all-local: $(SOFILES)
CLEANFILES += $(SOFILES)

if DARWIN
SOSUFFIX=dylib
else
SOSUFFIX=so
endif

$(SOFILES): %.so: %.la libEspresso.so
	$(AM_V_GEN)test -L $@ || $(LN_S) @objdir@/$@ $@

libEspresso.so:
	$(AM_V_GEN)test -L $@ || $(LN_S) $(top_builddir)/src/core/.libs/libEspresso.$(SOSUFFIX) $@

# build rule for myconfig.pxi
# create gen_pxiconfig.cpp from features.def
gen_pxiconfig.cpp: $(srcdir)/gen_pxiconfig.py $(top_srcdir)/src/features.def 
	$(AM_V_GEN)$(PYTHON) $(srcdir)/gen_pxiconfig.py \
	  $(top_srcdir)/src/features.def gen_pxiconfig.cpp

gen_pxiconfig_CPPFLAGS=-I$(top_srcdir)/src/core -I$(top_builddir)/src/core -I$(top_builddir)/src -DHAVE_CONFIG_H
gen_pxiconfig: gen_pxiconfig.cpp $(top_builddir)/src/core/config.hpp $(top_builddir)/src/core/myconfig-final.hpp
	$(AM_V_CXXLD)$(CXXLINK) $(gen_pxiconfig_CPPFLAGS) $<

CLEANFILES += gen_pxiconfig

# create myconfig.pxi from gen_pxiconfig
myconfig.pxi: gen_pxiconfig
	$(AM_V_GEN)./gen_pxiconfig > myconfig.pxi

CLEANFILES += gen_pxiconfig.cpp myconfig.pxi
BUILT_SOURCES = gen_pxiconfig.cpp myconfig.pxi

