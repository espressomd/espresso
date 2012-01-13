# Copyright (C) 2011 Stefan Kesselheim, Olaf Lenz
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

#CYTHON =cython
#CYTHONARG =

all: pythonlib

.PHONY: pythonlib

BUILT_SOURCES = $(abs_srcdir)/myconfig.pxi 

$(abs_srcdir)/myconfig.pxi: myconfig.pxi
	cp myconfig.pxi $(abs_srcdir)/myconfig.pxi

myconfig.pxi: $(abs_top_builddir)/src/myconfig-final.h
	$(PYTHON) $(abs_top_srcdir)/cython/make_myconfig.py $(abs_top_builddir)/src/myconfig-final.h

PYTHON_SETUP=CC='$(CC) $(CPPFLAGS)' $(PYTHON) setup.py

pythonlib:
	cd $(abs_srcdir) && \
	  $(PYTHON_SETUP) build --build-base=$(abs_builddir) build_ext -I$(abs_top_srcdir)/src:$(abs_top_builddir)/src -L$(abs_top_builddir)/src/.libs

install-exec-local: 
	cd $(abs_srcdir) && \
	  $(PYTHON_SETUP) \
	    build --build-base=$(abs_builddir) \
	    build_ext -I$(abs_top_srcdir)/src:$(abs_top_builddir)/src -L$(abs_top_builddir)/src/.libs \
	    install --install-lib=$(prefix)/lib/ 

clean-local:
	cd $(abs_srcdir) && \
	  $(PYTHON_SETUP) clean

#pyexec_LTLIBRARIES = libespresso.la libinteraction_data.la libglobal_variables.la libparticle_data.la
#
#libespresso_la_SOURCES = espresso.pyx 
#libinteraction_data_la_SOURCES = interaction_data.pyx 
#libglobal_variables_la_SOURCES = global_variables.pyx 
#libparticle_data_la_SOURCES = particle_data.pyx
#
#
#
#.pyx.c:
#	${CYTHON} ${CYTHONARGS} $< -o $@
#c.lo:
#	${LIBTOOL} --mode=compile ${CC} ${CPPFLAGS} $< -o $@ 
#

