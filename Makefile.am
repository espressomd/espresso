# Copyright (C) 2007,2008,2009,2010,2011 Olaf Lenz, Axel Arnold
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
SUBDIRS = src . scripts tools testsuite doc

ACLOCAL_AMFLAGS = -I config

# List files that should go into the distribution but are not required
# by any other means
EXTRA_DIST = bootstrap.sh samples packages $(extra)

#################################################################
# Installation of the main program
#################################################################
# Two binaries are generated: Espresso for the build dir,
# Espresso.install for the installation dir. 
# The ".install" suffix is removed upon installation.

noinst_PROGRAMS = Espresso
Espresso_CPPFLAGS = -D 'ESPRESSO_SCRIPTS_DEFAULT="$(buildscriptsdir)"'
Espresso_SOURCES = src/scriptsdir.c src/main.c
Espresso_LDADD = src/libEspresso.a

bin_PROGRAMS = Espresso.install
Espresso_install_CPPFLAGS = -D 'ESPRESSO_SCRIPTS_DEFAULT="$(scriptsdir)"'
Espresso_install_SOURCES = src/scriptsdir.c src/main.c
Espresso_install_LDADD = src/libEspresso.a

# rename Espresso after installation
install-exec-hook:
	$(am__mv) \
		$(DESTDIR)$(bindir)/Espresso.install \
		$(DESTDIR)$(bindir)/Espresso

uninstall-local:
	-rm -f $(DESTDIR)$(bindir)/Espresso


#################################################################
# Documentation
#################################################################
.PHONY: FORCE doc doxygen ug dg tutorials
doc doxygen ug dg tutorials: FORCE
	cd doc; $(MAKE) --print-directory $@
FORCE:
