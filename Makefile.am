# Copyright (C) 2012,2013 The ESPResSo project
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
SUBDIRS = config src . scripts tools testsuite doc

ACLOCAL_AMFLAGS = -I config
AUTOMAKE_OPTIONS = gnu subdir-objects dist-xz

# List files that should go into the distribution but are not required
# by any other means
EXTRA_DIST = bootstrap.sh samples packages $(extra)

all-local: Espresso
CLEANFILES = Espresso 

if !DEVEL_SRC
Espresso: FORCE
	test -e Espresso || $(LN_S) src/Espresso .
else
Espresso: version.txt
	test -e Espresso || $(LN_S) src/Espresso .
endif

# Sample myconfig generation
myconfig-sample.h: 
	./config.status myconfig-sample.h
CLEANFILES += myconfig-sample.h
# config.status deps in configure.ac

# Handle version
EXTRA_DIST += version.txt

if DEVEL_SRC
MOSTLYCLEANFILES = version.txt
version.txt: FORCE
	cd $(top_srcdir); sh config/genversion.sh -rd > $(abs_top_builddir)/version.txt

dist-hook:
	@VERSION=`cd $(top_srcdir); sh config/genversion.sh -rd`; \
	if test "$$VERSION" != "$(PACKAGE_VERSION)"; then \
		echo "***************************************************************************"; \
		echo "ERROR: Package version is not up-to-date."; \
		echo "ERROR: $$VERSION (current) != $(PACKAGE_VERSION) (configure)"; \
		echo "ERROR: Please run autoreconf -vf"; \
		echo "ERROR: before running 'make dist' again."; \
		echo "***************************************************************************"; \
		exit 1; \
	fi;
endif

#################################################################
# Documentation
#################################################################
if DEVEL_SRC

.PHONY: FORCE doc ug dg tutorials doxygen 
ug dg tutorials: FORCE
	cd doc; $(MAKE) --print-directory $@

else

.PHONY: FORCE doc doxygen 

endif

doc doxygen: FORCE
	cd src; $(MAKE) --print-directory $@
	cd doc; $(MAKE) --print-directory $@

FORCE:
