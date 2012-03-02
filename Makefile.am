# Copyright (C) 2012 The ESPResSo project
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

# List files that should go into the distribution but are not required
# by any other means
EXTRA_DIST = bootstrap.sh samples packages $(extra)

all-local: Espresso myconfig-sample.h
CLEANFILES = Espesso myconfig-sample.h 

Espresso: 
	$(LN_S) src/Espresso .

myconfig-sample.h: config/myconfig-sample-header.h config/myconfig-sample-body.h
	cat config/myconfig-sample-header.h > myconfig-sample.h
	cat config/myconfig-sample-body.h >> myconfig-sample.h

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
	fi; \
	echo $$VERSION > $(top_distdir)/version.txt
DISCLEANFILES = version.txt

#################################################################
# Documentation
#################################################################
.PHONY: FORCE doc ug dg tutorials doxygen 
ug dg tutorials: FORCE
	cd doc; $(MAKE) --print-directory $@
doc doxygen: FORCE
	cd src; $(MAKE) --print-directory $@
	cd doc; $(MAKE) --print-directory $@
FORCE:
