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
# All Tcl-scripts in the scripts-directory should be listed here

# make a copy of all script files
# PLEASE KEEP ALPHABETIC ORDER
ESSCRIPTS = \
	init.tcl \
	ABHmath.tcl \
	auxiliary.tcl \
	blockfile_support.tcl \
	bundle.tcl \
	countBonds.tcl \
	dielectrics.tcl \
	gibbs_ghmc.tcl \
	parallel_tempering.tcl \
	pdb.tcl \
	polymer.tcl \
	statistics.tcl \
	vtf.tcl \
	vtk.tcl \
	object_in_fluid.tcl

EXTRA_DIST=$(ESSCRIPTS)

# copy the scripts files from the source dir to the obj dir if they are different
# using INSTALL_DATA here as the scripts should not be executable
.PHONY: copy-scripts
FORCE:
copy-scripts: FORCE
	for file in $(ESSCRIPTS); do \
	  if test ! -e $$file || test ! $$file -ef $(abs_srcdir)/$$file; then \
	    $(INSTALL_DATA) $(abs_srcdir)/$$file $(abs_builddir); \
	  fi; \
	done

all-local: copy-scripts
check-local: copy-scripts

clean-local:
	for file in $(ESSCRIPTS); do \
	  if test -e $$file && test ! $$file -ef $(abs_srcdir)/$$file; then \
	    rm -f $$file; \
	  fi; \
	done

# install the scripts from the obj dir to the install dir
# using INSTALL_DATA here as the scripts should not be executable
install-exec-local:
	test -z "$(DESTDIR)$(pkgdatadir)" || $(MKDIR_P) "$(DESTDIR)$(scriptsdir)"
	echo " $(INSTALL_DATA) $(ESSCRIPTS) '$(DESTDIR)$(scriptsdir)'"
	$(INSTALL_DATA) $(ESSCRIPTS) "$(DESTDIR)$(scriptsdir)"

uninstall-local:
	-cd $(DESTDIR)$(scriptsdir) && rm $(ESSCRIPTS)
