# Copyright (C) 2007,2008,2009,2010 Olaf Lenz, Axel Arnold
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
dist_ESSCRIPTS = \
	init.tcl \
	ABHmath.tcl \
	auxiliary.tcl \
	blockfile_support.tcl \
	bundle.tcl \
	countBonds.tcl \
	pdb.tcl \
	polymer.tcl \
	statistics.tcl \
	vtf.tcl \
	tclline.tcl \
	parallel_tempering.tcl

# copy the scripts files from the source dir to the obj dir
# using INSTALL_DATA here as the scripts should not be executable
all:
	cd $(srcdir) && $(INSTALL_DATA) $(dist_ESSCRIPTS) $(abs_builddir)

# install the scripts from the obj dir to the install dir
# using INSTALL_DATA here as the scripts should not be executable
install-exec-local:
	test -z "$(DESTDIR)$(pkgdatadir)" || $(MKDIR_P) "$(DESTDIR)$(scriptsdir)"
	echo " $(INSTALL_DATA) $(dist_ESSCRIPTS) '$(DESTDIR)$(scriptsdir)'"
	$(INSTALL_DATA) $(dist_ESSCRIPTS) "$(DESTDIR)$(scriptsdir)"

uninstall-local:
	-cd $(DESTDIR)$(scriptsdir) && rm $(dist_ESSCRIPTS)
