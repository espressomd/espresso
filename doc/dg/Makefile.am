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
EXTRA_DIST=\
	background_errors.sh background_errors.awk \
	Doxyfile \
        pages/build.doc \
        pages/doxygen.doc \
        pages/header.doc \
        pages/new_files.doc \
        pages/progguide.doc \
        pages/ug.doc \
        figs/bond_angle.fig \
        figs/bond_angle.gif \
        figs/datastorage.gif \
        figs/dihedral_angle.fig \
        figs/dihedral_angle.gif \
        figs/dihedral_angle.pdf \
        figs/directions.fig \
        figs/directions.gif \
        figs/elc_errordist.gif \
        figs/ghost_cells.fig \
        figs/ghost_cells.gif \
        figs/ghost_communication.fig \
        figs/ghost_communication.gif \
        figs/linked_cells.fig \
        figs/linked_cells.gif \
        figs/logo.png \
        figs/move_to_p_buf.fig \
        figs/move_to_p_buf.gif \
        figs/particles.fig \
        figs/particles.gif

.PHONY: doc dg
doc: dg
dg:
if HAVE_DOXYGEN
	SRCDIR=$(srcdir) AWK=$(AWK)\
	  $(SHELL) $(srcdir)/background_errors.sh $(top_srcdir)/src/*.c $(top_srcdir)/src/*.h 
	SRCDIR=$(top_srcdir) PACKAGE="$(PACKAGE_STRING)"\
	  $(DOXYGEN) $(srcdir)/Doxyfile
	@echo "See doxygen's warnings in dg.log."
	@echo "***************************************************************************"
	@echo "* The developer's guide is now accessible at"
	@echo "*  file://`pwd`/html/index.html"
	@echo "***************************************************************************"
else
	@echo "doxygen was not found in your PATH."
	@echo "Can't build developer's guide without doxygen."
	@echo "Install doxygen and rerun configure."

endif

MOSTLYCLEANFILES=background_errors.doc background_errors.sorted background_errors.unsorted
clean-local:
	-rm -rf html
