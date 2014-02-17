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

SUBDIRS=01-lennard_jones 02-charged_system 03-object_in_fluid

if DEVEL_SRC

.PHONY: FORCE doc tutorials $(SUBDIRS)
doc: tutorials
tutorials: FORCE
	for dir in $(SUBDIRS); do \
		cd $$dir; $(MAKE) --print-directory $@; cd ..;\
	done
FORCE:

endif
