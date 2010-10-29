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
SUBDIRS = dg ug tutorials maggs_doc

.PHONY: FORCE doc dg ug tutorials
doc: ug dg tutorials
dg: FORCE
	cd dg; $(MAKE) --print-directory $@
ug: FORCE
	cd ug; $(MAKE) --print-directory $@
tutorials: FORCE
	cd tutorials; $(MAKE) --print-directory $@
FORCE:

