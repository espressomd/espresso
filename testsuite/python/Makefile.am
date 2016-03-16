# Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
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
# alphabetically sorted list of test scripts
tests = \
	bondedInteractions.py \
	cellsystem.py \
	electrostaticInteractions.py \
	ewald_gpu.py \
	magnetostaticInteractions.py \
	nonBondedInteractions.py \
	p3m_gpu.py \
	particle.py

# please keep the alphabetic ordering of the above list!

deactivated_tests =

# add data files for the tests here
EXTRA_DIST = $(tests) $(deactivated_tests)

# run the testsuite
# run the testsuite
check-local: runtest.sh
	if test -n "$(processors)"; then \
	  @builddir@/runtest.sh -p "${processors}" $(tests); \
	else \
	  @builddir@/runtest.sh $(tests); \
	fi
