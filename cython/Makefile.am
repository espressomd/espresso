# Copyright (C) 2013 The ESPResSo project
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

.pyx.cpp:
	cython -I $(builddir) $< -o $@

pyexec_LTLIBRARIES = global_variables.la
global_variables_la_SOURCES = global_variables.cpp
global_variables_la_LDFLAGS =  -module -avoid-version
global_variables_la_LIBADD = libEspresso.la

myconfig.pxi: ../src/myconfig-final.hpp
	$(PYTHON) $(srcdir)/make_myconfig.py ../src/myconfig-final.hpp
