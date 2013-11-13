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

AM_DEFAULT_SOURCE_EXT = .cpp
LIBS =  ../src/libEspresso.la $(PYTHON_EXTRA_LIBS) @LIBS@

if PYTHON_INTERFACE

myconfig.pxi: ../src/myconfig-final.hpp
	$(PYTHON) $(srcdir)/make_myconfig.py ../src/myconfig-final.hpp

AM_CPPFLAGS = $(PYTHON_CPPFLAGS) -I$(top_srcdir)/src
AM_LDFLAGS = $(PYTHON_LDFLAGS) -module -avoid-version -no-undefined

#SUFFIXES = .pyx
.pyx.cpp:
	cython -I $(builddir) $< -o $@

espresso.cpp: espresso.pxd myconfig.pxi
particle_data.cpp: particle_data.pxd myconfig.pxi
analyze.cpp: analyze.pxd myconfig.pxi
cellsystem.cpp: cellsystem.pxd myconfig.pxi
changeVolume.cpp: changeVolume.pxd myconfig.pxi
code_info.cpp: code_info.pxd myconfig.pxi
debye_hueckel.cpp: debye_hueckel.pxd myconfig.pxi
global_variables.cpp: global_variables.pxd myconfig.pxi
integrate.cpp: integrate.pxd myconfig.pxi
interaction_data.cpp: interaction_data.pxd myconfig.pxi
invalidateSystem.cpp: invalidateSystem.pxd myconfig.pxi
thermostat.cpp: thermostat.pxd myconfig.pxi
utils.cpp: utils.pxd myconfig.pxi

BUILT_SOURCES = 
pyexec_LTLIBRARIES = \
	espresso.la \
	particle_data.la \
	analyze.la \
	cellsystem.la \
	changeVolume.la \
	code_info.la \
	debye_hueckel.la \
	global_variables.la \
	integrate.la \
	interaction_data.la \
	invalidateSystem.la \
	thermostat.la \
	utils.la
#	
#	lb.la
#	cuda_init.la

endif
