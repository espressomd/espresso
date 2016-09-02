#
# Copyright (C) 2013,2014,2015,2016 The ESPResSo project
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
# Minimize Energy

from __future__ import print_function, absolute_import
include "myconfig.pxi"
from espressomd._system cimport *
cimport numpy as np
from espressomd.utils cimport *

cdef extern from "minimize_energy.hpp": 
    cdef void minimize_energy_init(const double f_max, const double gamma, const int max_steps, const double max_displacement) 

cdef extern from "communication.hpp":
    cdef int mpi_minimize_energy();
