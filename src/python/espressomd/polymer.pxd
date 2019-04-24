#
# Copyright (C) 2013-2018 The ESPResSo project
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


from __future__ import print_function, absolute_import
include "myconfig.pxi"

from libcpp.vector cimport vector
from .utils cimport Vector3d
from .c_analyze cimport PartCfg, partCfg

cdef extern from "polymer.hpp":
    vector[vector[Vector3d]] draw_polymer_positions(PartCfg &, int n_polymers, int beads_per_polymer, double bond_length, vector[Vector3d] & start_positions, double min_distance, int max_tries, int use_bond_angle, double bond_angle, int respect_constraints, int seed) except +
