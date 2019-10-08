# Copyright (C) 2010-2019 The ESPResSo project
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
include "myconfig.pxi"
from .analyze cimport PartCfg, partCfg

cdef extern from "diamond.hpp":
    int create_diamond(PartCfg &, double a, double bond_length, int MPC, int N_CI, double val_nodes, double val_cM, double val_CI, int cM_dist, int nonet)


cdef class Diamond:
    cdef _params
