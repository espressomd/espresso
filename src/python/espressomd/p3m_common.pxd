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

IF P3M == 1 or DP3M == 1:
    cdef extern from "electrostatics_magnetostatics/p3m-common.hpp":
        ctypedef struct P3MParameters:
            double alpha_L
            double r_cut_iL
            int    mesh[3]
            double mesh_off[3]
            int    cao
            double accuracy
            double epsilon
            double cao_cut[3]
            double a[3]
            double alpha
            double r_cut
            double additional_mesh[3]
