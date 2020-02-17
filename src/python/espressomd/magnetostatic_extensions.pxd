#
# Copyright (C) 2013-2019 The ESPResSo project
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

include "myconfig.pxi"

IF DIPOLES and DP3M:

    cdef extern from "electrostatics_magnetostatics/mdlc_correction.hpp":
        ctypedef struct dlc_struct "DLC_struct":
            double maxPWerror
            double gap_size
            double far_cut

        int mdlc_set_params(double maxPWerror, double gap_size, double far_cut)

        # links intern C-struct with python object
        cdef extern dlc_struct dlc_params
