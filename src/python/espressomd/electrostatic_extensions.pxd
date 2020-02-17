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
from .electrostatics cimport *
from libcpp.vector cimport vector
from libcpp cimport bool
from .utils cimport Vector3d

IF ELECTROSTATICS and P3M:

    cdef extern from "electrostatics_magnetostatics/elc.hpp":
        ctypedef struct ELC_struct:
            double maxPWerror
            double gap_size
            double far_cut
            bool neutralize
            double delta_mid_top,
            double delta_mid_bot,
            bool const_pot,
            double pot_diff

        int ELC_set_params(double maxPWerror, double min_dist, double far_cut,
                           bool neutralize, double delta_mid_top, double delta_mid_bot, bool const_pot, double pot_diff)

        # links intern C-struct with python object
        ELC_struct elc_params

    cdef extern from "electrostatics_magnetostatics/icc.hpp":
        ctypedef struct iccp3m_struct:
            int n_ic
            int num_iteration
            double eout
            vector[double] areas
            vector[double] ein
            vector[double] sigma
            double convergence
            vector[Vector3d] normals
            Vector3d ext_field
            double relax
            int citeration
            int first_id

        # links intern C-struct with python object
        iccp3m_struct iccp3m_cfg

        void iccp3m_alloc_lists()

    cdef extern from "communication.hpp":
        int mpi_iccp3m_init()
