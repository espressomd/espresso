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
# Handling of electrostatics

from __future__ import print_function, absolute_import
include "myconfig.pxi"
from espressomd.system cimport *
from espressomd.utils cimport *
from espressomd.electrostatics cimport *

IF ELECTROSTATICS and P3M:

    cdef extern from "elc.hpp":
        ctypedef struct ELC_struct:
            double maxPWerror
            double gap_size
            double far_cut
            int neutralize
            double delta_mid_top,
            double delta_mid_bot,
            int const_pot,
            double pot_diff

        int ELC_set_params(double maxPWerror, double min_dist, double far_cut,
        int neutralize, double delta_mid_top, double delta_mid_bot, int const_pot, double pot_diff)

        # links intern C-struct with python object
        ELC_struct elc_params

    cdef extern from "iccp3m.hpp":
        ctypedef struct iccp3m_struct:
            int n_ic
            int num_iteration
            double eout
            double * areas
            double * ein
            double * sigma
            double convergence
            double * nvectorx
            double * nvectory
            double * nvectorz
            double extx
            double exty
            double extz
            double relax
            int citeration
            int set_flag
            double * fx
            double * fy
            double * fz
            int first_id

        # links intern C-struct with python object
        iccp3m_struct iccp3m_cfg

        void iccp3m_set_initialized()
        void iccp3m_alloc_lists()

    cdef extern from "communication.hpp":
        int mpi_iccp3m_init(int dummy)
