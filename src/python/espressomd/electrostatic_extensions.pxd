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
from .utils cimport Vector3d

IF ELECTROSTATICS and P3M:

    cdef extern from "electrostatics_magnetostatics/icc.hpp":
        ctypedef struct icc_struct:
            int n_icc
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
        cdef extern icc_struct icc_cfg

        void icc_set_params(int n_icc, double convergence, double relaxation,
                            Vector3d & ext_field, int max_iterations,
                            int first_id, double eps_out,
                            vector[double] & areas,
                            vector[double] & e_in,
                            vector[double] & sigma,
                            vector[Vector3d] & normals) except +

        void icc_deactivate()
