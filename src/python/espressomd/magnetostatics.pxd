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

from libcpp cimport bool

include "myconfig.pxi"

cdef extern from "SystemInterface.hpp":
    cdef cppclass SystemInterface:
        pass
cdef extern from "EspressoSystemInterface.hpp":
    cdef cppclass EspressoSystemInterface(SystemInterface):
        @staticmethod
        EspressoSystemInterface & Instance()

IF DIPOLES == 1:
    cdef extern from "electrostatics_magnetostatics/common.hpp":
        void mpi_bcast_coulomb_params()

    cdef extern from "electrostatics_magnetostatics/dipole.hpp" namespace "Dipole":
        void set_Dprefactor(double prefactor) except +
        double get_Dprefactor()
        void disable_method_local()

    cdef extern from "electrostatics_magnetostatics/magnetic_non_p3m_methods.hpp":
        void dawaanr_set_params() except +
        void mdds_set_params(int n_replica) except +
        int mdds_get_n_replica()

    IF(CUDA == 1) and (ROTATION == 1):
        cdef extern from "actor/DipolarDirectSum.hpp":
            cdef cppclass DipolarDirectSum:
                DipolarDirectSum(SystemInterface & s) except +
                void set_params()
                void activate()
                void deactivate()

    IF(DIPOLAR_BARNES_HUT == 1):
        cdef extern from "actor/DipolarBarnesHut.hpp":
            cdef cppclass DipolarBarnesHut:
                DipolarBarnesHut(SystemInterface & s) except +
                void set_params(float epssq, float itolsq)
                void activate()
                void deactivate()

IF DP3M == 1:
    from p3m_common cimport P3MParameters

    cdef extern from "electrostatics_magnetostatics/p3m-dipolar.hpp":
        void dp3m_set_params(double r_cut, int mesh, int cao, double alpha, double accuracy) except +
        void dp3m_set_tune_params(double r_cut, int mesh, int cao, double accuracy)
        void dp3m_set_mesh_offset(double x, double y, double z) except +
        void dp3m_set_eps(double eps)
        int dp3m_adaptive_tune(int timings, bool verbose)
        void dp3m_deactivate()

        ctypedef struct dp3m_data_struct:
            P3MParameters params

        cdef extern dp3m_data_struct dp3m
