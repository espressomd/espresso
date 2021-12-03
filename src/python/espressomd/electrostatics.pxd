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
from libcpp cimport bool
from .analyze cimport PartCfg, partCfg

cdef extern from "SystemInterface.hpp":
    cdef cppclass SystemInterface:
        pass
cdef extern from "EspressoSystemInterface.hpp":
    cdef cppclass EspressoSystemInterface(SystemInterface):
        @staticmethod
        EspressoSystemInterface & Instance()


IF ELECTROSTATICS:
    cdef extern from "electrostatics_magnetostatics/common.hpp":
        void mpi_bcast_coulomb_params()

    IF P3M:
        from p3m_common cimport P3MParameters

    cdef extern from "electrostatics_magnetostatics/coulomb.hpp":

        cdef enum CoulombMethod:
            COULOMB_NONE, \
                COULOMB_DH, \
                COULOMB_P3M, \
                COULOMB_MMM1D, \
                COULOMB_ELC_P3M, \
                COULOMB_RF, \
                COULOMB_P3M_GPU, \
                COULOMB_MMM1D_GPU, \
                COULOMB_SCAFACOS

        ctypedef struct Coulomb_parameters:
            double prefactor
            CoulombMethod method

        cdef extern Coulomb_parameters coulomb

    cdef extern from "electrostatics_magnetostatics/coulomb.hpp" namespace "Coulomb":

        int set_prefactor(double prefactor) except +
        void deactivate_method()

    IF P3M:
        from p3m_common cimport P3MParameters

        cdef extern from "electrostatics_magnetostatics/p3m.hpp":
            void p3m_set_params(double r_cut, int * mesh, int cao, double alpha, double accuracy) except +
            void p3m_set_tune_params(double r_cut, int mesh[3], int cao, double accuracy)
            void p3m_set_mesh_offset(double x, double y, double z) except +
            void p3m_set_eps(double eps)
            int p3m_adaptive_tune(int timings, bool verbose)

            ctypedef struct p3m_data_struct:
                P3MParameters params

            # links intern C-struct with python object
            cdef extern p3m_data_struct p3m

        IF CUDA:
            cdef extern from "electrostatics_magnetostatics/p3m_gpu.hpp":
                void p3m_gpu_init(int cao, int * mesh, double alpha) except +

        cdef extern from "electrostatics_magnetostatics/elc.hpp":
            ctypedef struct ELCParameters:
                double maxPWerror
                double gap_size
                double far_cut
                bool neutralize
                double delta_mid_top
                double delta_mid_bot
                bool const_pot
                double pot_diff

            void ELC_set_params(double maxPWerror, double min_dist, double far_cut,
                                bool neutralize, double delta_mid_top,
                                double delta_mid_bot, bool const_pot, double pot_diff) except +

            # links intern C-struct with python object
            ELCParameters elc_params

    cdef extern from "electrostatics_magnetostatics/debye_hueckel.hpp":
        ctypedef struct DebyeHueckelParameters:
            double r_cut
            double kappa

        cdef extern DebyeHueckelParameters dh_params

        void dh_set_params(double kappa, double r_cut) except +

    cdef extern from "electrostatics_magnetostatics/reaction_field.hpp":
        ctypedef struct ReactionFieldParameters:
            double kappa
            double epsilon1
            double epsilon2
            double r_cut

        cdef extern ReactionFieldParameters rf_params

        void rf_set_params(double kappa, double epsilon1, double epsilon2,
                           double r_cut) except +

IF ELECTROSTATICS:
    cdef extern from "electrostatics_magnetostatics/mmm1d.hpp":
        ctypedef struct MMM1DParameters:
            double far_switch_radius_2
            double maxPWerror
            int    bessel_cutoff

        cdef extern MMM1DParameters mmm1d_params

        void MMM1D_set_params(double switch_rad, double maxPWerror)
        int MMM1D_init()
        int mmm1d_tune(int timings, bool verbose)

IF ELECTROSTATICS and MMM1D_GPU:

    cdef extern from "actor/Mmm1dgpuForce.hpp":
        cdef cppclass Mmm1dgpuForce:
            Mmm1dgpuForce(SystemInterface & s) except+
            void setup(SystemInterface & s)
            void tune(SystemInterface & s, float _maxPWerror, float _far_switch_radius, int _bessel_cutoff)
            void set_params(float _boxz, float _coulomb_prefactor, float _maxPWerror, float _far_switch_radius, int _bessel_cutoff)

            void activate()
            void deactivate()

cdef extern from "utils/checks/charge_neutrality.hpp" namespace "Utils":
    bool check_charge_neutrality[ParticleRange](ParticleRange & partCfg)
