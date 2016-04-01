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

include "myconfig.pxi"
from _system cimport *
cimport numpy as np
from utils cimport *

IF ELECTROSTATICS:
    IF P3M:
        from p3m_common cimport p3m_parameter_struct
    cdef extern from "interaction_data.hpp":
        cdef enum CoulombMethod:
            COULOMB_NONE, \
                COULOMB_DH, \
                COULOMB_P3M, \
                COULOMB_MMM1D, \
                COULOMB_MMM2D, \
                COULOMB_MAGGS, \
                COULOMB_ELC_P3M, \
                COULOMB_RF, \
                COULOMB_INTER_RF, \
                COULOMB_P3M_GPU, \
                COULOMB_MMM1D_GPU, \
                COULOMB_EWALD_GPU, \
                COULOMB_EK

        int coulomb_set_bjerrum(double bjerrum)

        ctypedef struct Coulomb_parameters:
            double bjerrum
            double prefactor
            CoulombMethod method

        cdef extern Coulomb_parameters coulomb

    IF P3M:
        cdef extern from "p3m-common.hpp":
            ctypedef struct p3m_parameter_struct:
                double alpha_L
                double r_cut_iL
                int    mesh[3]
                double mesh_off[3]
                int    cao
                int    inter
                double accuracy
                double epsilon
                double cao_cut[3]
                double a[3]
                double ai[3]
                double alpha
                double r_cut
                int    inter2
                int    cao3
                double additional_mesh[3]

        cdef extern from "p3m.hpp":
            int p3m_set_params(double r_cut, int * mesh, int cao, double alpha, double accuracy)
            void p3m_set_tune_params(double r_cut, int mesh[3], int cao, double alpha, double accuracy, int n_interpol)
            int p3m_set_mesh_offset(double x, double y, double z)
            int p3m_set_eps(double eps)
            int p3m_set_ninterpol(int n)
            int p3m_adaptive_tune(char ** log)

            ctypedef struct p3m_data_struct:
                p3m_parameter_struct params

            # links intern C-struct with python object
            cdef extern p3m_data_struct p3m

        IF CUDA:
            cdef extern from "p3m_gpu.hpp":
                void p3m_gpu_init(int cao, int * mesh, double alpha, double * box)

            cdef inline python_p3m_gpu_init(params):
                cdef int cao
                cdef int mesh[3]
                cdef double alpha
                cdef double box[3]
                cao = params["cao"]
                mesh = params["mesh"]
                alpha = params["alpha"]
                box = params["box"]
                p3m_gpu_init(cao, mesh, alpha, box)

        # Convert C arguments into numpy array
        cdef inline python_p3m_set_mesh_offset(mesh_off):
            cdef double mesh_offset[3]
            mesh_offset[0] = mesh_off[0]
            mesh_offset[1] = mesh_off[1]
            mesh_offset[2] = mesh_off[2]
            return p3m_set_mesh_offset(mesh_offset[0], mesh_offset[1], mesh_offset[2])

        cdef inline python_p3m_adaptive_tune():
            cdef char * log = NULL
            cdef int response
            response = p3m_adaptive_tune(& log)
            return response, log

        cdef inline python_p3m_set_params(p_r_cut, p_mesh, p_cao, p_alpha, p_accuracy):
            cdef int mesh[3]
            cdef double r_cut
            cdef int cao
            cdef double alpha
            cdef double accuracy
            r_cut = p_r_cut
            cao = p_cao
            alpha = p_alpha
            accuracy = p_accuracy
            if isinstance(p_mesh, int):
                mesh[0] = p_mesh
                mesh[1] = p_mesh
                mesh[2] = p_mesh
            else:
                mesh = p_mesh

            p3m_set_params(r_cut, mesh, cao, alpha, accuracy)

        cdef inline python_p3m_set_tune_params(p_r_cut, p_mesh, p_cao, p_alpha, p_accuracy, p_n_interpol):
            # cdef inline python_p3m_set_tune_params():
            cdef int mesh[3]
            cdef double r_cut
            cdef int cao
            cdef double alpha
            cdef double accuracy
            cdef int n_interpol
            r_cut = p_r_cut
            cao = p_cao
            alpha = p_alpha
            accuracy = p_accuracy
            n_interpol = p_n_interpol
            if isinstance(p_mesh, int):
                mesh[0] = p_mesh
                mesh[1] = p_mesh
                mesh[2] = p_mesh
            else:
                mesh = p_mesh

            p3m_set_tune_params(r_cut, mesh, cao, alpha, accuracy, n_interpol)

    cdef extern from "debye_hueckel.hpp":
        IF COULOMB_DEBYE_HUECKEL:
            ctypedef struct Debye_hueckel_params:
                double r_cut
                double kappa
                double eps_int, eps_ext
                double r0, r1
                double alpha
        ELSE:
            ctypedef struct Debye_hueckel_params:
                double r_cut
                double kappa

        cdef extern Debye_hueckel_params dh_params

        int dh_set_params(double kappa, double r_cut)
        int dh_set_params_cdh(double kappa, double r_cut, double eps_int, double eps_ext, double r0, double r1, double alpha)

IF ELECTROSTATICS and CUDA and EWALD_GPU:
    cdef extern from "SystemInterface.hpp":
        cdef cppclass SystemInterface:
            SystemInterface()

    cdef extern from "EspressoSystemInterface.hpp":
        cdef cppclass EspressoSystemInterface:
            @staticmethod
            EspressoSystemInterface * _Instance()

    cdef extern from "actor/EwaldGPU.hpp":
        cdef cppclass EwaldgpuForce:
            EwaldgpuForce(EspressoSystemInterface & s, double r_cut, int num_kx, int num_ky, int num_kz, double alpha)
            int set_params(double rcut, int num_kx, int num_ky, int num_kz, double alpha)
            int set_params_tune(double accuracy, double precision, int K_max, int time_calc_steps)
            int adaptive_tune(char ** log, EspressoSystemInterface & s)
            double tune_alpha(double accuracy, double precision, int K, double V, double q_sqr, int N)
            double tune_rcut(double accuracy, double precision, double alpha, double V, double q_sqr, int N)
            int determine_calc_time_steps()

        ctypedef struct Ewaldgpu_params:
            double rcut
            int num_kx
            int num_ky
            int num_kz
            double alpha
            double accuracy
            double precision
            bint isTuned  # Tuning is over
            bint isTunedFlag  # Flag tuning is over
            int K_max  # Maximal reciprocal K-vector in tuning
            int time_calc_steps  # Steps in time_force_calc function

        cdef extern Ewaldgpu_params ewaldgpu_params

        # ctypedef extern class EwaldgpuForce ewaldgpuForce
#    cdef extern from "EspressoSystemInterface.cpp":
#        cdef cppclass extern EspressoSystemInterface *EspressoSystemInterface;
