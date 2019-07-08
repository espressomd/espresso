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
# Handling of electrostatics

from __future__ import print_function, absolute_import
include "myconfig.pxi"
cimport numpy as np
from libcpp cimport bool

from core.SystemInterface cimport SystemInterface
from espressomd.utils import is_valid_type, to_str
from espressomd.utils cimport handle_errors
from core.p3m_common cimport P3MParameters
from core.p3m cimport p3m_set_mesh_offset, p3m_adaptive_tune, p3m_set_params, p3m_set_tune_params
from core.mmm1d cimport MMM1D_init, MMM1D_sanity_checks, mmm1d_tune


IF ELECTROSTATICS:
    IF P3M:
        IF CUDA:
            from core.p3m_gpu cimport p3m_gpu_init
            cdef inline python_p3m_gpu_init(params):
                cdef int cao
                cdef int mesh[3]
                cdef double alpha
                cao = params["cao"]
                # Mesh can be specified as single int, but here, an array is
                # needed
                if not hasattr(params["mesh"], "__getitem__"):
                    for i in range(3):
                        mesh[i] = params["mesh"]
                else:
                    mesh = params["mesh"]
                alpha = params["alpha"]
                p3m_gpu_init(cao, mesh, alpha)

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
            response = p3m_adaptive_tune( & log)
            handle_errors("Error in p3m_adaptive_tune")
            if log.strip():
                print(to_str(log))
            return response

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
            if is_valid_type(p_mesh, int):
                mesh[0] = p_mesh
                mesh[1] = p_mesh
                mesh[2] = p_mesh
            else:
                mesh = p_mesh

            return p3m_set_params(r_cut, mesh, cao, alpha, accuracy)

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
            if is_valid_type(p_mesh, int):
                mesh[0] = p_mesh
                mesh[1] = p_mesh
                mesh[2] = p_mesh
            else:
                mesh = p_mesh

            p3m_set_tune_params(r_cut, mesh, cao, alpha, accuracy, n_interpol)


IF ELECTROSTATICS:
    cdef inline pyMMM1D_tune():
        cdef char * log = NULL
        cdef int resp
        MMM1D_init();
        if MMM1D_sanity_checks() == 1:
            handle_errors(
                "MMM1D Sanity check failed: wrong periodicity or wrong cellsystem, PRTFM")
        resp = mmm1d_tune(& log)
        if resp:
            print(to_str(log))
        return resp
