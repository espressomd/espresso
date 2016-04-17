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
include "myconfig.pxi"
import numpy as np
from actors cimport Actor
from globals cimport temperature

IF DIPOLES == 1:
    cdef class MagnetostaticInteraction(Actor):

        def validate_params(self):
            if not (("bjerrum_length" in self._params) ^ ("prefactor" in self._params)):
                raise ValueError(
                    "Either the bjerrum length or the explicit prefactor has to be given")

            if "bjerrum_length" in self._params:
                if not (self._params["bjerrum_length"] > 0.0):
                    raise ValueError(
                        "Bjerrum_length should be a positive double")
            if "prefactor" in self._params:
                if not self._params["prefactor"] > 0:
                    raise ValueError("prefactor should be a positive double")

        def set_magnetostatics_prefactor(self):
            """changes the magnetostatics prefactor, using either the bjrerrum
               length or the explicit prefactor given in the _params dictionary
               of the class."""
            if "bjerrum_length" in self._params:
                if temperature == 0:
                    raise Exception(
                        "Bjerrum length is not defined, if temperature is zero")
                if dipolar_set_Dbjerrum(self._params["bjerrum_length"]):
                    raise Exception(
                        "Could not set magnetostatic bjerrum length")
                return True
            if "prefactor" in self._params:
                if temperature == 0.:
                    if dipolar_set_Dbjerrum(self._params["prefactor"]):
                        raise Exception(
                            "Could not set magnetostatic prefactor")
                else:
                    if dipolar_set_Dbjerrum(self._params["prefactor"] / temperature):
                        raise Exception(
                            "Could not set magnetostatic prefactor")
                    else:
                        self._params["bjerrum_length"] = self.params[
                            "prefactor"] / temperature

        def get_params(self):
            self._params = self._get_params_from_es_core()
            return self._params

        def _get_active_method_from_es_core(self):
            return coulomb.Dmethod


IF DP3M == 1:
    cdef class DipolarP3M(MagnetostaticInteraction):

        def validate_params(self):
            super(DipolarP3M, self).validate_params()
            default_params = self.default_params()

            if not (self._params["r_cut"] >= 0 or self._params["r_cut"] == default_params["r_cut"]):
                raise ValueError("P3M r_cut has to be >=0")

            if not (isinstance(self._params["mesh"], int) or len(self._params["mesh"])):
                raise ValueError(
                    "P3M mesh has to be an integer or integer list of length 3")

            if (isinstance(self._params["mesh"], basestring) and len(self._params["mesh"]) == 3):
                if (self._params["mesh"][0] % 2 != 0 and self._params["mesh"][0] != -1) or (self._params["mesh"][1] % 2 != 0 and self._params["mesh"][1] != -1) or (self._params["mesh"][2] % 2 != 0 and self._params["mesh"][2] != -1):
                    raise ValueError(
                        "P3M requires an even number of mesh points in all directions")

            if not (self._params["cao"] >= -1 and self._params["cao"] <= 7):
                raise ValueError(
                    "P3M cao has to be an integer between -1 and 7")

            if not (self._params["accuracy"] > 0):
                raise ValueError("P3M accuracy has to be positive")

            if self._params["epsilon"] == "metallic":
                self._params = 0.0

            if not (isinstance(self._params["epsilon"], float) or self._params["epsilon"] == "metallic"):
                raise ValueError("epsilon should be a double or 'metallic'")

            if not (self._params["inter"] == default_params["inter"] or self._params["inter"] > 0):
                raise ValueError("inter should be a positive integer")

            if not (self._params["mesh_off"] == default_params["mesh_off"] or len(self._params["mesh_off"]) == 3):
                raise ValueError(
                    "mesh_off should be a list of length 3 and values between 0.0 and 1.0")

        def valid_keys(self):
            return "prefactor", "alpha_L", "r_cut_iL", "mesh", "mesh_off", "cao", "inter", "accuracy", "epsilon", "cao_cut", "a", "ai", "alpha", "r_cut", "inter2", "cao3", "additional_mesh", "bjerrum_length", "tune"

        def required_keys(self):
            return ["accuracy", ]

        def default_params(self):
            return {"cao": -1,
                    "inter": -1,
                    "r_cut": -1,
                    "accuracy": -1,
                    "mesh": -1,
                    "epsilon": 0.0,
                    "mesh_off": [-1, -1, -1],
                    "tune": True}

        def _get_params_from_es_core(self):
            params = {}
            params.update(dp3m.params)
            params["prefactor"] = coulomb.Dprefactor
            params["tune"] = self._params["tune"]
            return params

        def _set_params_in_es_core(self):
            self.set_magnetostatics_prefactor()
            dp3m_set_eps(self._params["epsilon"])
            dp3m_set_ninterpol(self._params["inter"])
            self.python_dp3m_set_mesh_offset(self._params["mesh_off"])
            self.python_dp3m_set_params(self._params["r_cut"], self._params["mesh"], self._params[
                "cao"], self._params["alpha"], self._params["accuracy"])

        def _tune(self):
            dipolar_set_Dbjerrum(self._params["bjerrum_length"])
            dp3m_set_eps(self._params["epsilon"])
            self.python_dp3m_set_tune_params(self._params["r_cut"], self._params["mesh"], self._params[
                "cao"], -1.0, self._params["accuracy"], self._params["inter"])
            resp, log = self.python_dp3m_adaptive_tune()
            if resp:
                raise Exception(
                    "failed to tune dipolar P3M parameters to required accuracy")
            print log
            self._params.update(self._get_params_from_es_core())

        def _activate_method(self):
            if self._params["tune"]:
                self._tune()

            self._set_params_in_es_core()

        def python_dp3m_set_mesh_offset(self, mesh_off):
            cdef double mesh_offset[3]
            mesh_offset[0] = mesh_off[0]
            mesh_offset[1] = mesh_off[1]
            mesh_offset[2] = mesh_off[2]
            return dp3m_set_mesh_offset(mesh_offset[0], mesh_offset[1], mesh_offset[2])

        def python_dp3m_adaptive_tune(self):
            cdef char * log = NULL
            cdef int response
            response = dp3m_adaptive_tune( & log)
            return response, log

        def python_dp3m_set_params(self, p_r_cut, p_mesh, p_cao, p_alpha, p_accuracy):
            cdef int mesh
            cdef double r_cut
            cdef int cao
            cdef double alpha
            cdef double accuracy
            r_cut = p_r_cut
            cao = p_cao
            alpha = p_alpha
            accuracy = p_accuracy
            if hasattr(p_mesh, "__getitem__"):
                mesh = p_mesh[0]
            else:
                mesh = p_mesh
            dp3m_set_params(r_cut, mesh, cao, alpha, accuracy)

        def python_dp3m_set_tune_params(self, p_r_cut, p_mesh, p_cao, p_alpha, p_accuracy, p_n_interpol):
            # cdef  python_p3m_set_tune_params():
            cdef int mesh
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
            mesh = p_mesh
            dp3m_set_tune_params(r_cut, mesh, cao, alpha, accuracy, n_interpol)

IF DIPOLES == 1:
    cdef class DipolarDirectSumCpu(MagnetostaticInteraction):

        """Calculates magnetostatic interactions by direct summation over all
           pairs. If the system has periodic boundaries, the minimum image
           convention is applied."""

        def default_params(self):
            return {}

        def required_keys(self):
            return ()

        def valid_keys(self):
            return ("bjerrum_length", "prefactor")

        def _get_params_from_es_core(self):
            return {"prefactor": coulomb.Dprefactor}

        def _activate_method(self):
            self._set_params_in_es_core(self)

        def _set_params_in_es_core(self):
            self.set_magnetostatics_prefactor()
            if dawaanr_set_params():
                raise Exception(
                    "Could not activate magnetostatics method " + self.__class__.__name__)

    cdef class DipolarDirectSumWithReplicaCpu(MagnetostaticInteraction):

        """Calculates magnetostatic interactions by direct summation over all
           pairs. If the system has periodic boundaries, n_replica
           copies are attached on each side. Spherical cutoff is applied."""

        def default_params(self):
            return {}

        def required_keys(self):
            return ("n_replica",)

        def valid_keys(self):
            return ("bjerrum_length", "prefactor", "n_replica")

        def _get_params_from_es_core(self):
            return {"prefactor": coulomb.Dprefactor, "n_replica": Ncut_off_magnetic_dipolar_direct_sum}

        def _activate_method(self):
            self._set_params_in_es_core(self)

        def _set_params_in_es_core(self):
            self.set_magnetostatics_prefactor()
            if mdds_set_params(self._params["n_replica"]):
                raise Exception(
                    "Could not activate magnetostatics method " + self.__class__.__name__)
