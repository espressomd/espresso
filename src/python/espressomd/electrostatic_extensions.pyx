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

from __future__ import print_function, absolute_import
from . cimport utils
include "myconfig.pxi"
from espressomd cimport actors
from . import actors
import numpy as np

IF ELECTROSTATICS and P3M:
    cdef class ElectrostaticExtensions(actors.Actor):
        pass

    cdef class ELC(ElectrostaticExtensions):

        def validate_params(self):
            default_params = self.default_params()
            check_type_or_throw_except(
                self._params["maxPWerror"], 1, float, "")
            check_range_or_except(
                self._params, "maxPWerror", 0, False, "inf", True)
            check_type_or_throw_except(self._params["gap_size"], 1, float, "")
            check_range_or_except(
                self._params, "gap_size", 0, False, "inf", True)
            check_type_or_throw_except(self._params["far_cut"], 1, float, "")
            check_type_or_throw_except(
                self._params["neutralize"], 1, type(True), "")

        def valid_keys(self):
            return "maxPWerror", "gap_size", "far_cut", "neutralize"

        def required_keys(self):
            return ["maxPWerror", "gap_size"]

        def default_params(self):
            return {"maxPWerror": -1,
                    "gap_size": -1,
                    "far_cut": -1,
                    "neutralize": True}

        def _get_params_from_es_core(self):
            params = {}
            params.update(elc_params)
            return params

        def _set_params_in_es_core(self):
            if coulomb.method == COULOMB_P3M_GPU:
                raise Exception(
                    "ELC tuning failed, ELC is not set up to work with the GPU P3M")
            if ELC_set_params(self._params["maxPWerror"], self._params["gap_size"], self._params["far_cut"], int(self._params["neutralize"]), 0, 0, 0, 0):
                raise ValueError(
                    "Choose a 3d electrostatics method prior to ELC")

        def _activate_method(self):
            self._set_params_in_es_core()

        def _deactivateMethod(self):
            pass

    cdef class ICC(ElectrostaticExtensions):

        def validate_params(self):
            default_params = self.default_params()

            check_type_or_throw_except(self._params["n_icc"], 1, int, "")
            check_range_or_except(
                self._params, "n_icc", 1, True, "inf", True)

            check_type_or_throw_except(
                self._params["convergence"], 1, float, "")
            check_range_or_except(
                self._params, "convergence", 0, False, "inf", True)

            check_type_or_throw_except(
                self._params["relaxation"], 1, float, "")
            check_range_or_except(
                self._params, "relaxation", 0, False, "inf", True)

            check_type_or_throw_except(
                self._params["ext_field"], 3, float, "")

            check_type_or_throw_except(
                self._params["max_iterations"], 1, int, "")
            check_range_or_except(
                self._params, "max_iterations", 0, False, "inf", True)

            check_type_or_throw_except(
                self._params["first_id"], 1, int, "")
            check_range_or_except(
                self._params, "first_id", 0, True, "inf", True)

            check_type_or_throw_except(
                self._params["eps_out"], 1, float, "")

            # Required list input
            self._params["normals"] = np.array(self._params["normals"])
            if self._params["normals"].size != self._params["n_icc"] * 3:
                raise ValueError(
                    "Expecting normal list with " + self._params["n_icc"] * 3 + " entries.")
            check_type_or_throw_except(self._params["normals"], self._params[
                "n_icc"], np.ndarray, "Error in normal list.")

            check_type_or_throw_except(
                self._params["areas"], self._params["n_icc"], float, "Error in area list.")

            # Not Required
            if self._params.has_key("sigmas"):
                check_type_or_throw_except(
                    self._params["sigmas"], self._params["n_icc"], float, "Error in sigma list.")
            else:
                self._params["sigmas"] = np.zeros(self._params["n_icc"])

            if self._params.has_key("epsilons"):
                check_type_or_throw_except(
                    self._params["epsilons"], self._params["n_icc"], float, "Error in epsilon list.")
            else:
                self._params["epsilons"] = np.zeros(self._params["n_icc"])

        def valid_keys(self):
            return "n_icc", "convergence", "relaxation", "ext_field", "max_iterations", "first_id", "eps_out", "normals", "areas", "sigmas", "epsilons"

        def required_keys(self):
            return ["n_icc", "normals", "areas"]

        def default_params(self):
            return {"n_icc": 0,
                    "convergence": 1e-3,
                    "relaxation": 0.7,
                    "ext_field": [0, 0, 0],
                    "max_iterations": 100,
                    "first_id": 0,
                    "esp_out": 1,
                    "normals": [],
                    "areas": [],
                    "sigmas": [],
                    "epsilons": []}

        def _get_params_from_es_core(self):
            params = {}
            params["n_icc"] = iccp3m_cfg.n_ic

            # Fill Lists
            normals = []
            areas = []
            sigmas = []
            epsilons = []
            for i in range(iccp3m_cfg.n_ic):
                normals.append([iccp3m_cfg.nvectorx[i], iccp3m_cfg.nvectory[
                               i], iccp3m_cfg.nvectorz[i]])
                areas.append(iccp3m_cfg.areas[i])
                epsilons.append(iccp3m_cfg.ein[i])
                sigmas.append(iccp3m_cfg.sigma[i])

            params["normals"] = normals
            params["areas"] = areas
            params["epsilons"] = epsilons
            params["sigmas"] = sigmas

            params["ext_field"] = [iccp3m_cfg.extx,
                                   iccp3m_cfg.exty, iccp3m_cfg.extz]
            params["first_id"] = iccp3m_cfg.first_id
            params["max_iterations"] = iccp3m_cfg.num_iteration
            params["convergence"] = iccp3m_cfg.convergence
            params["relaxation"] = iccp3m_cfg.relax
            params["eps_out"] = iccp3m_cfg.eout

            return params

        def _set_params_in_es_core(self):

            # First set number of icc particles
            iccp3m_cfg.n_ic = self._params["n_icc"]
            # Allocate ICC lists
            iccp3m_alloc_lists()

            # Fill Lists
            for i in range(iccp3m_cfg.n_ic):
                iccp3m_cfg.nvectorx[i] = self._params["normals"][i][0]
                iccp3m_cfg.nvectory[i] = self._params["normals"][i][1]
                iccp3m_cfg.nvectorz[i] = self._params["normals"][i][2]
                iccp3m_cfg.areas[i] = self._params["areas"][i]
                iccp3m_cfg.ein[i] = self._params["epsilons"][i]
                iccp3m_cfg.sigma[i] = self._params["sigmas"][i]

            iccp3m_cfg.extx = self._params["ext_field"][0]
            iccp3m_cfg.exty = self._params["ext_field"][1]
            iccp3m_cfg.extz = self._params["ext_field"][2]
            iccp3m_cfg.first_id = self._params["first_id"]
            iccp3m_cfg.num_iteration = self._params["max_iterations"]
            iccp3m_cfg.convergence = self._params["convergence"]
            iccp3m_cfg.relax = self._params["relaxation"]
            iccp3m_cfg.eout = self._params["eps_out"]
            iccp3m_cfg.citeration = 0

            iccp3m_set_initialized()
            iccp3m_cfg.set_flag = 1

            # Broadcasts vars
            mpi_iccp3m_init(0)

        def _activate_method(self):
            self._set_params_in_es_core()

        def _deactivate_method(self):
            pass
