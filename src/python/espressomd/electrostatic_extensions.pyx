#
# Copyright (C) 2013,2014 The ESPResSo project
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

cimport utils
include "myconfig.pxi"
cimport actors
import actors
import numpy as np

IF ELECTROSTATICS and P3M:
	cdef class ElectrostaticExtensions(actors.Actor):
		pass

	cdef class ELC(ElectrostaticExtensions):

		def validateParams(self):
			default_params = self.defaultParams()
			checkTypeOrExcept(self._params["maxPWerror"], 1, float, "")
			checkRangeOrExcept(self._params["maxPWerror"], 0, False, "inf", True)
			checkTypeOrExcept(self._params["gap_size"], 1, float, "")
			checkRangeOrExcept(self._params["gap_size"], 0, False, "inf", True)
			checkTypeOrExcept(self._params["far_cut"], 1, float, "")
			checkTypeOrExcept(self._params["neutralize"], 1, type(True), "")

		def validKeys(self):
			return "maxPWerror", "gap_size", "far_cut", "neutralize"

		def requiredKeys(self):
			return ["maxPWerror", "gap_size"]

		def defaultParams(self):
			return {"maxPWerror": -1,
				"gap_size": -1,
				"far_cut": -1,
				"neutralize": True}

		def _getParamsFromEsCore(self):
			params = {}
			params.update(elc_params)
			return params

		def _setParamsInEsCore(self):
			if coulomb.method == COULOMB_P3M_GPU:
				raise Exception("ELC tuning failed, ELC is not set up to work with the GPU P3M")
			if ELC_set_params(self._params["maxPWerror"], self._params["gap_size"], self._params["far_cut"], int(self._params["neutralize"]), 0, 0, 0, 0):
				raise ValueError("Choose a 3d electrostatics method prior to ELC")

		def _activateMethod(self):
			self._setParamsInEsCore()

		def _deactivateMethod(self):
			pass

	cdef class ICC(ElectrostaticExtensions):

		def validateParams(self):
			default_params = self.defaultParams()

			checkTypeOrExcept(self._params["n_icc"], 1, int, "")
			checkRangeOrExcept(self._params["n_icc"], 1, True, "inf", True)
            
			checkTypeOrExcept(self._params["convergence"], 1, float, "")
			checkRangeOrExcept(self._params["convergence"], 0, False, "inf", True)
            
			checkTypeOrExcept(self._params["relaxation"], 1, float, "")
			checkRangeOrExcept(self._params["relaxation"], 0, False, "inf", True)
            
			checkTypeOrExcept(self._params["ext_field"], 3, float, "")

			checkTypeOrExcept(self._params["max_iterations"], 1, int, "")
			checkRangeOrExcept(self._params["max_iterations"], 0, False, "inf", True)

			checkTypeOrExcept(self._params["first_id"], 1, int, "")
			checkRangeOrExcept(self._params["first_id"], 0, True, "inf", True)

			checkTypeOrExcept(self._params["eps_out"], 1, float, "")

			#Required list input
			self._params["normals"]=np.array(self._params["normals"])
			if self._params["normals"].size != self._params["n_icc"] * 3:
				raise ValueError("Expecting normal list with " +  self._params["n_icc"] * 3 + " entries.")
			checkTypeOrExcept(self._params["normals"], self._params["n_icc"], np.ndarray, "Error in normal list.")

			checkTypeOrExcept(self._params["areas"], self._params["n_icc"], float, "Error in area list.")

			#Not Required
			if self._params.has_key("sigmas"):
				checkTypeOrExcept(self._params["sigmas"], self._params["n_icc"], float, "Error in sigma list.")
			else:
				self._params["sigmas"]=np.zeros(self._params["n_icc"])

			if self._params.has_key("epsilons"):
				checkTypeOrExcept(self._params["epsilons"], self._params["n_icc"], float, "Error in epsilon list.")
			else:
				self._params["epsilons"]=np.zeros(self._params["n_icc"])

		def validKeys(self):
			return "n_icc", "convergence", "relaxation", "ext_field", "max_iterations", "first_id", "eps_out", "normals", "areas", "sigmas", "epsilons"
		def requiredKeys(self):
			return ["n_icc", "normals", "areas"]

		def defaultParams(self):
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

		def _getParamsFromEsCore(self):
			params = {}
			params.update(iccp3m_cfg)
			return params

		def _setParamsInEsCore(self):
			iccp3m_cfg.n_ic = self._params["n_icc"]

			#Fill Lists
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
			iccp3m_cfg.num_iteration=self._params["first_id"]
			iccp3m_cfg.convergence=self._params["convergence"]
			iccp3m_cfg.relax=self._params["relaxation"]
			iccp3m_cfg.eout=self._params["eps_out"]
			iccp3m_cfg.citeration=0
			iccp3m_initialized=1
			iccp3m_cfg.set_flag = 1

			mpi_iccp3m_init(0)


		def _activateMethod(self):
			self._setParamsInEsCore()

		def _deactivateMethod(self):
			pass
