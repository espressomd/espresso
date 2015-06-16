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
from actors import Actor

IF DIPOLES == 1:
    class MagnetostaticExtension(Actor):
        pass

    class DLC(MagnetostaticExtension):

        def validateParams(self):
            default_params = self.defaultParams()
            checkTypeOrExcept(self._params["maxPWerror"], 1, float, "")
            checkRangeOrExcept(
                self._params["maxPWerror"], 0, False, "inf", True)
            checkTypeOrExcept(self._params["gap_size"], 1, float, "")
            checkRangeOrExcept(self._params["gap_size"], 0, False, "inf", True)
            checkTypeOrExcept(self._params["far_cut"], 1, float, "")

        def validKeys(self):
            return "maxPWerror", "gap_size", "far_cut"

        def requiredKeys(self):
            return ["maxPWerror", "gap_size"]

        def defaultParams(self):
            return {"maxPWerror": -1,
                    "gap_size": -1,
                    "far_cut": -1}

        def _getParamsFromEsCore(self):
            params = {}
            params.update(dlc_params)
            return params

        def _setParamsInEsCore(self):
            if mdlc_set_params(self._params["maxPWerror"], self._params["gap_size"], self._params["far_cut"]):
                raise ValueError(
                    "Choose a 3d magnetostatics method prior to DLC")

        def _activateMethod(self):
            self._setParamsInEsCore()

        def _deactivateMethod(self):
            pass
