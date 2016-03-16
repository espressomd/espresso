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

cimport utils
include "myconfig.pxi"
from actors import Actor

IF DIPOLES == 1:
    class MagnetostaticExtension(Actor):

        pass

    class DLC(MagnetostaticExtension):

        def validate_params(self):
            default_params = self.default_params()
            check_type_or_throw_except(
                self._params["maxPWerror"], 1, float, "")
            check_range_or_except(
                self._params["maxPWerror"], 0, False, "inf", True)
            check_type_or_throw_except(self._params["gap_size"], 1, float, "")
            check_range_or_except(
                self._params["gap_size"], 0, False, "inf", True)
            check_type_or_throw_except(self._params["far_cut"], 1, float, "")

        def valid_keys(self):
            return "maxPWerror", "gap_size", "far_cut"

        def required_keys(self):
            return ["maxPWerror", "gap_size"]

        def default_params(self):
            return {"maxPWerror": -1,
                    "gap_size": -1,
                    "far_cut": -1}

        def _get_params_from_es_core(self):
            params = {}
            params.update(dlc_params)
            return params

        def _set_params_in_es_core(self):
            if mdlc_set_params(self._params["maxPWerror"], self._params["gap_size"], self._params["far_cut"]):
                raise ValueError(
                    "Choose a 3d magnetostatics method prior to DLC")

        def _activate_method(self):
            self._set_params_in_es_core()

        def _deactivate_method(self):
            pass
