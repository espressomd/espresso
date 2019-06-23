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

from __future__ import print_function, absolute_import
from . cimport utils
include "myconfig.pxi"
from .actors import Actor

IF DIPOLES and DP3M:
    class MagnetostaticExtension(Actor):

        pass

    class DLC(MagnetostaticExtension):

        """Provide the Dipolar Layer Correction (DLC) method.

        DLC works like ELC for electrostatics
        (:class:`espressomd.electrostatic_extensions.ELC`),
        but applied to magnetic dipoles.

        Notes
        -----
        At present, the empty gap (volume without any particles), is assumed to be
        along the z-axis. As a reference for the DLC method, see :cite:`brodka04a`.

        Attributes
        ----------
        far_cut : :obj:`float`
                  Cutoff of the exponential sum.
        gap_size : :obj:`float`
                   Size of the empty gap. Note that DLC relies on the user to make sure that
                   this condition is fulfilled.
        maxPWerror : :obj:`float`
                     Maximal pairwise error of the potential and force.

        """

        def validate_params(self):
            """Check validity of class attributes.

            """
            default_params = self.default_params()
            check_type_or_throw_except(
                self._params["maxPWerror"], 1, float, "")
            check_range_or_except(
                self._params, "maxPWerror", 0, False, "inf", True)
            check_type_or_throw_except(self._params["gap_size"], 1, float, "")
            check_range_or_except(
                self._params, "gap_size", 0, False, "inf", True)
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
            handle_errors("mdlc tuning failed, gap size too small")

        def _activate_method(self):
            self._set_params_in_es_core()

        def _deactivate_method(self):
            pass
