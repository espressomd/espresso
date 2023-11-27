#
# Copyright (C) 2013-2022 The ESPResSo project
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

from . import utils
from .script_interface import ScriptInterfaceHelper, script_interface_register
from .code_features import has_features
import numpy as np


class ElectrostaticExtensions(ScriptInterfaceHelper):

    _so_creation_policy = "GLOBAL"

    def __init__(self, **kwargs):
        self._check_required_features()

        if 'sip' not in kwargs:
            params = self.default_params()
            params.update(kwargs)
            self.validate_params(params)
            super().__init__(**params)
        else:
            super().__init__(**kwargs)

    def _check_required_features(self):
        if not has_features("ELECTROSTATICS"):
            raise NotImplementedError("Feature ELECTROSTATICS not compiled in")

    def validate_params(self, params):
        raise NotImplementedError("Derived classes must implement this method")

    def default_params(self):
        raise NotImplementedError("Derived classes must implement this method")

    def valid_keys(self):
        raise NotImplementedError("Derived classes must implement this method")

    def required_keys(self):
        raise NotImplementedError("Derived classes must implement this method")


@script_interface_register
class ICC(ElectrostaticExtensions):
    """
    Interface to the induced charge calculation scheme for dielectric
    interfaces. See :ref:`Dielectric interfaces with the ICC algorithm`
    for more details.

    Parameters
    ----------
    n_icc : :obj:`int`
        Total number of ICC Particles.
    first_id : :obj:`int`, optional
        ID of the first ICC Particle.
    convergence : :obj:`float`, optional
        Abort criteria of the iteration. It corresponds to the maximum relative
        change of any of the interface particle's charge.
    relaxation : :obj:`float`, optional
        SOR relaxation parameter.
    ext_field : :obj:`float`, optional
        Homogeneous electric field added to the calculation of dielectric boundary forces.
    max_iterations : :obj:`int`, optional
        Maximal number of iterations.
    eps_out : :obj:`float`, optional
        Relative permittivity of the outer region (where the particles are).
    normals : (``n_icc``, 3) array_like :obj:`float`
        Normal vectors pointing into the outer region.
    areas : (``n_icc``, ) array_like :obj:`float`
        Areas of the discretized surface.
    sigmas : (``n_icc``, ) array_like :obj:`float`, optional
        Additional surface charge density in the absence of any charge
        induction.
    epsilons : (``n_icc``, ) array_like :obj:`float`
        Dielectric constant associated to the areas.

    """
    _so_name = "Coulomb::ICCStar"
    _so_creation_policy = "GLOBAL"

    def validate_params(self, params):
        utils.check_type_or_throw_except(
            params["n_icc"], 1, int, "Invalid parameter 'n_icc'")
        utils.check_type_or_throw_except(
            params["first_id"], 1, int, "Invalid parameter 'first_id'")
        utils.check_type_or_throw_except(
            params["convergence"], 1, float, "Invalid parameter 'convergence'")
        utils.check_type_or_throw_except(
            params["relaxation"], 1, float, "Invalid parameter 'relaxation'")
        utils.check_type_or_throw_except(
            params["ext_field"], 3, float, "Invalid parameter 'ext_field'")
        utils.check_type_or_throw_except(
            params["max_iterations"], 1, int, "Invalid parameter 'max_iterations'")
        utils.check_type_or_throw_except(
            params["eps_out"], 1, float, "Invalid parameter 'eps_out'")

        n_icc = params["n_icc"]
        if n_icc <= 0:
            raise ValueError("Parameter 'n_icc' must be >= 1")

        if n_icc:
            if np.shape(params["normals"]) != (n_icc, 3):
                raise ValueError("Parameter 'normals' has incorrect shape")
        utils.check_array_type_or_throw_except(
            np.reshape(params["normals"], (-1,)), 3 * n_icc, float,
            "Parameter 'normals' has incorrect type")

        if "sigmas" not in params:
            params["sigmas"] = np.zeros(n_icc)

        for key in ("areas", "sigmas", "epsilons"):
            if np.shape(params[key]) != (n_icc,):
                raise ValueError(f"Parameter '{key}' has incorrect shape")
            utils.check_array_type_or_throw_except(
                np.reshape(params[key], (-1,)), n_icc, float,
                f"Parameter '{key}' has incorrect type")

    def valid_keys(self):
        return {"n_icc", "convergence", "relaxation", "ext_field",
                "max_iterations", "first_id", "eps_out", "normals",
                "areas", "sigmas", "epsilons", "check_neutrality"}

    def required_keys(self):
        return {"n_icc", "normals", "areas", "epsilons"}

    def default_params(self):
        return {"convergence": 1e-3,
                "relaxation": 0.7,
                "ext_field": [0., 0., 0.],
                "max_iterations": 100,
                "first_id": 0,
                "eps_out": 1,
                "check_neutrality": True}

    def last_iterations(self):
        """
        Number of iterations needed in last relaxation to
        reach the convergence criterion.

        Returns
        -------
        iterations : :obj:`int`
            Number of iterations

        """
        return self.citeration
