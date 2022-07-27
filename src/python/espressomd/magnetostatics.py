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
from .__init__ import has_features


class MagnetostaticInteraction(ScriptInterfaceHelper):
    """
    Common interface for magnetostatics solvers.

    Parameters
    ----------
    prefactor : :obj:`float`
        Magnetostatics prefactor :math:`\\frac{\\mu_0\\mu}{4\\pi}`

    """
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
        if not has_features("DIPOLES"):
            raise NotImplementedError("Feature DIPOLES not compiled in")

    def validate_params(self, params):
        """Check validity of given parameters.
        """
        utils.check_type_or_throw_except(
            params["prefactor"], 1, float, "prefactor should be a double")

    def default_params(self):
        raise NotImplementedError("Derived classes must implement this method")

    def valid_keys(self):
        raise NotImplementedError("Derived classes must implement this method")

    def required_keys(self):
        raise NotImplementedError("Derived classes must implement this method")

    def _activate(self):
        self.call_method("activate")
        utils.handle_errors("Dipolar actor activation failed")

    def _deactivate(self):
        self.call_method("deactivate")
        utils.handle_errors("Dipolar actor deactivation failed")

    def get_magnetostatics_prefactor(self):
        """
        Get the magnetostatics prefactor

        """
        return self.prefactor


@script_interface_register
class DipolarP3M(MagnetostaticInteraction):
    """
    Calculate magnetostatic interactions using the dipolar P3M method.
    See :ref:`Dipolar P3M` for more details.

    Parameters
    ----------
    prefactor : :obj:`float`
        Magnetostatics prefactor (:math:`\\mu_0/(4\\pi)`)
    accuracy : :obj:`float`
        P3M tunes its parameters to provide this target accuracy.
    alpha : :obj:`float`
        Ewald parameter.
    cao : :obj:`int`
        Charge-assignment order, an integer between 1 and 7.
    mesh : :obj:`int` or (3,) array_like of :obj:`int`
        The number of mesh points in x, y and z direction. Use a single
        value for cubic boxes.
    mesh_off : (3,) array_like of :obj:`float`, optional
        Mesh offset.
    r_cut : :obj:`float`
        Real space cutoff.
    tune : :obj:`bool`, optional
        Activate/deactivate the tuning method on activation
        (default is ``True``, i.e., activated).
    timings : :obj:`int`
        Number of force calculations during tuning.

    """
    _so_name = "Dipoles::DipolarP3M"

    def _check_required_features(self):
        if not has_features("DP3M"):
            raise NotImplementedError("Feature DP3M not compiled in")

    def validate_params(self, params):
        """Check validity of parameters.

        """
        super().validate_params(params)

        if utils.is_valid_type(params["mesh"], int):
            params["mesh"] = 3 * [params["mesh"]]
        else:
            utils.check_type_or_throw_except(params["mesh"], 3, int,
                                             "DipolarP3M mesh has to be an integer or integer list of length 3")

        if params["epsilon"] == "metallic":
            params["epsilon"] = 0.0

        utils.check_type_or_throw_except(
            params["epsilon"], 1, float,
            "epsilon should be a double or 'metallic'")

        utils.check_type_or_throw_except(params["mesh_off"], 3, float,
                                         "mesh_off should be a (3,) array_like of values between 0.0 and 1.0")

        if not utils.is_valid_type(params["timings"], int):
            raise TypeError("DipolarP3M timings has to be an integer")
        if params["timings"] <= 0:
            raise ValueError("DipolarP3M timings must be > 0")
        if not utils.is_valid_type(params["tune"], bool):
            raise TypeError("DipolarP3M tune has to be a boolean")

    def valid_keys(self):
        return {"prefactor", "alpha_L", "r_cut_iL", "mesh", "mesh_off",
                "cao", "accuracy", "epsilon", "cao_cut", "a", "ai",
                "alpha", "r_cut", "cao3", "tune", "timings", "verbose"}

    def required_keys(self):
        return {"accuracy"}

    def default_params(self):
        return {"cao": -1,
                "r_cut": -1,
                "alpha": -1,
                "accuracy": -1,
                "mesh": [-1, -1, -1],
                "epsilon": 0.0,
                "mesh_off": [0.5, 0.5, 0.5],
                "prefactor": 0.,
                "tune": True,
                "timings": 10,
                "verbose": True}


@script_interface_register
class DipolarDirectSumCpu(MagnetostaticInteraction):
    """
    Calculate magnetostatic interactions by direct summation over all pairs.
    See :ref:`Dipolar direct sum` for more details.

    If the system has periodic boundaries, the minimum image convention is
    applied in the respective directions.

    Parameters
    ----------
    prefactor : :obj:`float`
        Magnetostatics prefactor (:math:`\\mu_0/(4\\pi)`)

    """
    _so_name = "Dipoles::DipolarDirectSum"

    def default_params(self):
        return {}

    def required_keys(self):
        return set()

    def valid_keys(self):
        return {"prefactor"}


@script_interface_register
class DipolarDirectSumWithReplicaCpu(MagnetostaticInteraction):
    """
    Calculate magnetostatic interactions by direct summation over all pairs.
    See :ref:`Dipolar direct sum` for more details.

    If the system has periodic boundaries, ``n_replica`` copies of the system are
    taken into account in the respective directions. Spherical cutoff is applied.

    Parameters
    ----------
    prefactor : :obj:`float`
        Magnetostatics prefactor (:math:`\\mu_0/(4\\pi)`)
    n_replica : :obj:`int`
        Number of replicas to be taken into account at periodic boundaries.

    """
    _so_name = "Dipoles::DipolarDirectSumWithReplica"

    def default_params(self):
        return {}

    def required_keys(self):
        return {"n_replica"}

    def valid_keys(self):
        return {"prefactor", "n_replica"}


@script_interface_register
class Scafacos(MagnetostaticInteraction):

    """
    Calculate the dipolar interaction using dipoles-capable methods
    from the ScaFaCoS library. See :ref:`ScaFaCoS magnetostatics` for
    more details.

    Parameters
    ----------
    prefactor : :obj:`float`
        Magnetostatics prefactor (:math:`\\mu_0/(4\\pi)`).
    method_name : :obj:`str`
        Name of the ScaFaCoS method to use.
    method_params : :obj:`dict`
        Dictionary with the key-value pairs of the method parameters as
        defined in ScaFaCoS. Note that the values are cast to strings
        to match ScaFaCoS' interface.

    Methods
    -------
    get_available_methods()
        List long-range methods available in the ScaFaCoS library.

    """
    _so_name = "Dipoles::DipolarScafacos"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = MagnetostaticInteraction._so_bind_methods + \
        ("get_available_methods", )

    def _check_required_features(self):
        if not has_features("DIPOLES"):
            raise NotImplementedError("Feature DIPOLES not compiled in")
        if not has_features("SCAFACOS_DIPOLES"):
            raise NotImplementedError(
                "Feature SCAFACOS_DIPOLES not compiled in")

    def validate_params(self, params):
        pass

    def default_params(self):
        return {}

    def valid_keys(self):
        return {"method_name", "method_params", "prefactor"}

    def required_keys(self):
        return {"method_name", "method_params", "prefactor"}


@script_interface_register
class DipolarDirectSumGpu(MagnetostaticInteraction):
    """
    Calculate magnetostatic interactions by direct summation over all
    pairs. See :ref:`Dipolar direct sum` for more details.

    If the system has periodic boundaries, the minimum image convention
    is applied in the respective directions.

    This is the GPU version of :class:`espressomd.magnetostatics.DipolarDirectSumCpu`
    but uses floating point precision.

    Requires feature ``DIPOLAR_DIRECT_SUM``, which depends on
    ``DIPOLES`` and ``CUDA``.

    Parameters
    ----------
    prefactor : :obj:`float`
        Magnetostatics prefactor (:math:`\\mu_0/(4\\pi)`)

    """
    _so_name = "Dipoles::DipolarDirectSumGpu"
    _so_creation_policy = "GLOBAL"

    def _check_required_features(self):
        if not has_features("DIPOLAR_DIRECT_SUM"):
            raise NotImplementedError(
                "Features CUDA and DIPOLES not compiled in")

    def default_params(self):
        return {}

    def required_keys(self):
        return set()

    def valid_keys(self):
        return {"prefactor"}


@script_interface_register
class DipolarBarnesHutGpu(MagnetostaticInteraction):

    """
    Calculates magnetostatic interactions by direct summation over all
    pairs. See :ref:`Barnes-Hut octree sum on GPU` for more details.

    TODO: If the system has periodic boundaries, the minimum image
    convention is applied.

    Requires feature ``DIPOLAR_BARNES_HUT``, which depends on
    ``DIPOLES`` and ``CUDA``.

    """
    _so_name = "Dipoles::DipolarBarnesHutGpu"
    _so_creation_policy = "GLOBAL"

    def _check_required_features(self):
        if not has_features("DIPOLAR_BARNES_HUT"):
            raise NotImplementedError(
                "Features CUDA and DIPOLES not compiled in")

    def default_params(self):
        return {"epssq": 100.0, "itolsq": 4.0}

    def required_keys(self):
        return set()

    def valid_keys(self):
        return {"prefactor", "epssq", "itolsq"}


@script_interface_register
class DLC(MagnetostaticInteraction):

    """
    Electrostatics solver for systems with two periodic dimensions.
    See :ref:`Dipolar Layer Correction (DLC)` for more details.

    Notes
    -----
    At present, the empty gap (volume without any particles), is assumed to be
    along the z-axis. As a reference for the DLC method, see :cite:`brodka04a`.

    Parameters
    ----------
    gap_size : :obj:`float`
        The gap size gives the height :math:`h` of the empty region between
        the system box and the neighboring artificial images. |es| checks
        that the gap is empty and will throw an error if it isn't. Therefore
        you should really make sure that the gap region is empty (e.g.
        with wall constraints).
    maxPWerror : :obj:`float`
        Maximal pairwise error of the potential and force.
    far_cut : :obj:`float`, optional
        Cutoff of the exponential sum.

    """
    _so_name = "Dipoles::DipolarLayerCorrection"
    _so_creation_policy = "GLOBAL"

    def validate_params(self, params):
        utils.check_type_or_throw_except(
            params["maxPWerror"], 1, float, "maxPWerror has to be a float")
        utils.check_type_or_throw_except(
            params["gap_size"], 1, float, "gap_size has to be a float")
        utils.check_type_or_throw_except(
            params["far_cut"], 1, float, "far_cut has to be a float")

    def default_params(self):
        return {"far_cut": -1.}

    def valid_keys(self):
        return {"actor", "maxPWerror", "gap_size", "far_cut"}

    def required_keys(self):
        return {"actor", "maxPWerror", "gap_size"}
