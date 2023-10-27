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


@script_interface_register
class Container(ScriptInterfaceHelper):
    _so_name = "Coulomb::Container"
    _so_features = ("ELECTROSTATICS",)
    _so_bind_methods = ("clear",)


class ElectrostaticInteraction(ScriptInterfaceHelper):
    """
    Common interface for electrostatics solvers.

    Parameters
    ----------
    prefactor : :obj:`float`
        Electrostatics prefactor :math:`\\frac{1}{4\\pi\\varepsilon_0\\varepsilon_r}`

    """
    _so_creation_policy = "GLOBAL"
    _so_features = ("ELECTROSTATICS",)

    def __init__(self, **kwargs):
        if 'sip' not in kwargs:
            for key in self.required_keys():
                if key not in kwargs:
                    raise RuntimeError(f"Parameter '{key}' is missing")
            params = self.default_params()
            params.update(kwargs)
            self.validate_params(params)
            super().__init__(**params)
            for key in params:
                if key not in self._valid_parameters():
                    raise RuntimeError(
                        f"Parameter '{key}' is not a valid parameter")
        else:
            super().__init__(**kwargs)

    def validate_params(self, params):
        """Check validity of given parameters.
        """
        utils.check_type_or_throw_except(
            params["prefactor"], 1, float, "Parameter 'prefactor' should be a float")

    def default_params(self):
        raise NotImplementedError("Derived classes must implement this method")

    def required_keys(self):
        raise NotImplementedError("Derived classes must implement this method")


@script_interface_register
class DH(ElectrostaticInteraction):
    """
    Electrostatics solver based on the Debye-Hueckel framework.
    See :ref:`Debye-Hückel potential` for more details.

    Parameters
    ----------
    prefactor : :obj:`float`
        Electrostatics prefactor (see :eq:`coulomb_prefactor`).
    kappa : :obj:`float`
        Inverse Debye screening length.
    r_cut : :obj:`float`
        Cutoff radius for this interaction.

    """
    _so_name = "Coulomb::DebyeHueckel"
    _so_creation_policy = "GLOBAL"

    def required_keys(self):
        return {"prefactor", "kappa", "r_cut"}

    def default_params(self):
        return {"check_neutrality": True}


@script_interface_register
class ReactionField(ElectrostaticInteraction):
    """
    Electrostatics solver based on the Reaction Field framework.
    See :ref:`Reaction Field method` for more details.

    Parameters
    ----------
    prefactor : :obj:`float`
        Electrostatics prefactor (see :eq:`coulomb_prefactor`).
    kappa : :obj:`float`
        Inverse Debye screening length.
    epsilon1 : :obj:`float`
        interior dielectric constant
    epsilon2 : :obj:`float`
        exterior dielectric constant
    r_cut : :obj:`float`
        Cutoff radius for this interaction.

    """
    _so_name = "Coulomb::ReactionField"
    _so_creation_policy = "GLOBAL"

    def required_keys(self):
        return {"prefactor", "kappa", "epsilon1", "epsilon2", "r_cut"}

    def default_params(self):
        return {"check_neutrality": True}


class _P3MBase(ElectrostaticInteraction):

    def required_keys(self):
        return {"prefactor", "accuracy"}

    def default_params(self):
        return {"cao": -1,
                "r_cut": -1.,
                "alpha": -1.,
                "mesh": [-1, -1, -1],
                "epsilon": 0.,
                "mesh_off": [-1., -1., -1.],
                "prefactor": 0.,
                "check_neutrality": True,
                "check_complex_residuals": True,
                "tune": True,
                "timings": 10,
                "verbose": True}

    def validate_params(self, params):
        super().validate_params(params)

        if utils.is_valid_type(params["mesh"], int):
            params["mesh"] = 3 * [params["mesh"]]
        utils.check_type_or_throw_except(
            params["mesh"], 3, int,
            "Parameter 'mesh' has to be an integer or integer list of length 3")
        if (params["mesh"][0] % 2 != 0 and params["mesh"][0] != -1) or \
           (params["mesh"][1] % 2 != 0 and params["mesh"][1] != -1) or \
           (params["mesh"][2] % 2 != 0 and params["mesh"][2] != -1):
            raise ValueError(
                "P3M requires an even number of mesh points in all directions")

        if params["epsilon"] == "metallic":
            params["epsilon"] = 0.0

        utils.check_type_or_throw_except(
            params["epsilon"], 1, float,
            "Parameter 'epsilon' has to be a float or 'metallic'")

        utils.check_type_or_throw_except(
            params["mesh_off"], 3, float,
            "Parameter 'mesh_off' has to be a (3,) array_like of values between 0.0 and 1.0")

        if not utils.is_valid_type(params["timings"], int):
            raise TypeError("Parameter 'timings' has to be an integer")
        if not utils.is_valid_type(params["tune"], bool):
            raise TypeError("Parameter 'tune' has to be a boolean")


@script_interface_register
class P3M(_P3MBase):
    """
    P3M electrostatics solver.

    Particle--Particle--Particle--Mesh (P3M) is a Fourier-based Ewald
    summation method to calculate potentials in N-body simulation.
    See :ref:`Coulomb P3M` for more details.

    Parameters
    ----------
    prefactor : :obj:`float`
        Electrostatics prefactor (see :eq:`coulomb_prefactor`).
    accuracy : :obj:`float`
        P3M tunes its parameters to provide this target accuracy.
    alpha : :obj:`float`, optional
        The Ewald parameter.
    cao : :obj:`float`, optional
        The charge-assignment order, an integer between 1 and 7.
    epsilon : :obj:`float` or :obj:`str`, optional
        A positive number for the dielectric constant of the
        surrounding medium. Use ``'metallic'`` to set the dielectric
        constant of the surrounding medium to infinity (default).
    mesh : :obj:`int` or (3,) array_like of :obj:`int`, optional
        The number of mesh points in x, y and z direction. Use a single
        value for cubic boxes.
    mesh_off : (3,) array_like of :obj:`float`, optional
        Mesh offset.
    r_cut : :obj:`float`, optional
        The real space cutoff.
    tune : :obj:`bool`, optional
        Used to activate/deactivate the tuning method on activation.
        Defaults to ``True``.
    timings : :obj:`int`
        Number of force calculations during tuning.
    verbose : :obj:`bool`, optional
        If ``False``, disable log output during tuning.
    check_neutrality : :obj:`bool`, optional
        Raise a warning if the system is not electrically neutral when
        set to ``True`` (default).
    check_complex_residuals: :obj:`bool`, optional
        Raise a warning if the backward Fourier transform has non-zero
        complex residuals when set to ``True`` (default).

    """
    _so_name = "Coulomb::CoulombP3M"
    _so_creation_policy = "GLOBAL"
    _so_features = ("P3M",)


@script_interface_register
class P3MGPU(_P3MBase):
    """
    P3M electrostatics solver with GPU support.

    Particle--Particle--Particle--Mesh (P3M) is a Fourier-based Ewald
    summation method to calculate potentials in N-body simulation.
    See :ref:`Coulomb P3M on GPU` for more details.

    Parameters
    ----------
    prefactor : :obj:`float`
        Electrostatics prefactor (see :eq:`coulomb_prefactor`).
    accuracy : :obj:`float`
        P3M tunes its parameters to provide this target accuracy.
    alpha : :obj:`float`, optional
        The Ewald parameter.
    cao : :obj:`float`, optional
        The charge-assignment order, an integer between 0 and 7.
    epsilon : :obj:`float` or :obj:`str`, optional
        A positive number for the dielectric constant of the
        surrounding medium. Use ``'metallic'`` to set the dielectric
        constant of the surrounding medium to infinity (default).
    mesh : :obj:`int` or (3,) array_like of :obj:`int`, optional
        The number of mesh points in x, y and z direction. Use a single
        value for cubic boxes.
    mesh_off : (3,) array_like of :obj:`float`, optional
        Mesh offset.
    r_cut : :obj:`float`, optional
        The real space cutoff
    tune : :obj:`bool`, optional
        Used to activate/deactivate the tuning method on activation.
        Defaults to ``True``.
    timings : :obj:`int`
        Number of force calculations during tuning.
    verbose : :obj:`bool`, optional
        If ``False``, disable log output during tuning.
    check_neutrality : :obj:`bool`, optional
        Raise a warning if the system is not electrically neutral when
        set to ``True`` (default).
    check_complex_residuals: :obj:`bool`, optional
        Raise a warning if the backward Fourier transform has non-zero
        complex residuals when set to ``True`` (default).

    """
    _so_name = "Coulomb::CoulombP3MGPU"
    _so_creation_policy = "GLOBAL"
    _so_features = ("P3M", "CUDA")


@script_interface_register
class ELC(ElectrostaticInteraction):
    """
    Electrostatics solver for systems with two periodic dimensions.
    See :ref:`Electrostatic Layer Correction (ELC)` for more details.

    Parameters
    ----------
    actor : :obj:`P3M`, required
        Base P3M actor.
    gap_size : :obj:`float`, required
        The gap size gives the height :math:`h` of the empty region between
        the system box and the neighboring artificial images. |es| checks
        that the gap is empty and will throw an error if it isn't. Therefore
        you should really make sure that the gap region is empty (e.g.
        with wall constraints).
    maxPWerror : :obj:`float`, required
        The maximal pairwise error sets the least upper bound (LUB) error
        of the force between any two charges without prefactors (see the
        papers). The algorithm tries to find parameters to meet this LUB
        requirements or will throw an error if there are none.
    delta_mid_top : :obj:`float`, optional
        Dielectric contrast :math:`\\Delta_t` between the upper boundary
        and the simulation box. Value between -1 and +1 (inclusive).
    delta_mid_bottom : :obj:`float`, optional
        Dielectric contrast :math:`\\Delta_b` between the lower boundary
        and the simulation box. Value between -1 and +1 (inclusive).
    const_pot : :obj:`bool`, optional
        Activate a constant electric potential between the top and bottom
        of the simulation box.
    pot_diff : :obj:`float`, optional
        If ``const_pot`` is enabled, this parameter controls the applied
        voltage between the boundaries of the simulation box in the
        *z*-direction (at :math:`z = 0` and :math:`z = L_z - h`).
    neutralize : :obj:`bool`, optional
        By default, *ELC* just as P3M adds a homogeneous neutralizing
        background to the system in case of a net charge. However, unlike
        in three dimensions, this background adds a parabolic potential
        across the slab :cite:`ballenegger09a`. Therefore, under normal
        circumstances, you will probably want to disable the neutralization
        for non-neutral systems. This corresponds then to a formal
        regularization of the forces and energies :cite:`ballenegger09a`.
        Also, if you add neutralizing walls explicitly as constraints, you
        have to disable the neutralization. When using a dielectric
        contrast or full metallic walls (``delta_mid_top != 0`` or
        ``delta_mid_bot != 0`` or ``const_pot=True``), ``neutralize`` is
        overwritten and switched off internally. Note that the special
        case of non-neutral systems with a *non-metallic* dielectric jump
        (e.g. ``delta_mid_top`` or ``delta_mid_bot`` in ``]-1,1[``) is not
        covered by the algorithm and will throw an error.
    far_cut : :obj:`float`, optional
        Cutoff radius, use with care, intended for testing purposes. When
        setting the cutoff directly, the maximal pairwise error is ignored.
    """
    _so_name = "Coulomb::ElectrostaticLayerCorrection"
    _so_creation_policy = "GLOBAL"
    _so_features = ("P3M",)

    def validate_params(self, params):
        utils.check_type_or_throw_except(
            params["maxPWerror"], 1, float, "maxPWerror has to be a float")
        utils.check_type_or_throw_except(
            params["gap_size"], 1, float, "gap_size has to be a float")
        utils.check_type_or_throw_except(
            params["far_cut"], 1, float, "far_cut has to be a float")
        utils.check_type_or_throw_except(
            params["neutralize"], 1, bool, "neutralize has to be a bool")

    def required_keys(self):
        return {"actor", "maxPWerror", "gap_size"}

    def default_params(self):
        return {"far_cut": -1.,
                "delta_mid_top": 0.,
                "delta_mid_bot": 0.,
                "const_pot": False,
                "pot_diff": 0.,
                "neutralize": True,
                "check_neutrality": True}


@script_interface_register
class MMM1D(ElectrostaticInteraction):
    """
    Electrostatics solver for systems with one periodic direction.
    See :ref:`MMM1D` for more details.

    Parameters
    ----------
    prefactor : :obj:`float`
        Electrostatics prefactor (see :eq:`coulomb_prefactor`).
    maxWPerror : :obj:`float`
        Maximal pairwise error.
    far_switch_radius : :obj:`float`, optional
        Radius where near-field and far-field calculation are switched.
    verbose : :obj:`bool`, optional
        If ``False``, disable log output during tuning.
    timings : :obj:`int`, optional
        Number of force calculations during tuning.
    check_neutrality : :obj:`bool`, optional
        Raise a warning if the system is not electrically neutral when
        set to ``True`` (default).

    """
    _so_name = "Coulomb::CoulombMMM1D"
    _so_creation_policy = "GLOBAL"

    def default_params(self):
        return {"far_switch_radius": -1.,
                "verbose": True,
                "timings": 15,
                "check_neutrality": True}

    def required_keys(self):
        return {"prefactor", "maxPWerror"}


@script_interface_register
class MMM1DGPU(ElectrostaticInteraction):
    """
    Electrostatics solver with GPU support for systems with one periodic
    direction. See :ref:`MMM1D on GPU` for more details.

    Parameters
    ----------
    prefactor : :obj:`float`
        Electrostatics prefactor (see :eq:`coulomb_prefactor`).
    maxWPerror : :obj:`float`
        Maximal pairwise error.
    far_switch_radius : :obj:`float`, optional
        Radius where near-field and far-field calculation are switched
    bessel_cutoff : :obj:`int`, optional
    timings : :obj:`int`, optional
        Number of force calculations during tuning.
    check_neutrality : :obj:`bool`, optional
        Raise a warning if the system is not electrically neutral when
        set to ``True`` (default).
    """
    _so_name = "Coulomb::CoulombMMM1DGpu"
    _so_creation_policy = "GLOBAL"
    _so_features = ("MMM1D_GPU",)

    def default_params(self):
        return {"far_switch_radius": -1.,
                "bessel_cutoff": -1,
                "check_neutrality": True}

    def required_keys(self):
        return {"prefactor", "maxPWerror"}


@script_interface_register
class Scafacos(ElectrostaticInteraction):

    """
    Calculate the Coulomb interaction using the ScaFaCoS library.
    See :ref:`ScaFaCoS electrostatics` for more details.

    Parameters
    ----------
    prefactor : :obj:`float`
        Coulomb prefactor as defined in :eq:`coulomb_prefactor`.
    method_name : :obj:`str`
        Name of the ScaFaCoS method to use.
    method_params : :obj:`dict`
        Dictionary containing the method-specific parameters.

    Methods
    -------
    get_available_methods()
        List long-range methods available in the ScaFaCoS library.

    set_near_field_delegation()
        Choose whether to delegate short-range calculation to ESPResSo
        (this is the default when the method supports it) or ScaFaCos.

        Parameters
        ----------
        delegate : :obj:`bool`
            Delegate to ESPResSo if ``True`` and the method supports it.

    get_near_field_delegation()
        Find whether the short-range calculation is delegated to ESPResSo
        (this is the default when the method supports it) or ScaFaCos.

        Returns
        -------
        delegate : :obj:`bool`
            Delegate to ESPResSo if ``True`` and the method supports it,
            ``False`` if delegated to ScaFaCoS or the method doesn't have a
            short-range kernel.

    """
    _so_name = "Coulomb::CoulombScafacos"
    _so_creation_policy = "GLOBAL"
    _so_features = ("ELECTROSTATICS", "SCAFACOS")
    _so_bind_methods = ElectrostaticInteraction._so_bind_methods + \
        ("get_available_methods",
         "get_near_field_delegation",
         "set_near_field_delegation")

    def default_params(self):
        return {"check_neutrality": True}

    def required_keys(self):
        return {"method_name", "method_params", "prefactor"}
