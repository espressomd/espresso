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

import abc
import enum
from . import code_features
from .script_interface import ScriptObjectMap, ScriptInterfaceHelper, script_interface_register


class NonBondedInteraction(ScriptInterfaceHelper, metaclass=abc.ABCMeta):
    """
    Represents an instance of a non-bonded interaction, such as Lennard-Jones.

    Methods
    -------
    deactivate()
        Reset parameters for the interaction.

    """
    _so_bind_methods = ("deactivate",)

    def __init__(self, **kwargs):
        code_features.assert_features(self.__class__.__dict__["_so_feature"])
        if "sip" in kwargs:
            super().__init__(**kwargs)
        else:
            params = self.default_params()
            params.update(kwargs)
            super().__init__(**params)

    def __str__(self):
        return f'{self.__class__.__name__}({self.get_params()})'

    def set_params(self, **kwargs):
        """Set new parameters.

        """
        params = self.default_params()
        params.update(kwargs)

        err_msg = f"setting {self.__class__.__name__} raised an error"
        self.call_method("set_params", handle_errors_message=err_msg, **params)

    def __reduce__(self):
        return (NonBondedInteraction._restore_object,
                (self.__class__, self.get_params()))

    @classmethod
    def _restore_object(cls, derived_class, kwargs):
        return derived_class(**kwargs)

    @abc.abstractmethod
    def default_params(self):
        pass


@script_interface_register
class LennardJonesInteraction(NonBondedInteraction):
    """
    Standard 6-12 Lennard-Jones potential.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        epsilon : :obj:`float`
            Magnitude of the interaction.
        sigma : :obj:`float`
            Interaction length scale.
        cutoff : :obj:`float`
            Cutoff distance of the interaction.
        shift : :obj:`float` or :obj:`str` {'auto'}
            Constant shift of the potential. If ``'auto'``, a default value
            is computed from ``sigma`` and ``cutoff``. The LJ potential
            will be shifted by :math:`4\\epsilon\\cdot\\text{shift}`.
        offset : :obj:`float`, optional
            Offset distance of the interaction.
        min : :obj:`float`, optional
            Restricts the interaction to a minimal distance.

    """

    _so_name = "Interactions::InteractionLJ"
    _so_feature = "LENNARD_JONES"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {"offset": 0., "min": 0.}


@script_interface_register
class WCAInteraction(NonBondedInteraction):
    """
    Standard 6-12 Weeks-Chandler-Andersen potential.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        epsilon : :obj:`float`
            Magnitude of the interaction.
        sigma : :obj:`float`
            Interaction length scale.

    """

    _so_name = "Interactions::InteractionWCA"
    _so_feature = "WCA"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {}

    @property
    def cutoff(self):
        return self.call_method("get_cutoff")


@script_interface_register
class GenericLennardJonesInteraction(NonBondedInteraction):
    """
    Generalized Lennard-Jones potential.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        epsilon : :obj:`float`
            Magnitude of the interaction.
        sigma : :obj:`float`
            Interaction length scale.
        cutoff : :obj:`float`
            Cutoff distance of the interaction.
        shift : :obj:`float` or :obj:`str` {'auto'}
            Constant shift of the potential. If ``'auto'``, a default value
            is computed from the other parameters. The LJ potential
            will be shifted by :math:`\\epsilon\\cdot\\text{shift}`.
        offset : :obj:`float`
            Offset distance of the interaction.
        e1 : :obj:`int`
            Exponent of the repulsion term.
        e2 : :obj:`int`
            Exponent of the attraction term.
        b1 : :obj:`float`
            Prefactor of the repulsion term.
        b2 : :obj:`float`
            Prefactor of the attraction term.
        delta : :obj:`float`, optional
            ``LJGEN_SOFTCORE`` parameter delta. Allows control over how
            smoothly the potential drops to zero as lambda approaches zero.
        lam : :obj:`float`, optional
            ``LJGEN_SOFTCORE`` parameter lambda. Tune the strength of the
            interaction.

    """

    _so_name = "Interactions::InteractionLJGen"
    _so_feature = "LENNARD_JONES_GENERIC"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        if code_features.has_features("LJGEN_SOFTCORE"):
            return {"delta": 0., "lam": 1.}
        return {}


@script_interface_register
class LennardJonesCosInteraction(NonBondedInteraction):
    """Lennard-Jones cosine interaction.

    Methods
    -------
    set_params()
        Set or update parameters for the interaction.
        Parameters marked as required become optional once the
        interaction has been activated for the first time;
        subsequent calls to this method update the existing values.

        Parameters
        ----------
        epsilon : :obj:`float`
            Magnitude of the interaction.
        sigma : :obj:`float`
            Interaction length scale.
        cutoff : :obj:`float`
            Cutoff distance of the interaction.
        offset : :obj:`float`, optional
            Offset distance of the interaction.
    """

    _so_name = "Interactions::InteractionLJcos"
    _so_feature = "LJCOS"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {"offset": 0.}


@script_interface_register
class LennardJonesCos2Interaction(NonBondedInteraction):
    """Second variant of the Lennard-Jones cosine interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        epsilon : :obj:`float`
            Magnitude of the interaction.
        sigma : :obj:`float`
            Interaction length scale.
        offset : :obj:`float`, optional
            Offset distance of the interaction.
        width : :obj:`float`
            Width of interaction.
    """

    _so_name = "Interactions::InteractionLJcos2"
    _so_feature = "LJCOS2"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {"offset": 0.}

    @property
    def cutoff(self):
        return self.call_method("get_cutoff")


@script_interface_register
class HatInteraction(NonBondedInteraction):
    """Hat interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        F_max : :obj:`float`
            Magnitude of the interaction.
        cutoff : :obj:`float`
            Cutoff distance of the interaction.

    """

    _so_name = "Interactions::InteractionHat"
    _so_feature = "HAT"

    def default_params(self):
        return {}


@script_interface_register
class GayBerneInteraction(NonBondedInteraction):
    """Gay--Berne interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        eps : :obj:`float`
            Potential well depth.
        sig : :obj:`float`
            Interaction range.
        cut : :obj:`float`
            Cutoff distance of the interaction.
        k1 : :obj:`float` or :obj:`str`
            Molecular elongation.
        k2 : :obj:`float`, optional
            Ratio of the potential well depths for the side-by-side
            and end-to-end configurations.
        mu : :obj:`float`, optional
            Adjustable exponent.
        nu : :obj:`float`, optional
            Adjustable exponent.

    """

    _so_name = "Interactions::InteractionGayBerne"
    _so_feature = "GAY_BERNE"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {}


@script_interface_register
class TabulatedNonBonded(NonBondedInteraction):
    """Tabulated interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        min : :obj:`float`,
            The minimal interaction distance.
        max : :obj:`float`,
            The maximal interaction distance.
        energy: array_like of :obj:`float`
            The energy table.
        force: array_like of :obj:`float`
            The force table.

    """

    _so_name = "Interactions::InteractionTabulated"
    _so_feature = "TABULATED"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {}

    @property
    def cutoff(self):
        return self.call_method("get_cutoff")


@script_interface_register
class DPDInteraction(NonBondedInteraction):
    """DPD interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        weight_function : :obj:`int`, {0, 1}
            The distance dependence of the parallel part,
            either 0 (constant) or 1 (linear)
        gamma : :obj:`float`
            Friction coefficient of the parallel part
        k : :obj:`float`
            Exponent in the modified weight function
        r_cut : :obj:`float`
            Cutoff of the parallel part
        trans_weight_function : :obj:`int`, {0, 1}
            The distance dependence of the orthogonal part,
            either 0 (constant) or 1 (linear)
        trans_gamma : :obj:`float`
            Friction coefficient of the orthogonal part
        trans_r_cut : :obj:`float`
            Cutoff of the orthogonal part

    """

    _so_name = "Interactions::InteractionDPD"
    _so_feature = "DPD"

    def default_params(self):
        return {
            "weight_function": 0,
            "gamma": 0.0,
            "k": 1.0,
            "r_cut": -1.0,
            "trans_weight_function": 0,
            "trans_gamma": 0.0,
            "trans_r_cut": -1.0}


@script_interface_register
class SmoothStepInteraction(NonBondedInteraction):
    """Smooth step interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        d : :obj:`float`
            Short range repulsion parameter.
        n : :obj:`int`, optional
            Exponent of short range repulsion.
        eps : :obj:`float`
            Magnitude of the second (soft) repulsion.
        k0 : :obj:`float`, optional
            Exponential factor in second (soft) repulsion.
        sig : :obj:`float`, optional
            Length scale of second (soft) repulsion.
        cutoff : :obj:`float`
            Cutoff distance of the interaction.

    """

    _so_name = "Interactions::InteractionSmoothStep"
    _so_feature = "SMOOTH_STEP"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {"n": 10, "k0": 0., "sig": 0.}


@script_interface_register
class BMHTFInteraction(NonBondedInteraction):
    """BMHTF interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        a : :obj:`float`
            Magnitude of exponential part of the interaction.
        b : :obj:`float`
            Exponential factor of the interaction.
        c : :obj:`float`
            Magnitude of the term decaying with the sixth power of r.
        d : :obj:`float`
            Magnitude of the term decaying with the eighth power of r.
        sig : :obj:`float`
            Shift in the exponent.
        cutoff : :obj:`float`
            Cutoff distance of the interaction.

    """

    _so_name = "Interactions::InteractionBMHTF"
    _so_feature = "BMHTF_NACL"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {}


@script_interface_register
class MorseInteraction(NonBondedInteraction):
    """Morse interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        eps : :obj:`float`
            The magnitude of the interaction.
        alpha : :obj:`float`
            Stiffness of the Morse interaction.
        rmin : :obj:`float`
            Distance of potential minimum
        cutoff : :obj:`float`, optional
            Cutoff distance of the interaction.

    """

    _so_name = "Interactions::InteractionMorse"
    _so_feature = "MORSE"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {"cutoff": 0.}


@script_interface_register
class BuckinghamInteraction(NonBondedInteraction):
    """Buckingham interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        a : :obj:`float`
            Magnitude of the exponential part of the interaction.
        b : :obj:`float`, optional
            Exponent of the exponential part of the interaction.
        c : :obj:`float`
            Prefactor of term decaying with the sixth power of distance.
        d : :obj:`float`
            Prefactor of term decaying with the fourth power of distance.
        discont : :obj:`float`
            Distance below which the potential is linearly continued.
        cutoff : :obj:`float`
            Cutoff distance of the interaction.
        shift: :obj:`float`, optional
            Constant potential shift.

    """

    _so_name = "Interactions::InteractionBuckingham"
    _so_feature = "BUCKINGHAM"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {"b": 0., "shift": 0.}


@script_interface_register
class SoftSphereInteraction(NonBondedInteraction):
    """Soft sphere interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        a : :obj:`float`
            Magnitude of the interaction.
        n : :obj:`float`
            Exponent of the power law.
        cutoff : :obj:`float`
            Cutoff distance of the interaction.
        offset : :obj:`float`, optional
            Offset distance of the interaction.

    """

    _so_name = "Interactions::InteractionSoftSphere"
    _so_feature = "SOFT_SPHERE"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {"offset": 0.}


@script_interface_register
class HertzianInteraction(NonBondedInteraction):
    """Hertzian interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        eps : :obj:`float`
            Magnitude of the interaction.
        sig : :obj:`float`
            Interaction length scale.

    """

    _so_name = "Interactions::InteractionHertzian"
    _so_feature = "HERTZIAN"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {}


@script_interface_register
class GaussianInteraction(NonBondedInteraction):
    """Gaussian interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        eps : :obj:`float`
            Overlap energy epsilon.
        sig : :obj:`float`
            Variance sigma of the Gaussian interaction.
        cutoff : :obj:`float`
            Cutoff distance of the interaction.

    """

    _so_name = "Interactions::InteractionGaussian"
    _so_feature = "GAUSSIAN"

    def default_params(self):
        """Python dictionary of default parameters.

        """
        return {}


@script_interface_register
class TholeInteraction(NonBondedInteraction):
    """Thole interaction.

    Methods
    -------
    set_params()
        Set new parameters for the interaction.

        Parameters
        ----------
        scaling_coeff : :obj:`float`
            The factor used in the Thole damping function between
            polarizable particles i and j. Usually calculated by
            the polarizabilities :math:`\\alpha_i`, :math:`\\alpha_j`
            and damping parameters :math:`a_i`, :math:`a_j` via
            :math:`s_{ij} = \\frac{(a_i+a_j)/2}{((\\alpha_i\\cdot\\alpha_j)^{1/2})^{1/3}}`
        q1q2: :obj:`float`
            Charge factor of the involved charges. Has to be set because
            it acts only on the portion of the Drude core charge that is
            associated to the dipole of the atom. For charged, polarizable
            atoms that charge is not equal to the particle charge property.

    """

    _so_name = "Interactions::InteractionThole"
    _so_feature = "THOLE"

    def default_params(self):
        return {}


@script_interface_register
class NonBondedInteractionHandle(ScriptInterfaceHelper):

    """
    Provides access to all non-bonded interactions between two particle types.

    """
    _so_name = "Interactions::NonBondedInteractionHandle"

    def __getattr__(self, key):
        obj = super().__getattr__(key)
        return globals()[obj.__class__.__name__](
            _types=self.call_method("get_types"), **obj.get_params())

    def _serialize(self):
        serialized = []
        for name, obj in self.get_params().items():
            serialized.append((name, obj.__reduce__()[1]))
        return serialized

    def reset(self):
        for key in self._valid_parameters():
            getattr(self, key).deactivate()

    @classmethod
    def _restore_object(cls, types, kwargs):
        objects = {}
        for name, (obj_class, obj_params) in kwargs:
            objects[name] = obj_class(**obj_params)
        return NonBondedInteractionHandle(_types=types, **objects)

    def __reduce__(self):
        return (NonBondedInteractionHandle._restore_object,
                (self.call_method("get_types"), self._serialize()))


@script_interface_register
class NonBondedInteractions(ScriptInterfaceHelper):
    """
    Access to non-bonded interaction parameters via ``[i,j]``, where ``i, j``
    are particle types. Returns a :class:`NonBondedInteractionHandle` object.

    Methods
    -------
    reset()
        Reset all interaction parameters to their default values.

    """
    _so_name = "Interactions::NonBondedInteractions"
    _so_creation_policy = "GLOBAL"
    _so_bind_methods = ("reset",)

    def keys(self):
        return [tuple(x) for x in self.call_method("keys")]

    def __getitem__(self, key):
        self.call_method("check_key", key=key)
        return NonBondedInteractionHandle(_types=key)

    def __setitem__(self, key, value):
        self.call_method("insert", key=key, object=value)

    def __getstate__(self):
        n_types = self.call_method("get_n_types")
        state = []
        for i in range(n_types):
            for j in range(i, n_types):
                handle = NonBondedInteractionHandle(_types=(i, j))
                state.append(((i, j), handle._serialize()))
        return {"state": state}

    def __setstate__(self, params):
        for types, kwargs in params["state"]:
            obj = NonBondedInteractionHandle._restore_object(types, kwargs)
            self.call_method("insert", key=types, object=obj)

    @classmethod
    def _restore_object(cls, so_callback, so_callback_args, state):
        so = so_callback(*so_callback_args)
        so.__setstate__(state)
        return so

    def __reduce__(self):
        so_callback, (so_name, so_bytestring) = super().__reduce__()
        return (NonBondedInteractions._restore_object,
                (so_callback, (so_name, so_bytestring), self.__getstate__()))


class BONDED_IA(enum.IntEnum):
    NONE = 0
    FENE = enum.auto()
    HARMONIC = enum.auto()
    QUARTIC = enum.auto()
    BONDED_COULOMB = enum.auto()
    BONDED_COULOMB_SR = enum.auto()
    ANGLE_HARMONIC = enum.auto()
    ANGLE_COSINE = enum.auto()
    ANGLE_COSSQUARE = enum.auto()
    DIHEDRAL = enum.auto()
    TABULATED_DISTANCE = enum.auto()
    TABULATED_ANGLE = enum.auto()
    TABULATED_DIHEDRAL = enum.auto()
    THERMALIZED_DIST = enum.auto()
    RIGID_BOND = enum.auto()
    IBM_TRIEL = enum.auto()
    IBM_VOLUME_CONSERVATION = enum.auto()
    IBM_TRIBEND = enum.auto()
    OIF_GLOBAL_FORCES = enum.auto()
    OIF_LOCAL_FORCES = enum.auto()
    VIRTUAL_BOND = enum.auto()


class BondedInteraction(ScriptInterfaceHelper, metaclass=abc.ABCMeta):

    """
    Base class for bonded interactions.

    Either called with an interaction id, in which case the interaction
    will represent the bonded interaction as it is defined in ESPResSo core,
    or called with keyword arguments describing a new interaction.

    """

    _so_name = "Interactions::BondedInteraction"
    _so_creation_policy = "GLOBAL"

    def __init__(self, **kwargs):
        feature = self.__class__.__dict__.get("_so_feature")
        if feature is not None:
            code_features.assert_features(feature)

        if "sip" not in kwargs:
            if "bond_id" in kwargs:
                # create a new script interface object for a bond that already
                # exists in the core via its id (BondedInteractions getter and
                # checkpointing constructor #1)
                bond_id = kwargs["bond_id"]
                super().__init__(bond_id=bond_id)
                # Check if the bond type in ESPResSo core matches this class
                if self.call_method("get_zero_based_type",
                                    bond_id=bond_id) != self._type_number:
                    raise RuntimeError(
                        f"The bond with id {bond_id} is not defined as a "
                        f"{self._type_number.name} bond in the ESPResSo core.")
                self._bond_id = bond_id
                self._ctor_params = self.get_params()
            else:
                # create a new script interface object from bond parameters
                # (normal bond creation and checkpointing constructor #2)
                params = self.get_default_params()
                params.update(kwargs)
                super().__init__(**params)
                self._ctor_params = params
                self._bond_id = -1
        else:
            # create a new bond based on a bond in the script interface
            # (checkpointing constructor #3)
            super().__init__(**kwargs)
            self._bond_id = -1
            self._ctor_params = self.get_params()

    def __reduce__(self):
        if self._bond_id != -1:
            # checkpointing constructor #1
            return (BondedInteraction._restore_object,
                    (self.__class__, {"bond_id": self._bond_id}))
        else:
            # checkpointing constructor #2
            return (BondedInteraction._restore_object,
                    (self.__class__, self._serialize()))

    def _serialize(self):
        return self._ctor_params.copy()

    @classmethod
    def _restore_object(cls, derived_class, kwargs):
        return derived_class(**kwargs)

    def __setattr__(self, attr, value):
        super().__setattr__(attr, value)

    @property
    def params(self):
        return self.get_params()

    @params.setter
    def params(self, p):
        raise RuntimeError("Bond parameters are immutable.")

    def __str__(self):
        return f'{self.__class__.__name__}({self._ctor_params})'

    @abc.abstractmethod
    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        pass

    def __repr__(self):
        return f'<{self}>'

    def __eq__(self, other):
        return self.__class__ == other.__class__ and self.call_method(
            "is_same_bond", bond=other)

    def __ne__(self, other):
        return not self.__eq__(other)


@script_interface_register
class FeneBond(BondedInteraction):

    """
    FENE bond.

    Parameters
    ----------
    k : :obj:`float`
        Magnitude of the bond interaction.
    d_r_max : :obj:`float`
        Maximum stretch and compression length of the bond.
    r_0 : :obj:`float`, optional
        Equilibrium bond length.

    """

    _so_name = "Interactions::FeneBond"
    _type_number = BONDED_IA.FENE

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {"r_0": 0.}


@script_interface_register
class HarmonicBond(BondedInteraction):

    """
    Harmonic bond.

    Parameters
    ----------
    k : :obj:`float`
        Magnitude of the bond interaction.
    r_0 : :obj:`float`
        Equilibrium bond length.
    r_cut : :obj:`float`, optional
        Maximum distance beyond which the bond is considered broken.

    """

    _so_name = "Interactions::HarmonicBond"
    _type_number = BONDED_IA.HARMONIC

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {"r_cut": 0.}


@script_interface_register
class BondedCoulomb(BondedInteraction):

    """
    Bonded Coulomb bond.

    Parameters
    ----------

    prefactor : :obj:`float`
        Coulomb prefactor of the bonded Coulomb interaction.
    """

    _so_name = "Interactions::BondedCoulomb"
    _so_feature = "ELECTROSTATICS"
    _type_number = BONDED_IA.BONDED_COULOMB

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class BondedCoulombSRBond(BondedInteraction):

    """
    Bonded Coulomb short range bond. Calculates the short range part of
    Coulomb interactions.

    Parameters
    ----------

    q1q2 : :obj:`float`
        Charge factor of the involved particle pair. Note the
        particle charges are used to allow e.g. only partial subtraction
        of the involved charges.
    """

    _so_name = "Interactions::BondedCoulombSR"
    _so_feature = "ELECTROSTATICS"
    _type_number = BONDED_IA.BONDED_COULOMB_SR

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class ThermalizedBond(BondedInteraction):

    """
    Thermalized bond.

    Parameters
    ----------

    temp_com : :obj:`float`
        Temperature of the Langevin thermostat for the center of mass of the
        particle pair.
    gamma_com: :obj:`float`
        Friction coefficient of the Langevin thermostat for the center of mass
        of the particle pair.
    temp_distance: :obj:`float`
        Temperature of the Langevin thermostat for the distance vector
        of the particle pair.
    gamma_distance: :obj:`float`
        Friction coefficient of the Langevin thermostat for the
        distance vector of the particle pair.
    r_cut: :obj:`float`, optional
        Maximum distance beyond which the bond is considered broken.
    seed : :obj:`int`
        Seed of the philox RNG. Must be positive.
        Required for the first thermalized bond in the system. Subsequent
        thermalized bonds don't need a seed; if one is provided nonetheless,
        it will overwrite the seed of all previously defined thermalized bonds,
        even if the new bond is not added to the system.

    """

    _so_name = "Interactions::ThermalizedBond"
    _type_number = BONDED_IA.THERMALIZED_DIST

    def __init__(self, *args, **kwargs):
        if kwargs and "sip" not in kwargs:
            kwargs["rng_state"] = kwargs.get("rng_state")
        super().__init__(*args, **kwargs)

    def _serialize(self):
        params = self._ctor_params.copy()
        params["rng_state"] = self.call_method("get_rng_state")
        return params

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {"r_cut": 0., "seed": None}


@script_interface_register
class RigidBond(BondedInteraction):

    """
    Rigid bond.

    Parameters
    ----------
    r : :obj:`float`
        Bond length.
    ptol : :obj:`float`, optional
        Tolerance for positional deviations.
    vtop : :obj:`float`, optional
        Tolerance for velocity deviations.

    """

    _so_name = "Interactions::RigidBond"
    _so_feature = "BOND_CONSTRAINT"
    _type_number = BONDED_IA.RIGID_BOND

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        # TODO rationality of Default Parameters has to be checked
        return {"ptol": 0.001, "vtol": 0.001}


@script_interface_register
class Dihedral(BondedInteraction):

    """
    Dihedral potential with phase shift.

    Parameters
    ----------
    mult : :obj:`int`
        Multiplicity of the potential (number of minima).
    bend : :obj:`float`
        Bending constant.
    phase : :obj:`float`
        Angle of the first local minimum in radians.

    """

    _so_name = "Interactions::DihedralBond"
    _type_number = BONDED_IA.DIHEDRAL

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class TabulatedDistance(BondedInteraction):

    """
    Tabulated bond length.

    Parameters
    ----------

    min : :obj:`float`
        The minimal interaction distance.
    max : :obj:`float`
        The maximal interaction distance.
    energy: array_like of :obj:`float`
        The energy table.
    force: array_like of :obj:`float`
        The force table.

    """

    _so_name = "Interactions::TabulatedDistanceBond"
    _so_feature = "TABULATED"
    _type_number = BONDED_IA.TABULATED_DISTANCE

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class TabulatedAngle(BondedInteraction):

    """
    Tabulated bond angle.

    Parameters
    ----------

    energy: array_like of :obj:`float`
        The energy table for the range :math:`0-\\pi`.
    force: array_like of :obj:`float`
        The force table for the range :math:`0-\\pi`.

    """

    _so_name = "Interactions::TabulatedAngleBond"
    _so_feature = "TABULATED"
    _type_number = BONDED_IA.TABULATED_ANGLE

    pi = 3.14159265358979

    def __init__(self, *args, **kwargs):
        if len(args) == 0 and "sip" not in kwargs:
            kwargs.update({"min": 0., "max": self.pi})
        super().__init__(*args, **kwargs)

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class TabulatedDihedral(BondedInteraction):

    """
    Tabulated bond dihedral.

    Parameters
    ----------

    energy: array_like of :obj:`float`
        The energy table for the range :math:`0-2\\pi`.
    force: array_like of :obj:`float`
        The force table for the range :math:`0-2\\pi`.

    """

    _so_name = "Interactions::TabulatedDihedralBond"
    _so_feature = "TABULATED"
    _type_number = BONDED_IA.TABULATED_DIHEDRAL

    pi = 3.14159265358979

    def __init__(self, *args, **kwargs):
        if len(args) == 0 and "sip" not in kwargs:
            kwargs.update({"min": 0., "max": 2. * self.pi})
        super().__init__(*args, **kwargs)

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class Virtual(BondedInteraction):

    """
    Virtual bond.
    """

    _so_name = "Interactions::VirtualBond"
    _type_number = BONDED_IA.VIRTUAL_BOND

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class AngleHarmonic(BondedInteraction):

    """
    Bond-angle-dependent harmonic potential.

    Parameters
    ----------
    phi0 : :obj:`float`
        Equilibrium bond angle in radians.
    bend : :obj:`float`
        Magnitude of the bond interaction.

    """

    _so_name = "Interactions::AngleHarmonicBond"
    _type_number = BONDED_IA.ANGLE_HARMONIC

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class AngleCosine(BondedInteraction):

    """
    Bond-angle-dependent cosine potential.

    Parameters
    ----------
    phi0 : :obj:`float`
        Equilibrium bond angle in radians.
    bend : :obj:`float`
        Magnitude of the bond interaction.

    """

    _so_name = "Interactions::AngleCosineBond"
    _type_number = BONDED_IA.ANGLE_COSINE

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class AngleCossquare(BondedInteraction):

    """
    Bond-angle-dependent cosine squared potential.

    Parameters
    ----------
    phi0 : :obj:`float`
        Equilibrium bond angle in radians.
    bend : :obj:`float`
        Magnitude of the bond interaction.

    """

    _so_name = "Interactions::AngleCossquareBond"
    _type_number = BONDED_IA.ANGLE_COSSQUARE

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class IBM_Triel(BondedInteraction):

    """
    IBM Triel bond.

    See Figure C.1 in :cite:`kruger12a`.

    Parameters
    ----------
    ind1, ind2, ind3 : :obj:`int`
        First, second and third bonding partner. Used for
        initializing reference state
    k1 : :obj:`float`
        Shear elasticity for Skalak and Neo-Hookean
    k2 : :obj:`float`, optional
        Area resistance for Skalak
    maxDist : :obj:`float`
        Gives an error if an edge becomes longer than maxDist
    elasticLaw : :obj:`str`, {'NeoHookean', 'Skalak'}
        Type of elastic bond

    """

    _so_name = "Interactions::IBMTriel"
    _type_number = BONDED_IA.IBM_TRIEL

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {"k2": 0}


@script_interface_register
class IBM_Tribend(BondedInteraction):

    """
    IBM Tribend bond.

    See Figure C.2 in :cite:`kruger12a`.

    Parameters
    ----------
    ind1, ind2, ind3, ind4 : :obj:`int`
        First, second, third and fourth bonding partner. Used for
        initializing reference state
    kb : :obj:`float`
        Bending modulus
    refShape : :obj:`str`, optional, {'Flat', 'Initial'}
        Reference shape, default is ``'Flat'``

    """

    _so_name = "Interactions::IBMTribend"
    _type_number = BONDED_IA.IBM_TRIBEND

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {"refShape": "Flat"}


@script_interface_register
class IBM_VolCons(BondedInteraction):

    """
    IBM volume conservation bond.

    See Figure C.3 in :cite:`kruger12a`.

    Parameters
    ----------
    softID : :obj:`int`
        Used to identify the object to which this bond belongs. Each object
        (cell) needs its own ID. For performance reasons, it is best to
        start from ``softID=0`` and increment by 1 for each subsequent bond.
    kappaV : :obj:`float`
        Modulus for volume force

    Methods
    -------
    current_volume()
        Query the current volume of the soft object associated to this bond.
        The volume is initialized once all :class:`IBM_Triel` bonds have
        been added and the forces have been recalculated.

    """

    _so_name = "Interactions::IBMVolCons"
    _so_bind_methods = ("current_volume",)
    _type_number = BONDED_IA.IBM_VOLUME_CONSERVATION

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class OifGlobalForces(BondedInteraction):

    """
    Characterize the distribution of the force of the global mesh deformation
    onto individual vertices of the mesh.

    Part of the :ref:`Object-in-fluid` method.

    Parameters
    ----------
    A0_g : :obj:`float`
        Relaxed area of the mesh
    ka_g : :obj:`float`
        Area coefficient
    V0 : :obj:`float`
        Relaxed volume of the mesh
    kv : :obj:`float`
        Volume coefficient

    """

    _so_name = "Interactions::OifGlobalForcesBond"
    _type_number = BONDED_IA.OIF_GLOBAL_FORCES

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class OifLocalForces(BondedInteraction):

    """
    Characterize the deformation of two triangles sharing an edge.

    Part of the :ref:`Object-in-fluid` method.

    Parameters
    ----------
    r0 : :obj:`float`
        Equilibrium bond length of triangle edges
    ks : :obj:`float`
        Non-linear stretching coefficient of triangle edges
    kslin : :obj:`float`
        Linear stretching coefficient of triangle edges
    phi0 : :obj:`float`
        Equilibrium angle between the two triangles
    kb : :obj:`float`
        Bending coefficient for the angle between the two triangles
    A01 : :obj:`float`
        Equilibrium surface of the first triangle
    A02 : :obj:`float`
        Equilibrium surface of the second triangle
    kal : :obj:`float`
        Stretching coefficient of a triangle surface
    kvisc : :obj:`float`
        Viscous coefficient of the triangle vertices

    """

    _so_name = "Interactions::OifLocalForcesBond"
    _type_number = BONDED_IA.OIF_LOCAL_FORCES

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class QuarticBond(BondedInteraction):

    """
    Quartic bond.

    Parameters
    ----------
    k0 : :obj:`float`
        Magnitude of the square term.
    k1 : :obj:`float`
        Magnitude of the fourth order term.
    r : :obj:`float`
        Equilibrium bond length.
    r_cut : :obj:`float`
        Maximum interaction length.
    """

    _so_name = "Interactions::QuarticBond"
    _type_number = BONDED_IA.QUARTIC

    def get_default_params(self):
        """Gets default values of optional parameters.

        """
        return {}


@script_interface_register
class BondedInteractions(ScriptObjectMap):

    """
    Represents the bonded interactions list.

    Individual interactions can be accessed using ``BondedInteractions[i]``,
    where ``i`` is the bond id.

    Methods
    -------
    remove()
        Remove a bond from the list.
        This is a no-op if the bond does not exist.

        Parameters
        ----------
        bond_id : :obj:`int`

    clear()
        Remove all bonds.

    """

    _so_name = "Interactions::BondedInteractions"
    _so_creation_policy = "GLOBAL"
    _bond_classes = {
        cls._type_number: cls for cls in globals().values()
        if isinstance(cls, type) and issubclass(cls, BondedInteraction) and cls != BondedInteraction
    }

    def add(self, *args):
        """
        Add a bond to the list.

        Parameters
        ----------
        bond: :class:`espressomd.interactions.BondedInteraction`
            Either a bond object...

        """

        if len(args) == 1 and isinstance(args[0], BondedInteraction):
            bonded_ia = args[0]
        else:
            raise TypeError("A BondedInteraction object needs to be passed.")
        bond_id = self._insert_bond(None, bonded_ia)
        return bond_id

    def __getitem__(self, bond_id):
        if self.call_method('has_bond', bond_id=bond_id):
            bond_obj = self.call_method('get_bond', bond_id=bond_id)
            bond_obj._bond_id = bond_id
            return bond_obj

        # Find out the type of the interaction from ESPResSo
        bond_type = self.call_method("get_zero_based_type", bond_id=bond_id)

        # Check if the bonded interaction exists in ESPResSo core
        if bond_type == BONDED_IA.NONE:
            raise ValueError(f"The bond with id {bond_id} is not yet defined.")

        # Find the appropriate class representing such a bond
        bond_class = self._bond_classes[bond_type]

        # Create a new script interface object (i.e. a copy of the shared_ptr)
        # which links to the bonded interaction object
        return bond_class(bond_id=bond_id)

    def __setitem__(self, bond_id, bond_obj):
        self._insert_bond(bond_id, bond_obj)

    def _insert_bond(self, bond_id, bond_obj):
        """
        Inserts a new bond. If a ``bond_id`` is given, the bond is inserted at
        that id. If no id is given, a new id is generated.

        Bonds can only be overwritten if the new bond is of the same type as the
        old one, e.g. a :class:`~espressomd.interactions.FeneBond` bond can only
        be overwritten by another :class:`~espressomd.interactions.FeneBond` bond.
        """

        # Validate arguments
        if not isinstance(bond_obj, BondedInteraction):
            raise ValueError(
                "Only subclasses of BondedInteraction can be assigned.")

        # Send the script interface object pointer to the core
        if bond_id is None:
            bond_id = self.call_method("insert", object=bond_obj)
        else:
            # Throw error if attempting to overwrite a bond of different type
            if self.call_method("contains", key=bond_id):
                old_type = self._bond_classes[
                    self.call_method("get_zero_based_type", bond_id=bond_id)]
                if not isinstance(bond_obj, old_type):
                    raise ValueError(
                        "Bonds can only be overwritten by bonds of equal type.")
            self.call_method("insert", key=bond_id, object=bond_obj)

        # Save the bond id in the BondedInteraction instance
        bond_obj._bond_id = bond_id

        return bond_id

    def __len__(self):
        return self.call_method('get_size')

    # Support iteration over active bonded interactions
    def __iter__(self):
        for bond_id in self.call_method('get_bond_ids'):
            if self.call_method("get_zero_based_type", bond_id=bond_id):
                yield self[bond_id]

    def __getstate__(self):
        params = {}
        for bond_id in self.call_method('get_bond_ids'):
            if self.call_method("get_zero_based_type", bond_id=bond_id):
                obj = self[bond_id]
                if hasattr(obj, "params"):
                    params[bond_id] = (obj._type_number, obj._serialize())
        return params

    def __setstate__(self, params):
        for bond_id, (type_number, bond_params) in params.items():
            self[bond_id] = self._bond_classes[type_number](**bond_params)

    def __reduce__(self):
        so_callback, (so_name, so_bytestring) = super().__reduce__()
        return (BondedInteractions._restore_object,
                (so_callback, (so_name, so_bytestring), self.__getstate__()))

    @classmethod
    def _restore_object(cls, so_callback, so_callback_args, state):
        so = so_callback(*so_callback_args)
        so.__setstate__(state)
        return so
