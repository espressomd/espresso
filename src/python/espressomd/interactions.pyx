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
include "myconfig.pxi"
from . import utils
# Non-bonded interactions

cdef class NonBondedInteraction(object):

    cdef public object _part_types
    cdef object _params

    # init dict to access all user defined nonbonded-inters via
    # user_interactions[type1][type2][parameter]
    user_interactions = {}

    def __init__(self, *args, **kwargs):
        """Represents an instance of a non-bonded interaction, such as lennard jones
        Either called with two particle type id, in which case, the interaction
        will represent the bonded interaction as it is defined in Espresso core
        Or called with keyword arguments describing a new interaction."""

        # Interaction id as argument
        if len(args) == 2 and isinstance(args[0], int) and isinstance(args[1], int):
            self._part_types = args

            # Load the parameters currently set in the Espresso core
            self._params = self._get_params_from_es_core()

        # Or have we been called with keyword args describing the interaction
        elif len(args) == 0:
            # Initialize default values
            self._params = self.default_params()
            self._part_types = [-1, -1]

            # Check if all required keys are given
            for k in self.required_keys():
                if k not in kwargs:
                    raise ValueError(
                        "At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__())

            self._params.update(kwargs)

            # Validation of parameters
            self.validate_params()

        else:
            raise Exception(
                "The constructor has to be called either with two particle type ids (as interger), or with a set of keyword arguments describing a new interaction")

    def is_valid(self):
        """Check, if the data stored in the instance still matches what is in Espresso"""

        # check, if the bond parameters saved in the class still match those
        # saved in Espresso
        temp_params = self._get_params_from_es_core()
        if self._params != temp_params:
            return False

        # If we're still here, the instance is valid
        return True

    def get_params(self):
        """Get interaction parameters"""
        # If this instance refers to an actual interaction defined in the es core, load
        # current parameters from there
        if self._part_types[0] >= 0 and self._part_types[1] >= 0:
            self._params = self._get_params_from_es_core()

        return self._params

    def __str__(self):
        return self.__class__.__name__ + "(" + str(self.get_params()) + ")"

    def set_params(self, **p):
        """Update parameters. Only given """
        # Check, if any key was passed, which is not known
        for k in p.keys():
            if k not in self.valid_keys():
                raise ValueError(
                    "Only the following keys are supported: " + self.valid_keys().__str__())

        # When an interaction is newly activated, all required keys must be
        # given
        if not self.is_active():
            for k in self.required_keys():
                if k not in p:
                    raise ValueError(
                        "At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__())

        # If this instance refers to an interaction defined in the espresso core,
        # load the parameters from there

        if self._part_types[0] >= 0 and self._part_types[1] >= 0:
            self._params = self._get_params_from_es_core()

        # Put in values given by the user
        self._params.update(p)

        if self._part_types[0] >= 0 and self._part_types[1] >= 0:
            self._set_params_in_es_core()

        # update interaction dict when user sets interaction
        if self._part_types[0] not in self.user_interactions:
            self.user_interactions[self._part_types[0]] = {}
        self.user_interactions[self._part_types[0]][self._part_types[1]] = {}
        new_params = self.get_params()
        for p_key in new_params:
            self.user_interactions[self._part_types[0]][
                self._part_types[1]][p_key] = new_params[p_key]
        self.user_interactions[self._part_types[0]][
            self._part_types[1]]['type_name'] = self.type_name()

    def validate_params(self):
        return True

    def __getattribute__(self, name):
        """Every time _set_params_in_es_core is called, the parameter dict is also updated."""
        attr = object.__getattribute__(self, name)
        if hasattr(attr, '__call__') and attr.__name__ == "_set_params_in_es_core":
            def sync_params(*args, **kwargs):
                result = attr(*args, **kwargs)
                self._params.update(self._get_params_from_es_core())
                return result
            return sync_params
        else:
            return attr

    def _get_params_from_es_core(self):
        raise Exception(
            "Subclasses of NonBondedInteraction must define the _get_params_from_es_core() method.")

    def _set_params_in_es_core(self):
        raise Exception(
            "Subclasses of NonBondedInteraction must define the _set_params_in_es_core() method.")

    def default_params(self):
        raise Exception(
            "Subclasses of NonBondedInteraction must define the default_params() method.")

    def is_active(self):
        # If this instance refers to an actual interaction defined in the es core, load
        # current parameters from there
        if self._part_types[0] >= 0 and self._part_types[1] >= 0:
            self._params = self._get_params_from_es_core()
        raise Exception(
            "Subclasses of NonBondedInteraction must define the is_active() method.")

    def type_name(self):
        raise Exception(
            "Subclasses of NonBondedInteraction must define the type_name() method.")

    def valid_keys(self):
        raise Exception(
            "Subclasses of NonBondedInteraction must define the valid_keys() method.")

    def required_keys(self):
        raise Exception(
            "Subclasses of NonBondedInteraction must define the required_keys() method.")

# Lennard Jones

IF LENNARD_JONES == 1:
    cdef class LennardJonesInteraction(NonBondedInteraction):

        def validate_params(self):
            if self._params["epsilon"] < 0:
                raise ValueError("Lennard-Jones eps has to be >=0")
            if self._params["sigma"] < 0:
                raise ValueError("Lennard-Jones sigma has to be >=0")
            if self._params["cutoff"] < 0:
                raise ValueError("Lennard-Jones cutoff has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef ia_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "epsilon": ia_params.LJ_eps,
                "sigma": ia_params.LJ_sig,
                "cutoff": ia_params.LJ_cut,
                "shift": ia_params.LJ_shift,
                "offset": ia_params.LJ_offset,
                "cap": ia_params.LJ_capradius,
                "min": ia_params.LJ_min}

        def is_active(self):
            return (self._params["epsilon"] > 0)

        def set_params(self, **kwargs):
            """ Set parameters for the Lennard-Jones interaction.

            Parameters
            ----------

            epsilon : float
                      The magnitude of the interaction.
            sigma : float
                    Determines the interaction length scale.
            cutoff : float
                     Cutoff distance of the interaction.
            shift : float, string
                    Constant shift of the potential. (4*epsilon*shift).
            offset : float, optional
                     Offset distance of the interaction.
            cap : float, optional
                  If individual force caps are used, determines the distance
                  at which the force is capped.
            min : float, optional
                  Restricts the interaction to a minimal distance.
            """
            super(LennardJonesInteraction, self).set_params(**kwargs)

        def _set_params_in_es_core(self):
            # Handle the case of shift="auto"
            if self._params["shift"] == "auto":
                # Calc shift
                self._params["shift"] = -((self._params["sigma"] / self._params["cutoff"])**12 - (
                    self._params["sigma"] / self._params["cutoff"])**6)

            if lennard_jones_set_params(
                self._part_types[0], self._part_types[1],
                                        self._params["epsilon"],
                                        self._params["sigma"],
                                        self._params["cutoff"],
                                        self._params["shift"],
                                        self._params["offset"],
                                        self._params["cap"],
                                        self._params["min"]):
                raise Exception("Could not set Lennard Jones parameters")

        def default_params(self):
            return {
                "epsilon": 0.,
                "sigma": 0.,
                "cutoff": 0.,
                "shift": 0.,
                "offset": 0.,
                "cap": 0.,
                "min": 0.}

        def type_name(self):
            return "LennardJones"

        def valid_keys(self):
            return "epsilon", "sigma", "cutoff", "shift", "offset", "cap", "min"

        def required_keys(self):
            return "epsilon", "sigma", "cutoff", "shift"

# Lennard Jones

IF GAY_BERNE:
    cdef class GayBerneInteraction(NonBondedInteraction):

        def validate_params(self):
            return True

        def _get_params_from_es_core(self):
            cdef ia_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "eps": ia_params.GB_eps,
                "sig": ia_params.GB_sig,
                "cut": ia_params.GB_cut,
                "k1": ia_params.GB_k1,
                "k2": ia_params.GB_k2,
                "mu": ia_params.GB_mu,
                "nu": ia_params.GB_nu}

        def is_active(self):
            return (self._params["eps"] > 0)

        def set_params(self, **kwargs):
            """ Set parameters for the Lennard-Jones interaction.

            Parameters
            ----------

            eps : float
                  Potential well depth.
            sig : float
                  Interaction range.
            cut : float
                  Cutoff distance of the interaction.
            k1  : float, string
                  Molecular elongation.
            k2  : float, optional
                  Ratio of the potential well depths for the side-by-side
                  and end-to-end configurations.
            mu  : float, optional
                  Adjustable exponent.
            nu  : float, optional
                  Adjustable exponent.
            """
            super(GayBerneInteraction, self).set_params(**kwargs)

        def _set_params_in_es_core(self):
            if gay_berne_set_params(self._part_types[0], self._part_types[1],
                                    self._params["eps"],
                                    self._params["sig"],
                                    self._params["cut"],
                                    self._params["k1"],
                                    self._params["k2"],
                                    self._params["mu"],
                                    self._params["nu"]):
                raise Exception("Could not set Gay Berne parameters")

        def default_params(self):
            return {
                "eps": 0.0,
                "sig": 0.0,
                "cut": 0.0,
                "k1": 0.0,
                "k2": 0.0,
                "mu": 0.0,
                "nu": 0.0}

        def type_name(self):
            return "GayBerne"

        def valid_keys(self):
            return "eps", "sig", "cut", "k1", "k2", "mu", "nu"

        def required_keys(self):
            return "eps", "sig", "cut", "k1", "k2", "mu", "nu"

# Generic Lennard Jones
IF LENNARD_JONES_GENERIC == 1:

    cdef class GenericLennardJonesInteraction(NonBondedInteraction):

        def validate_params(self):
            if self._params["epsilon"] < 0:
                raise ValueError("Generic Lennard-Jones eps has to be >=0")
            if self._params["sigma"] < 0:
                raise ValueError("Generic Lennard-Jones sigma has to be >=0")
            if self._params["cutoff"] < 0:
                raise ValueError("Generic Lennard-Jones cutoff has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef ia_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "epsilon": ia_params.LJGEN_eps,
                "sigma": ia_params.LJGEN_sig,
                "cutoff": ia_params.LJGEN_cut,
                "shift": ia_params.LJGEN_shift,
                "offset": ia_params.LJGEN_offset,
                "e1": ia_params.LJGEN_a1,
                "e2": ia_params.LJGEN_a2,
                "b1": ia_params.LJGEN_b1,
                "b2": ia_params.LJGEN_b2,
                "lam": ia_params.LJGEN_lambda,
                "delta": ia_params.LJGEN_softrad
            }

        def is_active(self):
            return (self._params["epsilon"] > 0)

        def _set_params_in_es_core(self):
            # Handle the case of shift="auto"
            if self._params["shift"] == "auto":
                # Calc shift
                self._params["shift"] = -(self._params["b1"] * (self._params["sigma"] / self._params["cutoff"])**self._params[
                                          "e1"] - self._params["b2"] * (self._params["sigma"] / self._params["cutoff"])**self._params["e2"])
            IF LJGEN_SOFTCORE:
                if ljgen_set_params(self._part_types[0], self._part_types[1],
                                    self._params["epsilon"],
                                    self._params["sigma"],
                                    self._params["cutoff"],
                                    self._params["shift"],
                                    self._params["offset"],
                                    self._params["e1"],
                                    self._params["e2"],
                                    self._params["b1"],
                                    self._params["b2"],
                                    0.0,
                                    self._params["lam"],
                                    self._params["delta"]):
                    raise Exception(
                        "Could not set Generic Lennard Jones parameters")
            ELSE:
                if ljgen_set_params(self._part_types[0], self._part_types[1],
                                    self._params["epsilon"],
                                    self._params["sigma"],
                                    self._params["cutoff"],
                                    self._params["shift"],
                                    self._params["offset"],
                                    self._params["e1"],
                                    self._params["e2"],
                                    self._params["b1"],
                                    self._params["b2"],
                                    0.0):
                    raise Exception(
                        "Could not set Generic Lennard Jones parameters")

        def default_params(self):
            return {
                "epsilon": 0.,
                "sigma": 0.,
                "cutoff": 0.,
                "shift": 0.,
                "offset": 0.,
                "e1": 0,
                "e2": 0,
                "b1": 0.,
                "b2": 0.,
                "delta": 0.,
                "lam": 0.}

        def type_name(self):
            return "GenericLennardJones"

        def set_params(self, **kwargs):
            """
            Set parameters for the generic Lennard-Jones interaction.

            Parameters
            ----------
            epsilon : float
                      The magnitude of the interaction.
            sigma : float
                    Determines the interaction length scale.
            cutoff : float
                     Cutoff distance of the interaction.
            shift : float, string
                    Constant shift of the potential.
            offset : float
                     Offset distance of the interaction.
            e1 : int
                 Exponent of the repulsion term.
            e2 : int
                 Exponent of the attraction term.
            b1 : float
                 Prefactor of the repulsion term.
            b2 : float
                 Prefactor of the attraction term.
            delta : float, optional
                    LJGEN_SOFTCORE parameter. Allows control over how smoothly
                    the potential drops to zero as lambda approaches zero.
            lam : float, optional
                     LJGEN_SOFTCORE parameter lambda. Tune the strength of the
                     interaction.
            """
            super(GenericLennardJonesInteraction, self).set_params(**kwargs)

        def valid_keys(self):
            return "epsilon", "sigma", "cutoff", "shift", "offset", "e1", "e2", "b1", "b2", "delta", "lam"

        def required_keys(self):
            return "epsilon", "sigma", "cutoff", "shift", "offset", "e1", "e2", "b1", "b2"


class NonBondedInteractionHandle(object):

    """Provides access to all Non-bonded interactions between
    two particle types."""

    type1 = -1
    type2 = -1

    # Here, one line per non-bonded ia
    lennard_jones = None
    generic_lennard_jones = None
    tabulated = None
    gay_berne = None

    def __init__(self, _type1, _type2):
        """Takes two particle types as argument"""
        if not (isinstance(_type1, int) and isinstance(_type2, int)):
            raise TypeError("The particle types have to be of type integer.")
        self.type1 = _type1
        self.type2 = _type2

        # Here, add one line for each nonbonded ia
        IF LENNARD_JONES:
            self.lennard_jones = LennardJonesInteraction(_type1, _type2)
        IF LENNARD_JONES_GENERIC:
            self.generic_lennard_jones = GenericLennardJonesInteraction(
                _type1, _type2)
        IF TABULATED == 1:
            self.tabulated = TabulatedNonBonded(_type1, _type2)
        IF GAY_BERNE:
            self.gay_berne = GayBerneInteraction(_type1, _type2)


cdef class NonBondedInteractions(object):

    """Access to non-bonded interaction parameters via [i,j], where i,j are particle
    types. Returns NonBondedInteractionHandle.
    Also: access to force capping
    """

    def __getitem__(self, key):
        if not isinstance(key, tuple):
            raise ValueError(
                "NonBondedInteractions[] expects two particle types as indices.")
        if len(key) != 2 or (not isinstance(key[0], int)) or (not isinstance(key[1], int)):
            raise ValueError(
                "NonBondedInteractions[] expects two particle types as indices.")
        return NonBondedInteractionHandle(key[0], key[1])

    def set_force_cap(self, cap):
        if forcecap_set_params(cap):
            raise Exception("Could not set forcecap")

    def get_force_cap(self):
        return force_cap

    def __getstate__(self):
        # contains info about ALL nonbonded interactions
        odict = NonBondedInteractionHandle(-1, -
                                           1).lennard_jones.user_interactions
        odict['force_cap'] = self.get_force_cap()
        return odict

    def __setstate__(self, odict):
        self.set_force_cap(odict['force_cap'])
        del odict['force_cap']
        for _type1 in odict:
            for _type2 in odict[_type1]:
                attrs = dir(NonBondedInteractionHandle(_type1, _type2))
                for a in attrs:
                    attr_ref = getattr(
                        NonBondedInteractionHandle(_type1, _type2), a)
                    type_name_ref = getattr(attr_ref, "type_name", None)
                    if callable(type_name_ref) and type_name_ref() == odict[_type1][_type2]['type_name']:
                        # found nonbonded inter, e.g.
                        # LennardJonesInteraction(_type1, _type2)
                        inter_instance = attr_ref
                        break
                    else:
                        continue

                del odict[_type1][_type2]['type_name']
                inter_instance.set_params(**odict[_type1][_type2])


cdef class BondedInteraction(object):

    # This means, the instance does not yet represent a bond in the simulation
    _bond_id = -1

    def __init__(self, *args, **kwargs):
        """
           Either called with an interaction id, in which case, the interaction
           will represent the bonded interaction as it is defined in Espresso core
           Or called with keyword arguments describing a new interaction.
        """
        # Interaction id as argument
        if len(args) == 1 and isinstance(args[0], int):
            bond_id = args[0]
            # Check, if the bond in Espresso core is really defined as a FENE
            # bond
            if bonded_ia_params[bond_id].type != self.type_number():
                raise Exception(
                    "The bond with this id is not defined as a " + self.type_name() + " bond in the Espresso core.")

            self._bond_id = bond_id

            # Load the parameters currently set in the Espresso core
            self._params = self._get_params_from_es_core()
            self._bond_id = bond_id

        # Or have we been called with keyword args describing the interaction
        elif len(args) == 0:
            # Check if all required keys are given
            for k in self.required_keys():
                if k not in kwargs:
                    raise ValueError(
                        "At least the following keys have to be given as keyword arguments: " + self.required_keys().__str__())

            self.params = kwargs

            # Validation of parameters
            self.validate_params()

        else:
            raise Exception(
                "The constructor has to be called either with a bond id (as interger), or with a set of keyword arguments describing a new interaction")

    def is_valid(self):
        """Check, if the data stored in the instance still matches what is in Espresso"""
        # Check if the bond type in Espresso still matches the bond type saved
        # in this class
        if bonded_ia_params[self._bond_id].type != self.type_number():
            return False

        # check, if the bond parameters saved in the class still match those
        # saved in Espresso
        temp_params = self._get_params_from_es_core()
        if self._params != temp_params:
            return False

        # If we're still here, the instance is valid
        return True

    property params:
        def __get__(self):
            return self._params

        def __set__(self, p):
            # Check, if any key was passed, which is not known
            for k in p.keys():
                if k not in self.valid_keys():
                    raise ValueError(
                        "Only the following keys are supported: " + self.valid_keys().__str__)

            # Initialize default values
            self.set_default_params()
            # Put in values given by the user
            self._params.update(p)

    def validate_params(self):
        return True

    def __getattribute__(self, name):
        """Every time _set_params_in_es_core is called, the parameter dict is also updated."""
        attr = object.__getattribute__(self, name)
        if hasattr(attr, '__call__') and attr.__name__ == "_set_params_in_es_core":
            def sync_params(*args, **kwargs):
                result = attr(*args, **kwargs)
                self._params.update(self._get_params_from_es_core())
                return result
            return sync_params
        else:
            return attr

    def _get_params_from_es_core(self):
        raise Exception(
            "Subclasses of BondedInteraction must define the _get_params_from_es_core() method.")

    def _set_params_in_es_core(self):
        raise Exception(
            "Subclasses of BondedInteraction must define the _set_params_in_es_core() method.")

    def __str__(self):
        return self.__class__.__name__ + "(" + str(self._params) + ")"

    def set_default_params(self):
        raise Exception(
            "Subclasses of BondedInteraction must define the set_default_params() method.")

    def type_number(self):
        raise Exception(
            "Subclasses of BondedInteraction must define the type_number() method.")

    def type_name(self):
        raise Exception(
            "Subclasses of BondedInteraction must define the type_name() method.")

    def valid_keys(self):
        raise Exception(
            "Subclasses of BondedInteraction must define the valid_keys() method.")

    def required_keys(self):
        raise Exception(
            "Subclasses of BondedInteraction must define the required_keys() method.")

    def __repr__(self):
        if self._bond_id == -1:
            id_str = "inactive"
        else:
            id_str = str(self._bond_id)

        return self.__class__.__name__ + "(" + id_str + "): " + self._params.__str__()

    def __richcmp__(self, other, i):
        if i != 2:
            raise Exception("only == supported")
        if self.__class__ != other.__class__:
            return False
        if self._bond_id != other._bond_id:
            return False
        return self._params == other._params


class BondedInteractionNotDefined(object):

    def __init__(self, *args, **kwargs):
        raise Exception(
            self.__class__.__name__ + " not compiled into Espresso core")

    def type_number(self):
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

    def type_name(self):
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

    def valid_keys(self):
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

    def required_keys(self):
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

    def set_default_params(self):
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

    def _get_params_from_es_core(self):
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

    def _set_params_in_es_core(self):
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)


class FeneBond(BondedInteraction):

    def __init__(self, *args, **kwargs):
        """
        FeneBond initialiser. Used to instatiate a FeneBond identifier
        with a given set of parameters.

        Parameters
        ----------
        k : float
            Specifies the magnitude of the bond interaction.
        d_r_max : float
                  Specifies the maximum stretch and compression length of the
                  bond.
        r_0 : float, optional
              Specifies the equilibrium length of the bond.
        """
        super(FeneBond, self).__init__(*args, **kwargs)

    def type_number(self):
        return BONDED_IA_FENE

    def type_name(self):
        return "FENE"

    def valid_keys(self):
        return "k", "d_r_max", "r_0"

    def required_keys(self):
        return "k", "d_r_max"

    def set_default_params(self):
        self._params = {"r_0": 0.}
        # Everything else has to be supplied by the user, anyway

    def _get_params_from_es_core(self):
        return \
            {"k": bonded_ia_params[self._bond_id].p.fene.k,
             "d_r_max": bonded_ia_params[self._bond_id].p.fene.drmax,
             "r_0": bonded_ia_params[self._bond_id].p.fene.r0}

    def _set_params_in_es_core(self):
        fene_set_params(
            self._bond_id, self._params["k"], self._params["d_r_max"], self._params["r_0"])


class HarmonicBond(BondedInteraction):

    def __init__(self, *args, **kwargs):
        """
        HarmonicBond initialiser. Used to instatiate a HarmonicBond identifier
        with a given set of parameters.

        Parameters
        ----------
        k : float
            Specifies the magnitude of the bond interaction.
        r_0 : float
              Specifies the equilibrium length of the bond.
        r_cut : float, optional
                Specifies maximum distance beyond which the bond is considered
                broken.
        """
        super(HarmonicBond, self).__init__(*args, **kwargs)

    def type_number(self):
        return BONDED_IA_HARMONIC

    def type_name(self):
        return "HARMONIC"

    def valid_keys(self):
        return "k", "r_0", "r_cut"

    def required_keys(self):
        return "k", "r_0"

    def set_default_params(self):
        self._params = {"k'": 0., "r_0": 0., "r_cut": 0.}

    def _get_params_from_es_core(self):
        return \
            {"k": bonded_ia_params[self._bond_id].p.harmonic.k,
             "r_0": bonded_ia_params[self._bond_id].p.harmonic.r,
             "r_cut": bonded_ia_params[self._bond_id].p.harmonic.r_cut}

    def _set_params_in_es_core(self):
        harmonic_set_params(
            self._bond_id, self._params["k"], self._params["r_0"], self._params["r_cut"])


IF ROTATION:
    class HarmonicDumbbellBond(BondedInteraction):

        def __init__(self, *args, **kwargs):
            """
            HarmonicDumbbellBond initialiser. Used to instatiate a
            HarmonicDumbbellBond identifier with a given set of parameters.

            Parameters
            ----------
            k1 : float
                Specifies the magnitude of the bond interaction.
            k2 : float
                Specifies the magnitude of the angular interaction.
            r_0 : float
                  Specifies the equilibrium length of the bond.
            r_cut : float, optional
                    Specifies maximum distance beyond which the bond is considered
                    broken.
            """
            super(HarmonicDumbbellBond, self).__init__(*args, **kwargs)

        def type_number(self):
            return BONDED_IA_HARMONIC_DUMBBELL

        def type_name(self):
            return "HARMONIC_DUMBBELL"

        def valid_keys(self):
            return "k1", "k2", "r_0", "r_cut"

        def required_keys(self):
            return "k1", "k2", "r_0"

        def set_default_params(self):
            self._params = {"r_cut": 0.}

        def _get_params_from_es_core(self):
            return \
                {"k1": bonded_ia_params[self._bond_id].p.harmonic_dumbbell.k1,
                 "k2": bonded_ia_params[self._bond_id].p.harmonic_dumbbell.k2,
                 "r_0": bonded_ia_params[self._bond_id].p.harmonic_dumbbell.r,
                 "r_cut": bonded_ia_params[self._bond_id].p.harmonic_dumbbell.r_cut}

        def _set_params_in_es_core(self):
            harmonic_dumbbell_set_params(
                self._bond_id, self._params["k1"], self._params["k2"],
                self._params["r_0"], self._params["r_cut"])

IF ROTATION != 1:
    class HarmonicDumbbellBond(BondedInteraction):

        def __init__(self, *args, **kwargs):
            """
            HarmonicDumbbellBond initialiser. Used to instatiate a
            HarmonicDumbbellBond identifier with a given set of parameters.

            Parameters
            ----------
            k1 : float
                Specifies the magnitude of the bond interaction.
            k2 : float
                Specifies the magnitude of the angular interaction.
            r_0 : float
                  Specifies the equilibrium length of the bond.
            r_cut : float, optional
                    Specifies maximum distance beyond which the bond is considered
                    broken.
            """
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def type_number(self):
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def type_name(self):
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def valid_keys(self):
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def required_keys(self):
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def set_default_params(self):
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def _get_params_from_es_core(self):
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def _set_params_in_es_core(self):
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")


IF BOND_CONSTRAINT == 1:
    class RigidBond(BondedInteraction):

        def __init__(self, *args, **kwargs):
            """
            RigidBond initialiser. Used to instantiate a RigidBond identifier
            with a given set of parameters.

            Parameters
            ----------
            r : float
                Specifies the length of the rigid bond.
            ptol : float, optional
                   Specifies the tolerance for positional deviations.
            vtop : float, optional
                   Specifies the tolerance for velocity deviations.
            """
            super(RigidBond, self).__init__(*args, **kwargs)

        def type_number(self):
            return BONDED_IA_RIGID_BOND

        def type_name(self):
            return "RIGID"

        def valid_keys(self):
            return "r", "ptol", "vtol"

        def required_keys(self):
            return "r"

        def set_default_params(self):
            # TODO rationality of Default Parameters has to be checked
            self._params = {"r": 0.,
                            "ptol": 0.001,
                            "vtol": 0.001}

        def _get_params_from_es_core(self):
            return {"r": bonded_ia_params[self._bond_id].p.rigid_bond.d2**0.5, "ptol": bonded_ia_params[self._bond_id].p.rigid_bond.p_tol, "vtol": bonded_ia_params[self._bond_id].p.rigid_bond.v_tol}

        def _set_params_in_es_core(self):
            rigid_bond_set_params(
                self._bond_id, self._params["r"], self._params["ptol"], self._params["vtol"])
ELSE:
    class RigidBond(BondedInteractionNotDefined):
        name = "RIGID"


class Dihedral(BondedInteraction):

    def type_number(self):
        return BONDED_IA_DIHEDRAL

    def type_name(self):
        return "DIHEDRAL"

    def valid_keys(self):
        return "mult", "bend", "phase"

    def required_keys(self):
        return "mult", "bend", "phase"

    def set_default_params(self):
        self._params = {"mult'": 1., "bend": 0., "phase": 0.}

    def _get_params_from_es_core(self):
        return \
            {"mult": bonded_ia_params[self._bond_id].p.dihedral.mult,
             "bend": bonded_ia_params[self._bond_id].p.dihedral.bend,
             "phase": bonded_ia_params[self._bond_id].p.dihedral.phase}

    def _set_params_in_es_core(self):
        dihedral_set_params(
            self._bond_id, self._params["mult"], self._params["bend"], self._params["phase"])


IF TABULATED == 1:
    class Tabulated(BondedInteraction):

        def __init__(self, *args, **kwargs):
            """
            RigidBond initialiser. Used to instantiate a RigidBond identifier
            with a given set of parameters.

            Parameters
            ----------
            type : str
                   Specifies the type of bonded interaction. Possible inputs:
                   'distance', 'angle' and 'dihedral'.
            filename : str
                       Filename of the tabular.
            """
            super(Tabulated, self).__init__(*args, **kwargs)

        def type_number(self):
            return BONDED_IA_TABULATED

        def type_name(self):
            return "TABULATED"

        def valid_keys(self):
            return "type", "filename", "npoints", "minval", "maxval", "invstepsize"

        def required_keys(self):
            return "type", "filename"

        def set_default_params(self):
            self._params = {"type": "bond", "filename": ""}

        def _get_params_from_es_core(self):
            make_bond_type_exist(self._bond_id)
            res = \
                {"type": bonded_ia_params[self._bond_id].p.tab.type,
                 "filename":
                     utils.to_str(
                         bonded_ia_params[self._bond_id].p.tab.filename),
                 "npoints": bonded_ia_params[self._bond_id].p.tab.npoints,
                 "minval": bonded_ia_params[self._bond_id].p.tab.minval,
                 "maxval": bonded_ia_params[self._bond_id].p.tab.maxval,
                 "invstepsize": bonded_ia_params[self._bond_id].p.tab.invstepsize}
            if res["type"] == 1:
                res["type"] = "distance"
            if res["type"] == 2:
                res["type"] = "angle"
            if res["type"] == 3:
                res["type"] = "dihedral"
            return res

        def _set_params_in_es_core(self):
            if self._params["type"] == "distance":
                type_num = 1
            else:
                if self._params["type"] == "angle":
                    type_num = 2
                else:
                    if self._params["type"] == "dihedral":
                        type_num = 3
                    else:
                        raise ValueError(
                            "Tabulated type needs to be distance, angle, or diherdal")

            res = tabulated_bonded_set_params(
                self._bond_id, < TabulatedBondedInteraction > type_num, utils.to_char_pointer(self._params["filename"]))
            msg = ""
            if res == 1:
                msg = "unknon bond type"
            if res == 3:
                msg = "cannot open file"
            if res == 4:
                msg = "file too short"
            if msg == 5:
                msg = "file broken"
            if msg == 6:
                msg = "parameter out of bound"
            if res:
                raise Exception("Could not setup tabulated bond. " + msg)
            # Retrieve some params, Es calculates.
            self._params = self._get_params_from_es_core()

    cdef class TabulatedNonBonded(NonBondedInteraction):

        cdef int state

        def __init__(self, *args, **kwargs):
            self.state = -1
            super(TabulatedNonBonded, self).__init__(*args, **kwargs)

        def type_number(self):
            return "TABULATED_NONBONDED"

        def type_name(self):
            return "TABULATED"

        def valid_keys(self):
            return "filename"

        def required_keys(self):
            return ["filename", ]

        def set_default_params(self):
            self._params = {"filename": ""}

        def _get_params_from_es_core(self):
            cdef ia_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "filename": utils.to_str(ia_params.TAB_filename)}

        def _set_params_in_es_core(self):
            self.state = tabulated_set_params(self._part_types[0], self._part_types[
                                              1], utils.to_char_pointer(self._params["filename"]))

        def is_active(self):
            if self.state == 0:
                return True

IF TABULATED != 1:
    class Tabulated(BondedInteraction):

        def type_number(self):
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def type_name(self):
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def valid_keys(self):
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def required_keys(self):
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def set_default_params(self):
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def _get_params_from_es_core(self):
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def _set_params_in_es_core(self):
            raise Exception("TABULATED has to be defined in myconfig.hpp.")


IF LENNARD_JONES == 1:
    class Subt_Lj(BondedInteraction):

        def type_number(self):
            return BONDED_IA_SUBT_LJ

        def type_name(self):
            return "SUBT_LJ"

        def valid_keys(self):
            return "r", "k"

        def required_keys(self):
            return "r", "k"

        def set_default_params(self):
            self._params = {"k": 0, "r": 0}

        def _get_params_from_es_core(self):
            return \
                {"k": bonded_ia_params[self._bond_id].p.subt_lj.k,
                 "r": bonded_ia_params[self._bond_id].p.subt_lj.r}

        def _set_params_in_es_core(self):
            subt_lj_set_params(
                self._bond_id, self._params["k"], self._params["r"])

IF BOND_VIRTUAL == 1:
    class Virtual(BondedInteraction):

        def __init__(self, *args, **kwargs):
            """
            VirtualBond initialiser. Used to instantiate a VirtualBond identifier.
            """
            super(Virtual, self).__init__(*args, **kwargs)

        def type_number(self):
            return BONDED_IA_VIRTUAL_BOND

        def type_name(self):
            return "VIRTUAL"

        def valid_keys(self):
            return {}

        def required_keys(self):
            return []

        def set_default_params(self):
            self._params = {}

        def _get_params_from_es_core(self):
            return {}

        def _set_params_in_es_core(self):
            virtual_set_params(self._bond_id)

ELSE:
    class Virtual(BondedInteractionNotDefined):
        name = "BOND_VIRTUAL"

IF BOND_ENDANGLEDIST == 1:
    class Endangledist(BondedInteraction):

        def type_number(self):
            return BONDED_IA_ENDANGLEDIST

        def type_name(self):
            return "ENDANGLEDIST"

        def valid_keys(self):
            return "bend", "phi0", "distmin", "distmax"

        def required_keys(self):
            return "bend", "phi0", "distmin", "distmax"

        def set_default_params(self):
            self._params = {"bend": 0, "phi0": 0, "distmin": 0, "distmax": 1}

        def _get_params_from_es_core(self):
            return \
                {"bend": bonded_ia_params[self._bond_id].p.endangledist.bend,
                 "phi0": bonded_ia_params[self._bond_id].p.endangledist.phi0,
                 "distmin":
                     bonded_ia_params[self._bond_id].p.endangledist.distmin,
                 "distmax": bonded_ia_params[self._bond_id].p.endangledist.distmax}

        def _set_params_in_es_core(self):
            endangledist_set_params(
                self._bond_id, self._params["bend"], self._params[
                    "phi0"], self._params["distmin"],
                                    self._params["distmax"])

ELSE:
    class Endangledist(BondedInteractionNotDefined):
        name = "BOND_ENDANGLEDIST"

IF OVERLAPPED == 1:
    class Overlapped(BondedInteraction):

        def type_number(self):
            return BONDED_IA_OVERLAPPED

        def type_name(self):
            return "OVERLAPPED"

        def valid_keys(self):
            return "overlap_type", "filename"

        def required_keys(self):
            return "overlap_type", "filename"

        def set_default_params(self):
            self._params = {"overlap_type": 0, "filename": ""}

        def _get_params_from_es_core(self):
            make_bond_type_exist(self._bond_id)
            return \
                {"bend": bonded_ia_params[self._bond_id].p.overlap.type,
                 "phi0": utils.to_str(bonded_ia_params[self._bond_id].p.overlap.filename)}

        def _set_params_in_es_core(self):
            overlapped_bonded_set_params(
                self._bond_id, self._params["overlap_type"], utils.to_char_pointer(self._params["filename"]))

ELSE:
    class Overlapped(BondedInteractionNotDefined):
        name = "OVERLAPPED"

IF BOND_ANGLE == 1:
    class Angle_Harmonic(BondedInteraction):

        def type_number(self):
            return BONDED_IA_ANGLE_HARMONIC

        def type_name(self):
            return "ANGLE_HARMONIC"

        def valid_keys(self):
            return "bend", "phi0"

        def required_keys(self):
            return "bend", "phi0"

        def set_default_params(self):
            self._params = {"bend": 0, "phi0": 0}

        def _get_params_from_es_core(self):
            return \
                {"bend": bonded_ia_params[self._bond_id].p.angle_harmonic.bend,
                 "phi0": bonded_ia_params[self._bond_id].p.angle_harmonic.phi0}

        def _set_params_in_es_core(self):
            angle_harmonic_set_params(
                self._bond_id, self._params["bend"], self._params["phi0"])
ELSE:
    class Angle_Harmonic(BondedInteractionNotDefined):
        name = "BOND_ANGLE"

IF BOND_ANGLE == 1:
    class Angle_Cosine(BondedInteraction):

        def type_number(self):
            return BONDED_IA_ANGLE_COSINE

        def type_name(self):
            return "ANGLE_COSINE"

        def valid_keys(self):
            return "bend", "phi0"

        def required_keys(self):
            return "bend", "phi0"

        def set_default_params(self):
            self._params = {"bend": 0, "phi0": 0}

        def _get_params_from_es_core(self):
            return \
                {"bend": bonded_ia_params[self._bond_id].p.angle_cosine.bend,
                 "phi0": bonded_ia_params[self._bond_id].p.angle_cosine.phi0}

        def _set_params_in_es_core(self):
            angle_cosine_set_params(
                self._bond_id, self._params["bend"], self._params["phi0"])
ELSE:
    class Angle_Cosine(BondedInteractionNotDefined):
        name = "BOND_ANGLE"

IF BOND_ANGLE == 1:
    class Angle_Cossquare(BondedInteraction):

        def type_number(self):
            return BONDED_IA_ANGLE_COSSQUARE

        def type_name(self):
            return "ANGLE_COSSQUARE"

        def valid_keys(self):
            return "bend", "phi0"

        def required_keys(self):
            return "bend", "phi0"

        def set_default_params(self):
            self._params = {"bend": 0, "phi0": 0}

        def _get_params_from_es_core(self):
            return \
                {"bend": bonded_ia_params[self._bond_id].p.angle_cossquare.bend,
                 "phi0": bonded_ia_params[self._bond_id].p.angle_cossquare.phi0}

        def _set_params_in_es_core(self):
            angle_cossquare_set_params(
                self._bond_id, self._params["bend"], self._params["phi0"])
ELSE:
    class Angle_Cossquare(BondedInteractionNotDefined):
        name = "BOND_ANGLE"


class Oif_Global_Forces(BondedInteraction):

    def type_number(self):
        return BONDED_IA_OIF_GLOBAL_FORCES

    def type_name(self):
        return "OIF_GLOBAL_FORCES"

    def valid_keys(self):
        return "A0_g", "ka_g", "V0", "kv"

    def required_keys(self):
        return "A0_g", "ka_g", "V0", "kv"

    def set_default_params(self):
        self._params = {"A0_g": 1., "ka_g": 0., "V0": 1., "kv": 0.}

    def _get_params_from_es_core(self):
        return \
            {"A0_g": bonded_ia_params[self._bond_id].p.oif_global_forces.A0_g,
             "ka_g": bonded_ia_params[self._bond_id].p.oif_global_forces.ka_g,
             "V0": bonded_ia_params[self._bond_id].p.oif_global_forces.V0,
             "kv": bonded_ia_params[self._bond_id].p.oif_global_forces.kv}

    def _set_params_in_es_core(self):
        oif_global_forces_set_params(
            self._bond_id, self._params["A0_g"], self._params["ka_g"], self._params["V0"], self._params["kv"])


class Oif_Local_Forces(BondedInteraction):

    def type_number(self):
        return BONDED_IA_OIF_LOCAL_FORCES

    def type_name(self):
        return "OIF_LOCAL_FORCES"

    def valid_keys(self):
        return "r0", "ks", "kslin", "phi0", "kb", "A01", "A02", "kal"

    def required_keys(self):
        return "r0", "ks", "kslin", "phi0", "kb", "A01", "A02", "kal"

    def set_default_params(self):
        self._params = {"r0": 1., "ks": 0., "kslin": 0.,
                        "phi0": 0., "kb": 0., "A01": 0., "A02": 0., "kal": 0.}

    def _get_params_from_es_core(self):
        return \
            {"r0": bonded_ia_params[self._bond_id].p.oif_local_forces.r0,
             "ks": bonded_ia_params[self._bond_id].p.oif_local_forces.ks,
             "kslin": bonded_ia_params[self._bond_id].p.oif_local_forces.kslin,
             "phi0": bonded_ia_params[self._bond_id].p.oif_local_forces.phi0,
             "kb": bonded_ia_params[self._bond_id].p.oif_local_forces.kb,
             "A01": bonded_ia_params[self._bond_id].p.oif_local_forces.A01,
             "A02": bonded_ia_params[self._bond_id].p.oif_local_forces.A02,
             "kal": bonded_ia_params[self._bond_id].p.oif_local_forces.kal}

    def _set_params_in_es_core(self):
        oif_local_forces_set_params(
            self._bond_id, self._params["r0"], self._params["ks"], self._params["kslin"], self._params["phi0"], self._params["kb"], self._params["A01"], self._params["A02"], self._params["kal"])


bonded_interaction_classes = {
    int(BONDED_IA_FENE): FeneBond,
    int(BONDED_IA_HARMONIC): HarmonicBond,
    int(BONDED_IA_HARMONIC_DUMBBELL): HarmonicDumbbellBond,
    int(BONDED_IA_RIGID_BOND): RigidBond,
    int(BONDED_IA_DIHEDRAL): Dihedral,
    int(BONDED_IA_TABULATED): Tabulated,
    int(BONDED_IA_VIRTUAL_BOND): Virtual,
    int(BONDED_IA_ENDANGLEDIST): Endangledist,
    int(BONDED_IA_OVERLAPPED): Overlapped,
    int(BONDED_IA_ANGLE_HARMONIC): Angle_Harmonic,
    int(BONDED_IA_ANGLE_COSINE): Angle_Cosine,
    int(BONDED_IA_ANGLE_COSSQUARE): Angle_Cossquare,
    int(BONDED_IA_OIF_GLOBAL_FORCES): Oif_Global_Forces,
    int(BONDED_IA_OIF_LOCAL_FORCES): Oif_Local_Forces,
}
IF LENNARD_JONES:
    bonded_interaction_classes[int(BONDED_IA_SUBT_LJ)] = Subt_Lj


class BondedInteractions(object):

    """Represents the bonded interactions. Individual interactions can be accessed using
    BondedInteractions[i], where i is the bond id. Will return a bonded interaction
    from bonded_interaction_classes"""

    def __getitem__(self, key):
        if not isinstance(key, int):
            raise ValueError(
                "Index to BondedInteractions[] has to be an integer referring to a bond id")

        # Find out the type of the interaction from Espresso
        if key >= n_bonded_ia:
            raise IndexError(
                "Index to BondedInteractions[] out of range")
        bond_type = bonded_ia_params[key].type

        # Check if the bonded interaction exists in Espresso core
        if bond_type == -1:
            raise ValueError(
                "The bonded interaction with the id " + str(key) + " is not yet defined.")

        # Find the appropriate class representing such a bond
        bond_class = bonded_interaction_classes[bond_type]

        # And return an instance of it, which refers to the bonded interaction
        # id in Espresso
        return bond_class(key)

    def __setitem__(self, key, value):
        # Validate arguments

        # type of key must be int
        if not isinstance(key, int):
            raise ValueError(
                "Index to BondedInteractions[] has to ba an integer referring to a bond id")

        # Value must be subclass off BondedInteraction
        if not isinstance(value, BondedInteraction):
            raise ValueError(
                "Only subclasses of BondedInteraction can be assigned.")

        # Save the bond id in the BondedInteraction instance
        value._bond_id = key

        # Set the parameters of the BondedInteraction instance in the Es core
        value._set_params_in_es_core()

    def __len__(self):
        return n_bonded_ia

    # Support iteration over active bonded interactions
    def __iter__(self):
        for i in range(n_bonded_ia):
            if bonded_ia_params[i].type != -1:
                yield self[i]

    def add(self, bonded_ia):
        """Add a bonded ia to the simulation>"""
        self[n_bonded_ia] = bonded_ia

    def __getstate__(self):
        params = {}
        for i, bonded_instance in enumerate(self):
            if hasattr(bonded_instance, 'params'):
                params[i] = bonded_instance.params
                params[i]['bond_type'] = bonded_instance.type_number()
            else:
                params[i] = None
        return params

    def __setstate__(self, params):
        for i in params:
            if params[i] is not None:
                bond_type = params[i]['bond_type']
                del params[i]['bond_type']
                self[i] = bonded_interaction_classes[bond_type](**params[i])
