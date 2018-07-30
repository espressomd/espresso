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

from libcpp.string cimport string
import collections

include "myconfig.pxi"
from . import utils
from espressomd.utils import is_valid_type
from globals cimport immersed_boundaries



# Non-bonded interactions

cdef class NonBondedInteraction(object):
    """
    Represents an instance of a non-bonded interaction, such as lennard jones
    Either called with two particle type id, in which case, the interaction
    will represent the bonded interaction as it is defined in Espresso core
    Or called with keyword arguments describing a new interaction.

    """

    cdef public object _part_types
    cdef public object _params

    # init dict to access all user defined nonbonded-inters via
    # user_interactions[type1][type2][parameter]
    cdef public object user_interactions

    def __init__(self, *args, **kwargs):
        if self.user_interactions == None:
            self.user_interactions = {}
        # Interaction id as argument
        if len(args) == 2 and is_valid_type(args[0], int) and is_valid_type(args[1], int):
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
            self.validate_params()
        else:
            raise Exception(
                "The constructor has to be called either with two particle type ids (as integer), or with a set of keyword arguments describing a new interaction")

    def is_valid(self):
        """Check, if the data stored in the instance still matches what is in Espresso.

        """
        temp_params = self._get_params_from_es_core()
        if self._params != temp_params:
            return False
        return True

    def get_params(self):
        """Get interaction parameters.

        """
        # If this instance refers to an actual interaction defined in the es core, load
        # current parameters from there
        if self._part_types[0] >= 0 and self._part_types[1] >= 0:
            self._params = self._get_params_from_es_core()
        return self._params

    def __str__(self):
        return self.__class__.__name__ + "(" + str(self.get_params()) + ")"

    def __getstate__(self):
        odict = collections.OrderedDict()
        odict['user_interactions'] = self.user_interactions
        odict['_part_types'] = self._part_types
        odict['params'] = self.get_params()
        return odict

    def __setstate__(self, state):
        self.user_interactions = state['user_interactions']
        self._part_types = state['_part_types']
        self._params = state['params']

    def set_params(self, **p):
        """Update the given parameters.

        """
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
        """Check that parameters are valid.

        """
        return True

    def __getattribute__(self, name):
        """Every time _set_params_in_es_core is called, the parameter dict is also updated.

        """
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
        """Virtual method.

        """
        raise Exception(
            "Subclasses of NonBondedInteraction must define the default_params() method.")

    def is_active(self):
        """Virtual method.

        """
        # If this instance refers to an actual interaction defined in the es core, load
        # current parameters from there
        if self._part_types[0] >= 0 and self._part_types[1] >= 0:
            self._params = self._get_params_from_es_core()
        raise Exception(
            "Subclasses of NonBondedInteraction must define the is_active() method.")

    def type_name(self):
        """Virtual method.

        """
        raise Exception(
            "Subclasses of NonBondedInteraction must define the type_name() method.")

    def valid_keys(self):
        """Virtual method.

        """
        raise Exception(
            "Subclasses of NonBondedInteraction must define the valid_keys() method.")

    def required_keys(self):
        """Virtual method.

        """
        raise Exception(
            "Subclasses of NonBondedInteraction must define the required_keys() method.")

# Lennard-Jones

IF LENNARD_JONES == 1:
    cdef class LennardJonesInteraction(NonBondedInteraction):

        def validate_params(self):
            """Check that parameters are valid.

            Raises
            ------
            ValueError
                If not true.
            """
            if self._params["epsilon"] < 0:
                raise ValueError("Lennard-Jones eps has to be >=0")
            if self._params["sigma"] < 0:
                raise ValueError("Lennard-Jones sigma has to be >=0")
            if self._params["cutoff"] < 0:
                raise ValueError("Lennard-Jones cutoff has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "epsilon": ia_params.LJ_eps,
                "sigma": ia_params.LJ_sig,
                "cutoff": ia_params.LJ_cut,
                "shift": ia_params.LJ_shift,
                "offset": ia_params.LJ_offset,
                "min": ia_params.LJ_min}

        def is_active(self):
            """Check if interaction is active.

            """
            return (self._params["epsilon"] > 0)

        def set_params(self, **kwargs):
            """ Set parameters for the Lennard-Jones interaction.

            Parameters
            ----------

            epsilon : :obj:`float`
                      The magnitude of the interaction.
            sigma : :obj:`float`
                    Determines the interaction length scale.
            cutoff : :obj:`float`
                     Cutoff distance of the interaction.
            shift : :obj:`float` or :obj:`str`
                    Constant shift of the potential. (4*epsilon*shift).
            offset : :obj:`float`, optional
                     Offset distance of the interaction.
            min : :obj:`float`, optional
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
                    self._params["min"]):
                raise Exception("Could not set Lennard Jones parameters")

        def default_params(self):
            """Python dictionary of default parameters.

            """
            return {
                "epsilon": 0.,
                "sigma": 0.,
                "cutoff": 0.,
                "shift": 0.,
                "offset": 0.,
                "min": 0.}

        def type_name(self):
            """Name of interaction type.

            """
            return "LennardJones"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "epsilon", "sigma", "cutoff", "shift", "offset", "min"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "epsilon", "sigma", "cutoff", "shift"

# Generic Lennard Jones

IF LENNARD_JONES_GENERIC == 1:

    cdef class GenericLennardJonesInteraction(NonBondedInteraction):

        def validate_params(self):
            """Check that parameters are valid.

            Raises
            ------
            ValueError
                If not true.
            """
            if self._params["epsilon"] < 0:
                raise ValueError("Generic Lennard-Jones eps has to be >=0")
            if self._params["sigma"] < 0:
                raise ValueError("Generic Lennard-Jones sigma has to be >=0")
            if self._params["cutoff"] < 0:
                raise ValueError("Generic Lennard-Jones cutoff has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
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
            """Check if interaction is active.

            """
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
                                    ):
                    raise Exception(
                        "Could not set Generic Lennard Jones parameters")

        def default_params(self):
            """Python dictionary of default parameters.

            """
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
            """Name of interaction type.

            """
            return "GenericLennardJones"

        def set_params(self, **kwargs):
            """
            Set parameters for the generic Lennard-Jones interaction.

            Parameters
            ----------
            epsilon : :obj:`float`
                      The magnitude of the interaction.
            sigma : :obj:`float`
                    Determines the interaction length scale.
            cutoff : :obj:`float`
                     Cutoff distance of the interaction.
            shift : :obj:`float`, string
                    Constant shift of the potential.
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
                    LJGEN_SOFTCORE parameter. Allows control over how smoothly
                    the potential drops to zero as lambda approaches zero.
            lam : :obj:`float`, optional
                     LJGEN_SOFTCORE parameter lambda. Tune the strength of the
                     interaction.

            """
            super(GenericLennardJonesInteraction, self).set_params(**kwargs)

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "epsilon", "sigma", "cutoff", "shift", "offset", "e1", "e2", "b1", "b2", "delta", "lam"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "epsilon", "sigma", "cutoff", "shift", "offset", "e1", "e2", "b1", "b2"

# Lennard-Jones cosine

IF LJCOS:
    cdef class LennardJonesCosInteraction(NonBondedInteraction):

        def validate_params(self):
            if self._params["epsilon"] < 0:
                raise ValueError("Lennard-Jones eps has to be >=0")
            if self._params["sigma"] < 0:
                raise ValueError("Lennard-Jones sigma has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "epsilon": ia_params.LJCOS_eps,
                "sigma": ia_params.LJCOS_sig,
                "cutoff": ia_params.LJCOS_cut,
                "offset": ia_params.LJCOS_offset,
            }

        def is_active(self):
            return(self._params["epsilon"] > 0)

        def set_params(self, **kwargs):
            """ Set parameters for the Lennard-Jones Cosine2 interaction.

            Parameters
            ----------

            epsilon : :obj:`float`
                      The magnitude of the interaction.
            sigma : :obj:`float`
                    Determines the interaction length scale.
            cutoff : :obj:`float`
                     Cutoff distance of the interaction.
            offset : :obj:`float`
                     Offset distance of the interaction.
            """
            super(LennardJonesCosInteraction, self).set_params(**kwargs)

        def _set_params_in_es_core(self):
            if ljcos_set_params(self._part_types[0], self._part_types[1],
                                self._params["epsilon"],
                                self._params["sigma"],
                                self._params["cutoff"],
                                self._params["offset"]):
                raise Exception("Could not set Lennard Jones Cosine parameters")

        def default_params(self):
            return {
                "epsilon": 0.,
                "sigma": 0.,
                "cutoff": 0.,
                "offset": 0.,
            }

        def type_name(self):
            """Name of interaction type.

            """
            return "LennardJonesCos"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "epsilon", "sigma", "cutoff", "offset"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "epsilon", "sigma", "cutoff"

# Lennard-Jones cosine^2

IF LJCOS2:
    cdef class LennardJonesCos2Interaction(NonBondedInteraction):

        def validate_params(self):
            if self._params["epsilon"] < 0:
                raise ValueError("Lennard-Jones eps has to be >=0")
            if self._params["sigma"] < 0:
                raise ValueError("Lennard-Jones sigma has to be >=0")
            if self._params["width"] < 0:
                raise ValueError("Parameter width has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return{
                "epsilon": ia_params.LJCOS2_eps,
                "sigma": ia_params.LJCOS2_sig,
                "offset": ia_params.LJCOS2_offset,
                "width": ia_params.LJCOS2_w}

        def is_active(self):
            return(self._params["epsilon"] > 0)

        def set_params(self, **kwargs):
            """ Set parameters for the Lennard-Jones Cosine2 interaction.

            Parameters
            ----------

            epsilon : :obj:`float`
                      The magnitude of the interaction.
            sigma : :obj:`float`
                    Determines the interaction length scale.
            offset : :obj:`float`
                     Offset distance of the interaction.
            width : :obj:`float`
                     Width of interaction. 
            """
            super(LennardJonesCos2Interaction, self).set_params(**kwargs)

        def _set_params_in_es_core(self):
            if ljcos2_set_params(self._part_types[0],
                                 self._part_types[1],
                                 self._params["epsilon"],
                                 self._params["sigma"],
                                 self._params["offset"],
                                 self._params["width"]):
                raise Exception("Could not set Lennard Jones Cosine2 parameters")

        def default_params(self):
            return {
                "epsilon": 0.,
                "sigma": 0.,
                "offset": 0.,
                "width": 0.}

        def type_name(self):
            """Name of interaction type.

            """
            return "LennardJonesCos2"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "epsilon", "sigma", "offset", "width"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "epsilon", "sigma", "width"

IF HAT == 1:
    cdef class HatInteraction(NonBondedInteraction):
        def validate_params(self):
            if self._params["F_max"] < 0:
                raise ValueError("Hat max force has to be >=0")
            if self._params["cutoff"] < 0:
                raise ValueError("Hat cutoff has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0], self._part_types[1])
            return {
                "F_max": ia_params.HAT_Fmax,
                "cutoff": ia_params.HAT_r,
            }

        def is_active(self):
            return (self._params["F_max"] > 0)

        def set_params(self, **kwargs):
            """ Set parameters for the Hat interaction.

            Parameters
            ----------
            F_max : :obj:`float`
                      The magnitude of the interaction.
            cutoff : :obj:`float`
                     Cutoff distance of the interaction.

            """
            super(HatInteraction, self).set_params(**kwargs)

        def _set_params_in_es_core(self):
            if hat_set_params(self._part_types[0], self._part_types[1],
                              self._params["F_max"],
                              self._params["cutoff"]):
                raise Exception("Could not set Hat parameters")

        def default_params(self):
            return {
                "F_max": 0.,
                "cutoff": 0.,
            }

        def type_name(self):
            return "Hat"

        def valid_keys(self):
            return "F_max", "cutoff"

        def required_keys(self):
            return "F_max", "cutoff"

# Gay-Berne

IF GAY_BERNE:
    cdef class GayBerneInteraction(NonBondedInteraction):

        def validate_params(self):
            """Check that parameters are valid.

            """
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
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
            """Check if interaction is active.

            """
            return (self._params["eps"] > 0)

        def set_params(self, **kwargs):
            """Set parameters for the Gay-Berne interaction.

            Parameters
            ----------
            eps : :obj:`float`
                  Potential well depth.
            sig : :obj:`float`
                  Interaction range.
            cut : :obj:`float`
                  Cutoff distance of the interaction.
            k1 : :obj:`float` or :obj:`string`
                  Molecular elongation.
            k2 : :obj:`float`, optional
                  Ratio of the potential well depths for the side-by-side
                  and end-to-end configurations.
            mu : :obj:`float`, optional
                  Adjustable exponent.
            nu  : :obj:`float`, optional
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
                raise Exception("Could not set Gay-Berne parameters")

        def default_params(self):
            """Python dictionary of default parameters.

            """
            return {
                "eps": 0.0,
                "sig": 0.0,
                "cut": 0.0,
                "k1": 0.0,
                "k2": 0.0,
                "mu": 0.0,
                "nu": 0.0}

        def type_name(self):
            """Name of interaction type.

            """
            return "GayBerne"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "eps", "sig", "cut", "k1", "k2", "mu", "nu"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "eps", "sig", "cut", "k1", "k2", "mu", "nu"

IF DPD:
    cdef class DPDInteraction(NonBondedInteraction):
        def validate_params(self):
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0], self._part_types[1])
            return {
                "weight_function": ia_params.dpd_wf,
                "gamma": ia_params.dpd_gamma,
                "r_cut": ia_params.dpd_r_cut,
                "trans_weight_function": ia_params.dpd_twf,
                "trans_gamma": ia_params.dpd_tgamma,
                "trans_r_cut": ia_params.dpd_tr_cut
            }

        def is_active(self):
            return (self._params["r_cut"] > 0) or (self._params["trans_r_cut"] > 0)

        def set_params(self, **kwargs):
            """ Set parameters for the DPD interaction.

            Parameters
            ----------
            weight_function : :obj:`float`
                The distance dependence of the parallel part,
                either 0 (constant) or 1 (linear)
            gamma : :obj:`float`
                Friction coefficient of the parallel part
            r_cut : :obj:`float`
                Cutoff of the parallel part
            trans_weight_function : :obj:`float`
                The distance dependence of the orthogonal part,
                either 0 (constant) or 1 (linear)
            trans_gamma : :obj:`float`
                Friction coefficient of the orthogonal part
            trans_r_cut : :obj:`float`
                Cutoff of the orthogonal part

            """
            super(DPDInteraction, self).set_params(**kwargs)

        def _set_params_in_es_core(self):
            if dpd_set_params(self._part_types[0], self._part_types[1],
                              self._params["gamma"],
                              self._params["r_cut"],
                              self._params["weight_function"],
                              self._params["trans_gamma"],
                              self._params["trans_r_cut"],
                              self._params["trans_weight_function"]):
                raise Exception("Could not set DPD parameters")

        def default_params(self):
            return {
                "weight_function": 0,
                "gamma": 0.0,
                "r_cut": -1.0,
                "trans_weight_function": 0,
                "trans_gamma": 0.0,
                "trans_r_cut": -1.0}

        def type_name(self):
            return "DPD"

        def valid_keys(self):
            return self.default_params().keys()

        def required_keys(self):
            return []

# Smooth-step

IF SMOOTH_STEP == 1:

    cdef class SmoothStepInteraction(NonBondedInteraction):

        def validate_params(self):
            """Check that parameters are valid.

            """
            if self._params["eps"] < 0:
                raise ValueError("Smooth-step eps has to be >=0")
            if self._params["offset"] < 0:
                raise ValueError("Smooth-step offset has to be >=0")
            if self._params["cutoff"] < 0:
                raise ValueError("Smooth-step cutoff has to be >=0")
            if self._params["cap"] < 0:
                raise ValueError("Smooth-step cap has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "d": ia_params.SmSt_d,
                "n": ia_params.SmSt_n,
                "eps": ia_params.SmSt_eps,
                "k0": ia_params.SmSt_k0,
                "sig": ia_params.SmSt_sig,
                "cutoff": ia_params.SmSt_cut
            }

        def is_active(self):
            """Check if interaction is active.

            """
            return ((self._params["eps"] > 0) and (self._params["sig"] > 0))

        def _set_params_in_es_core(self):
            if smooth_step_set_params(self._part_types[0], self._part_types[1],
                                      self._params["d"],
                                      self._params["n"],
                                      self._params["eps"],
                                      self._params["k0"],
                                      self._params["sig"],
                                      self._params["cutoff"]):
                raise Exception("Could not set smooth-step parameters")

        def default_params(self):
            """Python dictionary of default parameters.

            """
            return {
                "d": 0.,
                "n": 10,
                "eps": 0.,
                "k0": 0.,
                "sig": 0.,
                "cutoff": 0.}

        def type_name(self):
            """Name of interaction type.

            """
            return "SmoothStep"

        def set_params(self, **kwargs):
            """
            Set parameters for the smooth-step interaction.

            Parameters
            ----------
            d : :obj:`float`
                Short range repulsion parameter.
            n : :obj:`int`
                Exponent of short range repulsion.
            eps : :obj:`float`
                  The magnitude of the second (soft) repulsion.
            k0 : :obj:`float`
                 Exponential factor in second (soft) repulsion.
            sig : :obj:`float`
                  Length scale of second (soft) repulsion.
            cutoff : :obj:`float`
                Cutoff distance of the interaction.

            """
            super(SmoothStepInteraction, self).set_params(**kwargs)

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "d", "n", "eps", "k0", "sig", "cutoff"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "d", "eps", "cutoff"

# BMHTF

IF BMHTF_NACL == 1:

    cdef class BMHTFInteraction(NonBondedInteraction):

        def validate_params(self):
            """Check that parameters are valid.

            """
            if self._params["a"] < 0:
                raise ValueError("BMHTF a has to be >=0")
            if self._params["c"] < 0:
                raise ValueError("BMHTF c has to be >=0")
            if self._params["d"] < 0:
                raise ValueError("BMHTF d has to be >=0")
            if self._params["cutoff"] < 0:
                raise ValueError("BMHTF cutoff has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "a": ia_params.BMHTF_A,
                "b": ia_params.BMHTF_B,
                "c": ia_params.BMHTF_C,
                "d": ia_params.BMHTF_D,
                "sig": ia_params.BMHTF_sig,
                "cutoff": ia_params.BMHTF_cut,
            }

        def is_active(self):
            """Check if interaction is active.

            """
            return ((self._params["a"] > 0) and (self._params["c"] > 0) and (self._params["d"] > 0))

        def _set_params_in_es_core(self):
            if BMHTF_set_params(self._part_types[0], self._part_types[1],
                                self._params["a"],
                                self._params["b"],
                                self._params["c"],
                                self._params["d"],
                                self._params["sig"],
                                self._params["cutoff"]):
                raise Exception("Could not set BMHTF parameters")

        def default_params(self):
            """Python dictionary of default parameters.

            """
            return {
                "a": 0.,
                "b": 0.,
                "c": 0.,
                "d": 0.,
                "sig": 0.,
                "cutoff": 0.}

        def type_name(self):
            """Name of interaction type.

            """
            return "BMHTF"

        def set_params(self, **kwargs):
            """
            Set parameters for the BMHTF interaction.

            Parameters
            ----------
            a : :obj:`float`
                The magnitude of exponential part of the interaction.
            b : :obj:`float`
                Exponential factor of the interaction.
            c : :obj:`float`
                The magnitude of the term decaying with the sixth power of r.
            d : :obj:`float`
                The magnitude of the term decaying with the eighth power of r.
            sig : :obj:`float`
                Shift in the exponent.
            cutoff : :obj:`float`
                Cutoff distance of the interaction.

            """
            super(BMHTFInteraction, self).set_params(**kwargs)

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "a", "b", "c", "d", "sig", "cutoff"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "a", "b", "c", "d", "sig", "cutoff"

# Morse

IF MORSE == 1:

    cdef class MorseInteraction(NonBondedInteraction):

        def validate_params(self):
            """Check that parameters are valid.

            """
            if self._params["eps"] < 0:
                raise ValueError("Morse eps has to be >=0")
            if self._params["offset"] < 0:
                raise ValueError("Morse offset has to be >=0")
            if self._params["cutoff"] < 0:
                raise ValueError("Morse cutoff has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "eps": ia_params.MORSE_eps,
                "alpha": ia_params.MORSE_alpha,
                "rmin": ia_params.MORSE_rmin,
                "cutoff": ia_params.MORSE_cut
            }

        def is_active(self):
            """Check if interaction is active.

            """
            return (self._params["eps"] > 0)

        def _set_params_in_es_core(self):
            if morse_set_params(self._part_types[0], self._part_types[1],
                                self._params["eps"],
                                self._params["alpha"],
                                self._params["rmin"],
                                self._params["cutoff"]):
                raise Exception("Could not set Morse parameters")

        def default_params(self):
            """Python dictionary of default parameters.

            """
            return {
                "eps": 0.,
                "alpha": 0.,
                "rmin": 0.,
                "cutoff": 0.}

        def type_name(self):
            """Name of interaction type.

            """
            return "Morse"

        def set_params(self, **kwargs):
            """
            Set parameters for the Morse interaction.

            Parameters
            ----------
            eps : :obj:`float`
                The magnitude of the interaction.
            alpha : :obj:`float`
                Stiffness of the Morse interaction.
            rmin : :obj:`float`
                Distance of potential minimum
            cutoff : :obj:`float`
                Cutoff distance of the interaction.

            """
            super(MorseInteraction, self).set_params(**kwargs)

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "eps", "alpha", "rmin", "cutoff"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "eps", "alpha", "rmin"

# Buckingham

IF BUCKINGHAM == 1:

    cdef class BuckinghamInteraction(NonBondedInteraction):

        def validate_params(self):
            """Check that parameters are valid.

            """
            if self._params["a"] < 0:
                raise ValueError("Buckingham a has to be >=0")
            if self._params["b"] < 0:
                raise ValueError("Buckingham b has to be >=0")
            if self._params["c"] < 0:
                raise ValueError("Buckingham c has to be >=0")
            if self._params["d"] < 0:
                raise ValueError("Buckingham d has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "a": ia_params.BUCK_A,
                "b": ia_params.BUCK_B,
                "c": ia_params.BUCK_C,
                "d": ia_params.BUCK_D,
                "cutoff": ia_params.BUCK_cut,
                "discont": ia_params.BUCK_discont,
                "shift": ia_params.BUCK_shift
            }

        def is_active(self):
            """Check if interaction is active.

            """
            return ((self._params["a"] > 0) or (self._params["b"] > 0) or (self._params["d"] > 0) or (self._params["shift"] > 0))

        def _set_params_in_es_core(self):
            if buckingham_set_params(self._part_types[0], self._part_types[1],
                                     self._params["a"],
                                     self._params["b"],
                                     self._params["c"],
                                     self._params["d"],
                                     self._params["cutoff"],
                                     self._params["discont"],
                                     self._params["shift"]):
                raise Exception("Could not set Buckingham parameters")

        def default_params(self):
            """Python dictionary of default parameters.

            """
            return {
                "a": 0.,
                "b": 0.,
                "c": 0.,
                "d": 0.,
                "cutoff": 0.,
                "discont": 0.,
                "shift": 0.}

        def type_name(self):
            """Name of interaction type.

            """
            return "Buckingham"

        def set_params(self, **kwargs):
            """
            Set parameters for the Buckingham interaction.

            Parameters
            ----------
            a : :obj:`float`
                Magnitude of the exponential part of the interaction.
            b : :obj`float`
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
            super(BuckinghamInteraction, self).set_params(**kwargs)

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "a", "b", "c", "d", "discont", "cutoff", "shift"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "a", "c", "d", "discont", "cutoff"

# Soft-sphere

IF SOFT_SPHERE == 1:

    cdef class SoftSphereInteraction(NonBondedInteraction):

        def validate_params(self):
            """Check that parameters are valid.

            """
            if self._params["a"] < 0:
                raise ValueError("Soft-sphere a has to be >=0")
            if self._params["offset"] < 0:
                raise ValueError("Soft-sphere offset has to be >=0")
            if self._params["cutoff"] < 0:
                raise ValueError("Soft-sphere cutoff has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "a": ia_params.soft_a,
                "n": ia_params.soft_n,
                "cutoff": ia_params.soft_cut,
                "offset": ia_params.soft_offset
            }

        def is_active(self):
            """Check if interaction is active.

            """
            return (self._params["a"] > 0)

        def _set_params_in_es_core(self):
            if soft_sphere_set_params(self._part_types[0], self._part_types[1],
                                      self._params["a"],
                                      self._params["n"],
                                      self._params["cutoff"],
                                      self._params["offset"]):
                raise Exception("Could not set Soft-sphere parameters")

        def default_params(self):
            """Python dictionary of default parameters.

            """
            return {
                "a": 0.,
                "n": 0.,
                "cutoff": 0.,
                "offset": 0.}

        def type_name(self):
            """Name of interaction type.

            """
            return "SoftSphere"

        def set_params(self, **kwargs):
            """
            Set parameters for the Soft-sphere interaction.

            Parameters
            ----------
            a : :obj:`float`
                The magnitude of the interaction.
            n : :obj:`float`
                Exponent of the power law.
            cutoff : :obj:`float`
                     Cutoff distance of the interaction.
            offset : :obj:`float`, optional
                     Offset distance of the interaction.

            """
            super(SoftSphereInteraction, self).set_params(**kwargs)

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "a", "n", "cutoff", "offset"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "a", "n", "cutoff"

IF AFFINITY == 1:
    cdef class AffinityInteraction(NonBondedInteraction):

        def validate_params(self):
            if self._params["affinity_cut"] < 0:
                raise ValueError("Affinity cutoff has to be >=0")
            if self._params["affinity_type"] < 0:
                raise ValueError("Affinity type has to be >=0")
            if self._params["affinity_kappa"] < 0:
                raise ValueError("Affinity kappa has to be >=0")
            if self._params["affinity_r0"] < 0:
                raise ValueError("Affinity r0 has to be >=0")
            if self._params["affinity_Kon"] < 0:
                raise ValueError("Affinity Kon has to be >=0")
            if self._params["affinity_Koff"] < 0:
                raise ValueError("Affinity Koff has to be >=0")
            if self._params["affinity_maxBond"] < 0:
                raise ValueError("Affinity maxBond has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(self._part_types[0], self._part_types[1])
            return {
                "affinity_type": ia_params.affinity_type,
                "affinity_kappa": ia_params.affinity_kappa,
                "affinity_r0": ia_params.affinity_r0,
                "affinity_Kon": ia_params.affinity_Kon,
                "affinity_Koff": ia_params.affinity_Koff,
                "affinity_maxBond": ia_params.affinity_maxBond,
                "affinity_cut": ia_params.affinity_cut}

        def is_active(self):
            return (self._params["affinity_kappa"] > 0)


        def _set_params_in_es_core(self):
            if affinity_set_params(self._part_types[0], self._part_types[1],
                                        self._params["affinity_type"],
                                        self._params["affinity_kappa"],
                                        self._params["affinity_r0"],
                                        self._params["affinity_Kon"],
                                        self._params["affinity_Koff"],
                                        self._params["affinity_maxBond"],
                                        self._params["affinity_cut"]):
                raise Exception("Could not set Affinity parameters")

        def default_params(self):
            return {
                "affinity_type": 0,
                "affinity_kappa": 0.,
                "affinity_r0": 0.,
                "affinity_Kon": 0.,
                "affinity_Koff": 0.,
                "affinity_maxBond": 0.,
                "affinity_cut": 0.}

        def type_name(self):
            return "Affinity"

        def valid_keys(self):
            return "affinity_type", "affinity_kappa", "affinity_r0", "affinity_Kon", "affinity_Koff", "affinity_maxBond", "affinity_cut"

        def required_keys(self):
            return "affinity_type", "affinity_kappa", "affinity_r0", "affinity_Kon", "affinity_Koff", "affinity_maxBond", "affinity_cut"



IF MEMBRANE_COLLISION == 1:
    cdef class MembraneCollisionInteraction(NonBondedInteraction):

        def validate_params(self):
            if self._params["cutoff"] < 0:
                raise ValueError("Membrane Collision cutoff has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(self._part_types[0], self._part_types[1])
            return {
                "a": ia_params.membrane_a,
                "n": ia_params.membrane_n,
                "cutoff": ia_params.membrane_cut,
                "offset": ia_params.membrane_offset}

        def is_active(self):
            return (self._params["a"] > 0)

        def _set_params_in_es_core(self):
            if membrane_collision_set_params(self._part_types[0], self._part_types[1],
                                        self._params["a"],
                                        self._params["n"],
                                        self._params["cutoff"],
                                        self._params["offset"]):
                raise Exception("Could not set Membrane Collision parameters")

        def default_params(self):
            return {
                "a": 0.,
                "n": 1.,
                "cutoff": 0.,
                "offset": 0.}

        def type_name(self):
            return "MembraneCollision"

        def valid_keys(self):
            return "a", "n", "cutoff", "offset"

        def required_keys(self):
            return "a", "n", "cutoff"

# Gay-Berne


# Hertzian

IF HERTZIAN == 1:

    cdef class HertzianInteraction(NonBondedInteraction):

        def validate_params(self):
            """Check that parameters are valid.

            """
            if self._params["eps"] < 0:
                raise ValueError("Hertzian eps a has to be >=0")
            if self._params["sig"] < 0:
                raise ValueError("Hertzian sig has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "eps": ia_params.Hertzian_eps,
                "sig": ia_params.Hertzian_sig
            }

        def is_active(self):
            """Check if interaction is active.

            """
            return (self._params["eps"] > 0)

        def _set_params_in_es_core(self):
            if hertzian_set_params(self._part_types[0], self._part_types[1],
                                   self._params["eps"],
                                   self._params["sig"]):
                raise Exception("Could not set Hertzian parameters")

        def default_params(self):
            """Python dictionary of default parameters.

            """
            return {
                "eps": 0.,
                "sig": 1.}

        def type_name(self):
            """Name of interaction type.

            """
            return "Hertzian"

        def set_params(self, **kwargs):
            """
            Set parameters for the Hertzian interaction.

            Parameters
            ----------
            eps : :obj:`float`
                  The magnitude of the interaction.
            sig : :obj:`float`
                  Parameter sigma which determines the length over which the potential decays.

            """
            super(HertzianInteraction, self).set_params(**kwargs)

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "eps", "sig"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "eps", "sig"

# Gaussian

IF GAUSSIAN == 1:

    cdef class GaussianInteraction(NonBondedInteraction):

        def validate_params(self):
            """Check that parameters are valid.

            """
            if self._params["eps"] < 0:
                raise ValueError("Gaussian eps a has to be >=0")
            if self._params["sig"] < 0:
                raise ValueError("Gaussian sig has to be >=0")
            if self._params["offset"] < 0:
                raise ValueError("Gaussian offset has to be >=0")
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])
            return {
                "eps": ia_params.Gaussian_eps,
                "sig": ia_params.Gaussian_sig,
                "cutoff": ia_params.Gaussian_cut
            }

        def is_active(self):
            """Check if interaction is active.

            """
            return (self._params["eps"] > 0)

        def _set_params_in_es_core(self):
            if gaussian_set_params(self._part_types[0], self._part_types[1],
                                   self._params["eps"],
                                   self._params["sig"],
                                   self._params["cutoff"]):
                raise Exception(
                    "Could not set Gaussian interaction parameters")

        def default_params(self):
            """Python dictionary of default parameters.

            """
            return {
                "eps": 0.,
                "sig": 1.,
                "cutoff": 0.}

        def type_name(self):
            """Name of interaction type.

            """
            return "Gaussian"

        def set_params(self, **kwargs):
            """
            Set parameters for the Gaussian interaction.

            Parameters
            ----------
            eps : :obj:`float`
                  Overlap energy epsilon.
            sig : :obj:`float`
                  Variance sigma of the Gaussian interactin.
            cutoff : :obj:`float`
                     Cutoff distance of the interaction.

            """
            super(GaussianInteraction, self).set_params(**kwargs)

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "eps", "sig", "cutoff"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "eps", "sig", "cutoff"


class NonBondedInteractionHandle(object):
    """
    Provides access to all Non-bonded interactions between
    two particle types.

    """

    type1 = -1
    type2 = -1

    # Here, one line per non-bonded ia
    lennard_jones = None
    lennard_jones_cos = None
    lennard_jones_cos2 = None
    generic_lennard_jones = None
    smooth_step = None
    bmhtf = None
    morse = None
    buckingham = None
    soft_sphere = None
    hertzian = None
    gaussian = None
    tabulated = None
    soft_sphere = None
    membrane_collision = None
    gay_berne = None
    dpd = None
    hat = None
    thole = None

    def __init__(self, _type1, _type2):
        """Takes two particle types as argument"""
        if not (is_valid_type(_type1, int) and is_valid_type(_type2, int)):
            raise TypeError("The particle types have to be of type integer.")
        self.type1 = _type1
        self.type2 = _type2

        # Here, add one line for each nonbonded ia
        IF LENNARD_JONES:
            self.lennard_jones = LennardJonesInteraction(_type1, _type2)
        IF MEMBRANE_COLLISION:
            self.membrane_collision = MembraneCollisionInteraction(_type1, _type2)
        IF SOFT_SPHERE:
            self.soft_sphere = SoftSphereInteraction(_type1, _type2)
        IF AFFINITY:
            self.affinity = AffinityInteraction(_type1, _type2)
        IF LENNARD_JONES_GENERIC:
            self.generic_lennard_jones = GenericLennardJonesInteraction(
                _type1, _type2)
        IF LJCOS:
            self.lennard_jones_cos = LennardJonesCosInteraction(_type1, _type2)
        IF LJCOS2:
            self.lennard_jones_cos2 = LennardJonesCos2Interaction(
                _type1, _type2)
        IF SMOOTH_STEP:
            self.smooth_step = SmoothStepInteraction(_type1, _type2)
        IF BMHTF_NACL:
            self.bmhtf = BMHTFInteraction(_type1, _type2)
        IF MORSE:
            self.morse = MorseInteraction(_type1, _type2)
        IF BUCKINGHAM:
            self.buckingham = BuckinghamInteraction(_type1, _type2)
        IF SOFT_SPHERE:
            self.soft_sphere = SoftSphereInteraction(_type1, _type2)
        IF HERTZIAN:
            self.hertzian = HertzianInteraction(_type1, _type2)
        IF GAUSSIAN:
            self.gaussian = GaussianInteraction(_type1, _type2)
        IF TABULATED:
            self.tabulated = TabulatedNonBonded(_type1, _type2)
        IF GAY_BERNE:
            self.gay_berne = GayBerneInteraction(_type1, _type2)
        IF DPD:
            self.dpd = DPDInteraction(_type1, _type2)
        IF HAT:
            self.hat = HatInteraction(_type1, _type2)
        IF THOLE:
            self.thole = TholeInteraction(_type1, _type2)


cdef class NonBondedInteractions(object):
    """
    Access to non-bonded interaction parameters via [i,j], where i,j are particle
    types. Returns NonBondedInteractionHandle.
    Also: access to force capping.

    """
    def __getitem__(self, key):
        if not isinstance(key, tuple):
            raise ValueError(
                "NonBondedInteractions[] expects two particle types as indices.")
        if len(key) != 2 or (not is_valid_type(key[0], int)) or (not is_valid_type(key[1], int)):
            raise ValueError(
                "NonBondedInteractions[] expects two particle types as indices.")
        return NonBondedInteractionHandle(key[0], key[1])

    def __getstate__(self):
        cdef string core_state
        core_state = ia_params_get_state()
        return core_state

    def __setstate__(self, core_state):
        cdef string state = core_state
        ia_params_set_state(state)


cdef class BondedInteraction(object):
    """Base class for bonded interactions.

    """

    # This means, the instance does not yet represent a bond in the simulation
    _bond_id = -1

    def __init__(self, *args, **kwargs):
        """
        Either called with an interaction id, in which case, the interaction
        will represent the bonded interaction as it is defined in Espresso core
        Or called with keyword arguments describing a new interaction.

        """
        # Interaction id as argument
        if len(args) == 1 and is_valid_type(args[0], int):
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


    def __reduce__(self):
        return (self.__class__, (self._bond_id,))

    def is_valid(self):
        """Check, if the data stored in the instance still matches what is in Espresso.

        """
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
                        "Key '{}' invalid! Only the following keys are supported: {}"
                        .format(k, ", ".join(self.valid_keys())))

            # Initialize default values
            self.set_default_params()
            # Put in values given by the user
            self._params.update(p)

    def validate_params(self):
        """Check that parameters are valid.

        """
        return True

    def __getattribute__(self, name):
        """Every time _set_params_in_es_core is called, the parameter dict is also updated.

        """
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
        """Sets parameters that are not required to their default value.

        """
        raise Exception(
            "Subclasses of BondedInteraction must define the set_default_params() method.")

    def type_number(self):
        raise Exception(
            "Subclasses of BondedInteraction must define the type_number() method.")

    def type_name(self):
        """Name of interaction type.

        """
        raise Exception(
            "Subclasses of BondedInteraction must define the type_name() method.")

    def valid_keys(self):
        """All parameters that can be set.

        """
        raise Exception(
            "Subclasses of BondedInteraction must define the valid_keys() method.")

    def required_keys(self):
        """Parameters that have to be set.

        """
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
        """Name of interaction type.

        """
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

    def valid_keys(self):
        """All parameters that can be set.

        """
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

    def required_keys(self):
        """Parameters that have to be set.

        """
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

    def set_default_params(self):
        """Sets parameters that are not required to their default value.

        """
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

    def _get_params_from_es_core(self):
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)

    def _set_params_in_es_core(self):
        raise Exception(("%s has to be defined in myconfig.hpp.") % self.name)


class FeneBond(BondedInteraction):

    def __init__(self, *args, **kwargs):
        """
        FeneBond initializer. Used to instatiate a FeneBond identifier
        with a given set of parameters.

        Parameters
        ----------
        k : :obj:`float`
            Specifies the magnitude of the bond interaction.
        d_r_max : :obj:`float`
                  Specifies the maximum stretch and compression length of the
                  bond.
        r_0 : :obj:`float`, optional
              Specifies the equilibrium length of the bond.

        """
        super(FeneBond, self).__init__(*args, **kwargs)

    def type_number(self):
        return BONDED_IA_FENE

    def type_name(self):
        """Name of interaction type.

        """
        return "FENE"

    def valid_keys(self):
        """All parameters that can be set.

        """
        return "k", "d_r_max", "r_0"

    def required_keys(self):
        """Parameters that have to be set.

        """
        return "k", "d_r_max"

    def set_default_params(self):
        """Sets parameters that are not required to their default value.

        """
        """Sets parameters that are not required to their default value.

        """
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
        HarmonicBond initializer. Used to instatiate a HarmonicBond identifier
        with a given set of parameters.

        Parameters
        ----------
        k : :obj:`float`
            Specifies the magnitude of the bond interaction.
        r_0 : :obj:`float`
              Specifies the equilibrium length of the bond.
        r_cut : :obj:`float`, optional
                Specifies maximum distance beyond which the bond is considered
                broken.

        """
        super(HarmonicBond, self).__init__(*args, **kwargs)

    def type_number(self):
        return BONDED_IA_HARMONIC

    def type_name(self):
        """Name of interaction type.

        """
        return "HARMONIC"

    def valid_keys(self):
        """All parameters that can be set.

        """
        return "k", "r_0", "r_cut"

    def required_keys(self):
        """Parameters that have to be set.

        """
        return "k", "r_0"

    def set_default_params(self):
        """Sets parameters that are not required to their default value.

        """
        self._params = {"k": 0., "r_0": 0., "r_cut": 0.}

    def _get_params_from_es_core(self):
        return \
            {"k": bonded_ia_params[self._bond_id].p.harmonic.k,
             "r_0": bonded_ia_params[self._bond_id].p.harmonic.r,
             "r_cut": bonded_ia_params[self._bond_id].p.harmonic.r_cut}

    def _set_params_in_es_core(self):
        harmonic_set_params(
            self._bond_id, self._params["k"], self._params["r_0"], self._params["r_cut"])

if ELECTROSTATICS:
    class BondedCoulomb(BondedInteraction):

        def __init__(self, *args, **kwargs):
            """ 
            BondedCoulombBond initialiser. Used to instatiate a BondedCoulombBond identifier
            with a given set of parameters. 

            Parameters
            ----------

            prefactor : :obj:`float`
                        Sets the coulomb prefactor of the bonded coulomb interaction.
            """
            super(BondedCoulomb, self).__init__(*args, **kwargs)


        def type_number(self):
            return BONDED_IA_BONDED_COULOMB

        def type_name(self):
            return "BONDED_COULOMB"

        def valid_keys(self):
            return {"prefactor"}

        def required_keys(self):
            return {"prefactor"}

        def set_default_params(self):
            self._params = {"prefactor": 1.}

        def _get_params_from_es_core(self):
            return \
                {"prefactor": bonded_ia_params[self._bond_id].p.bonded_coulomb.prefactor}

        def _set_params_in_es_core(self):
            bonded_coulomb_set_params(
                self._bond_id,  self._params["prefactor"])

if P3M:
    class BondedCoulombP3MSRBond(BondedInteraction):

        def __init__(self, *args, **kwargs):
            """ 
            BondedCoulombP3MSRBond initialiser. Used to instatiate a BondedCoulombP3MSRBond identifier
            with a given set of parameters. Calculates ony the P3M shortrange part.

            Parameters
            ----------

            q1q2 : :obj:`float`
                   Sets the charge factor of the involved particle pair. Not the particle charges are used to allow e.g. only partial subtraction of the involved charges.
            """
            super(BondedCoulombP3MSRBond, self).__init__(*args, **kwargs)


        def type_number(self):
            return BONDED_IA_BONDED_COULOMB_P3M_SR

        def type_name(self):
            return "BONDED_COULOMB_P3M_SR"

        def valid_keys(self):
            return {"q1q2"}

        def required_keys(self):
            return {"q1q2"}

        def set_default_params(self):
            self._params = {"q1q2": 1.}

        def _get_params_from_es_core(self):
            return \
                {"q1q2": bonded_ia_params[self._bond_id].p.bonded_coulomb_p3m_sr.q1q2}

        def _set_params_in_es_core(self):
            bonded_coulomb_p3m_sr_set_params(
                self._bond_id,  self._params["q1q2"])
    
class ThermalizedBond(BondedInteraction):

    def __init__(self, *args, **kwargs):
        """ 
        ThermalizedBond initialiser. Used to instatiate a ThermalizedBond identifier
        with a given set of parameters. 

        Parameters
        ----------

        temp_com : :obj:`float`
                    Sets the temerature of the Langevin thermostat for the com of the particle pair.
        gamma_com: :obj:`float`
                    Sets the friction coefficient of the Langevin thermostat for the com of the particle pair.
        temp_distance: :obj:`float` 
                    Sets the temerature of the Langevin thermostat for the distance vector of the particle pair.
        gamma_distance: :obj:`float` 
                     Sets the friction coefficient of the Langevin thermostat for the distance vector of the particle pair.
        r_cut: :obj:`float`, optional
                Specifies maximum distance beyond which the bond is considered broken.
        """
        super(ThermalizedBond, self).__init__(*args, **kwargs)


    def type_number(self):
        return BONDED_IA_THERMALIZED_DIST

    def type_name(self):
        return "THERMALIZED_DIST"

    def valid_keys(self):
        return "temp_com", "gamma_com", "temp_distance", "gamma_distance", "r_cut"

    def required_keys(self):
        return "temp_com", "gamma_com", "temp_distance", "gamma_distance"

    def set_default_params(self):
        self._params = {"temp_com": 1., "gamma_com": 1., "temp_distance": 1., "gamma_distance": 1., "r_cut": 0.}

    def _get_params_from_es_core(self):
        return \
            {"temp_com": bonded_ia_params[self._bond_id].p.thermalized_bond.temp_com,
             "gamma_com": bonded_ia_params[self._bond_id].p.thermalized_bond.gamma_com,
             "temp_distance": bonded_ia_params[self._bond_id].p.thermalized_bond.temp_distance,
             "gamma_distance": bonded_ia_params[self._bond_id].p.thermalized_bond.gamma_distance,
             "r_cut": bonded_ia_params[self._bond_id].p.thermalized_bond.r_cut}

    def _set_params_in_es_core(self):
        thermalized_bond_set_params(
            self._bond_id,  self._params["temp_com"],  self._params["gamma_com"],  self._params["temp_distance"],  self._params["gamma_distance"], self._params["r_cut"])

IF THOLE:
    cdef class TholeInteraction(NonBondedInteraction):

        def validate_params(self):
            return True

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params
            ia_params = get_ia_param_safe(self._part_types[0], self._part_types[1])
            return {
                "scaling_coeff": ia_params.THOLE_scaling_coeff,
                "q1q2": ia_params.THOLE_q1q2
            }

        def is_active(self):
            return (self._params["scaling_coeff"] != 0)

        def set_params(self, **kwargs):
            """ Set parameters for the Thole interaction.

            Parameters
            ----------
            scaling_coeff : :obj:`float`
                            The facor used in the thole damping function between 
                            polarizable particles i and j. Usually caluclated by 
                            the polarizabilities alpha_i, alpha_j and damping 
                            parameters  a_i, a_j via
                            scaling_coeff = (a_i+a_j)/2 / ((alpha_i*alpha_j)^(1/2))^(1/3)
            q1q2: :obj:`float`
                  charge factor of the involved charges. Has to be set because 
                  it acts only on the portion of the drude core charge that is 
                  associated to the dipole of the atom. For charged, polarizable 
                  atoms that charge is not equal to the particle charge property.

            """
            super(TholeInteraction, self).set_params(**kwargs)

        def _set_params_in_es_core(self):
            # Handle the case of shift="auto"
            if thole_set_params(self._part_types[0], self._part_types[1],
                                        self._params["scaling_coeff"],
                                        self._params["q1q2"]):
                raise Exception("Could not set Thole parameters")

        def default_params(self):
            return {"scaling_coeff": 0., "q1q2": 0.}

        def type_name(self):
            return "Thole"

        def valid_keys(self):
            return "scaling_coeff", "q1q2"

        def required_keys(self):
            return "scaling_coeff", "q1q2"

IF ROTATION:
    class HarmonicDumbbellBond(BondedInteraction):

        def __init__(self, *args, **kwargs):
            """
            HarmonicDumbbellBond initializer. Used to instatiate a
            HarmonicDumbbellBond identifier with a given set of parameters.

            Parameters
            ----------
            k1 : :obj:`float`
                Specifies the magnitude of the bond interaction.
            k2 : :obj:`float`
                Specifies the magnitude of the angular interaction.
            r_0 : :obj:`float`
                  Specifies the equilibrium length of the bond.
            r_cut : :obj:`float`, optional
                    Specifies maximum distance beyond which the bond is considered
                    broken.

            """
            super(HarmonicDumbbellBond, self).__init__(*args, **kwargs)

        def type_number(self):
            return BONDED_IA_HARMONIC_DUMBBELL

        def type_name(self):
            """Name of interaction type.

            """
            return "HARMONIC_DUMBBELL"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "k1", "k2", "r_0", "r_cut"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "k1", "k2", "r_0"

        def set_default_params(self):
            """Sets parameters that are not required to their default value.

            """
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
            HarmonicDumbbellBond initializer. Used to instatiate a
            HarmonicDumbbellBond identifier with a given set of parameters.

            Parameters
            ----------
            k1 : :obj:`float`
                Specifies the magnitude of the bond interaction.
            k2 : :obj:`float`
                Specifies the magnitude of the angular interaction.
            r_0 : :obj:`float`
                  Specifies the equilibrium length of the bond.
            r_cut : :obj:`float`, optional
                    Specifies maximum distance beyond which the bond is considered
                    broken.

            """
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def type_number(self):
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def type_name(self):
            """Name of interaction type.

            """
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def valid_keys(self):
            """All parameters that can be set.

            """
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def required_keys(self):
            """Parameters that have to be set.

            """
            raise Exception(
                "HarmonicDumbbellBond: ROTATION has to be defined in myconfig.hpp.")

        def set_default_params(self):
            """Sets parameters that are not required to their default value.

            """
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
            RigidBond initializer. Used to instantiate a RigidBond identifier
            with a given set of parameters.

            Parameters
            ----------
            r : :obj:`float`
                Specifies the length of the rigid bond.
            ptol : :obj:`float`, optional
                   Specifies the tolerance for positional deviations.
            vtop : :obj:`float`, optional
                   Specifies the tolerance for velocity deviations.

            """
            super(RigidBond, self).__init__(*args, **kwargs)

        def type_number(self):
            return BONDED_IA_RIGID_BOND

        def type_name(self):
            """Name of interaction type.

            """
            return "RIGID"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "r", "ptol", "vtol"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "r"

        def set_default_params(self):
            """Sets parameters that are not required to their default value.

            """
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
        """Name of interaction type.

        """
        return "DIHEDRAL"

    def valid_keys(self):
        """All parameters that can be set.

        """
        return "mult", "bend", "phase"

    def required_keys(self):
        """Parameters that have to be set.

        """
        return "mult", "bend", "phase"

    def set_default_params(self):
        """Sets parameters that are not required to their default value.

        """
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
            Tabulated bond initializer. Used to instantiate a Tabulated bond identifier
            with a given set of parameters.

            Parameters
            ----------

            type : :obj:`string`,
                The type of bond, one of 'distance', 'angle' or
                'dihedral'.
            min : :obj:`float`,
                The minimal interaction distance. Has to be 0 if
                type is 'angle' or 'dihedral'
            max : :obj:`float`,
                The maximal interaction distance. Has to be pi if
                type is 'angle' or 2pi if 'dihedral'
            energy: array_like :obj:`float`
                The energy table.
            force: array_like :obj:`float`
                The force table.

            """
            super(Tabulated, self).__init__(*args, **kwargs)

        def type_number(self):
            return BONDED_IA_TABULATED

        def type_name(self):
            """Name of interaction type.

            """
            return "TABULATED"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "type", "min", "max", "energy", "force"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "type", "min", "max", "energy", "force"

        def set_default_params(self):
            """Sets parameters that are not required to their default value.

            """
            self._params = {'min': -1., 'max': -1., 'energy': [], 'force': []}

        def _get_params_from_es_core(self):
            make_bond_type_exist(self._bond_id)
            res = \
                {"type": bonded_ia_params[self._bond_id].p.tab.type,
                 "min": bonded_ia_params[self._bond_id].p.tab.pot.minval,
                 "max": bonded_ia_params[self._bond_id].p.tab.pot.maxval,
                 "energy": bonded_ia_params[self._bond_id].p.tab.pot.energy_tab,
                 "force": bonded_ia_params[self._bond_id].p.tab.pot.force_tab
                }
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
            elif self._params["type"] == "angle":
                type_num = 2
            elif self._params["type"] == "dihedral":
                type_num = 3
            else:
                raise ValueError(
                    "Tabulated type needs to be distance, angle, or diherdal")

            res = tabulated_bonded_set_params(
                self._bond_id, < TabulatedBondedInteraction > type_num,
                self._params["min"],
                self._params["max"],
                self._params["energy"],
                self._params["force"])

            if res == 1:
                raise Exception("Could not setup tabulated bond. Invalid bond type.")
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
            """Name of the potential.

            """
            return "TABULATED"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "min", "max", "energy", "force"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return ["min", "max", "energy", "force"]

        def set_params(self, **kwargs):
            """ Set parameters for the TabulatedNonBonded interaction.

            Parameters
            ----------

            min : :obj:`float`,
                  The minimal interaction distance.
            max : :obj:`float`,
                  The maximal interaction distance.
            energy: array_like :obj:`float`
                  The energy table.
            force: array_like :obj:`float`
                  The force table.

            """
            super(TabulatedNonBonded, self).set_params(**kwargs)

        def set_default_params(self):
            """Sets parameters that are not required to their default value.

            """
            self._params = {'min': -1., 'max': -1, 'energy': [], 'force': []}

        def _get_params_from_es_core(self):
            cdef IA_parameters * ia_params = get_ia_param_safe(
                self._part_types[0],
                self._part_types[1])

            return {'min': ia_params.TAB.minval,
                    'max': ia_params.TAB.maxval,
                    'energy': ia_params.TAB.energy_tab,
                    'force': ia_params.TAB.force_tab}

        def _set_params_in_es_core(self):
            self.state = tabulated_set_params(self._part_types[0],
                                              self._part_types[1],
                                              self._params["min"],
                                              self._params["max"],
                                              self._params["energy"],
                                              self._params["force"])

        def is_active(self):
            """Check if interaction is active.

            """
            if self.state == 0:
                return True

IF TABULATED != 1:
    class Tabulated(BondedInteraction):
        def type_number(self):
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def type_name(self):
            """Name of interaction type.

            """
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def valid_keys(self):
            """All parameters that can be set.

            """
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def required_keys(self):
            """Parameters that have to be set.

            """
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def set_default_params(self):
            """Sets parameters that are not required to their default value.

            """
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def _get_params_from_es_core(self):
            raise Exception("TABULATED has to be defined in myconfig.hpp.")

        def _set_params_in_es_core(self):
            raise Exception("TABULATED has to be defined in myconfig.hpp.")


IF LENNARD_JONES == 1:
    class SubtLJ(BondedInteraction):
        def __init__(self, *args, **kwargs):
            super(SubtLJ, self).__init__(*args, **kwargs)

        def type_number(self):
            return BONDED_IA_SUBT_LJ

        def type_name(self):
            """Name of interaction type.

            """
            return "SUBT_LJ"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return {}

        def required_keys(self):
            """Parameters that have to be set.

            """
            return {}

        def set_default_params(self):
            """Sets parameters that are not required to their default value.

            """
            self._params = {}

        def _get_params_from_es_core(self):
            return {}

        def _set_params_in_es_core(self):
            subt_lj_set_params(self._bond_id)


class Virtual(BondedInteraction):
    def __init__(self, *args, **kwargs):
        """
        VirtualBond initializer. Used to instantiate a VirtualBond identifier.
        """
        super(Virtual, self).__init__(*args, **kwargs)

    def type_number(self):
        return BONDED_IA_VIRTUAL_BOND

    def type_name(self):
        """Name of interaction type.

        """
        return "VIRTUAL"

    def valid_keys(self):
        """All parameters that can be set.

        """
        return {}

    def required_keys(self):
        """Parameters that have to be set.

        """
        return []

    def set_default_params(self):
        """Sets parameters that are not required to their default value.

        """
        self._params = {}

    def _get_params_from_es_core(self):
        return {}

    def _set_params_in_es_core(self):
        virtual_set_params(self._bond_id)

IF BOND_ANGLE == 1:
    class AngleHarmonic(BondedInteraction):
        """
        Bond angle dependent harmonic potential.

        See :ref:`Bond-angle interactions`

        """

        def type_number(self):
            return BONDED_IA_ANGLE_HARMONIC

        def type_name(self):
            """Name of interaction type.

            """
            return "ANGLE_HARMONIC"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "bend", "phi0"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "bend", "phi0"

        def set_default_params(self):
            """Sets parameters that are not required to their default value.

            """
            self._params = {"bend": 0, "phi0": 0}

        def _get_params_from_es_core(self):
            return \
                {"bend": bonded_ia_params[self._bond_id].p.angle_harmonic.bend,
                 "phi0": bonded_ia_params[self._bond_id].p.angle_harmonic.phi0}

        def _set_params_in_es_core(self):
            angle_harmonic_set_params(
                self._bond_id, self._params["bend"], self._params["phi0"])
ELSE:
    class AngleHarmonic(BondedInteractionNotDefined):
        name = "AngleHarmonic"

IF BOND_ANGLE == 1:
    class AngleCosine(BondedInteraction):
        """
        Bond angle dependent ine potential.

        See :ref:`Bond-angle interactions`

        """

        def type_number(self):
            return BONDED_IA_ANGLE_COSINE

        def type_name(self):
            """Name of interaction type.

            """
            return "ANGLE_COSINE"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "bend", "phi0"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "bend", "phi0"

        def set_default_params(self):
            """Sets parameters that are not required to their default value.

            """
            self._params = {"bend": 0, "phi0": 0}

        def _get_params_from_es_core(self):
            return \
                {"bend": bonded_ia_params[self._bond_id].p.angle_cosine.bend,
                 "phi0": bonded_ia_params[self._bond_id].p.angle_cosine.phi0}

        def _set_params_in_es_core(self):
            angle_cosine_set_params(
                self._bond_id, self._params["bend"], self._params["phi0"])
ELSE:
    class AngleCosine(BondedInteractionNotDefined):
        name = "AngelCosine"

IF BOND_ANGLE == 1:
    class AngleCossquare(BondedInteraction):
        """
        Bond angle dependent cos^2 potential.

        See :ref:`Bond-angle interactions`

        """

        def type_number(self):
            return BONDED_IA_ANGLE_COSSQUARE

        def type_name(self):
            """Name of interaction type.

            """
            return "ANGLE_COSSQUARE"

        def valid_keys(self):
            """All parameters that can be set.

            """
            return "bend", "phi0"

        def required_keys(self):
            """Parameters that have to be set.

            """
            return "bend", "phi0"

        def set_default_params(self):
            """Sets parameters that are not required to their default value.

            """
            self._params = {"bend": 0, "phi0": 0}

        def _get_params_from_es_core(self):
            return \
                {"bend": bonded_ia_params[self._bond_id].p.angle_cossquare.bend,
                 "phi0": bonded_ia_params[self._bond_id].p.angle_cossquare.phi0}

        def _set_params_in_es_core(self):
            angle_cossquare_set_params(
                self._bond_id, self._params["bend"], self._params["phi0"])
ELSE:
    class AngleCossquare(BondedInteractionNotDefined):
        name = "AngleCossquare"

# IBM triel
IF IMMERSED_BOUNDARY:
    class IBM_Triel(BondedInteraction):

        def __init__(self, *args, **kwargs):
            """
            IBM_Triel initializer. Used to instatiate an IBM_Triel identifier
            with a given set of parameters.

            Parameters
            ----------
            indX : :obj:`int`
                  Specify first, second and third bonding partner. Used for initializing reference state
            k1 : :obj:`float`
                  Specifies shear elasticity for Skalak and Neo-Hookean
            k2 : :obj:`float`
                  Specifies area resistance for Skalak
            maxDist : :obj:`float`
                  Gives an error if an edge becomes longer than maxDist
            elasticLaw : :obj:`string`
                  Specify NeoHookean or Skalak
    
            """
            super(IBM_Triel, self).__init__(*args, **kwargs)

        def type_number(self):
            return BONDED_IA_IBM_TRIEL

        def type_name(self):
            return "IBM_Triel"

        def valid_keys(self):
            return "ind1", "ind2", "ind3", "k1","k2", "maxDist", "elasticLaw"

        def required_keys(self):
            return "ind1", "ind2", "ind3","k1", "maxDist", "elasticLaw"

        def set_default_params(self):
            self._params = {"k2":0}

        def _get_params_from_es_core(self):
            return \
               {"maxDist": bonded_ia_params[self._bond_id].p.ibm_triel.maxDist,
                 "k1": bonded_ia_params[self._bond_id].p.ibm_triel.k1,
                 "k2": bonded_ia_params[self._bond_id].p.ibm_triel.k2,
                 "elasticLaw": <int>bonded_ia_params[self._bond_id].p.ibm_triel.elasticLaw}

        def _set_params_in_es_core(self):
            cdef tElasticLaw el
            if self._params["elasticLaw"] == "NeoHookean":
                el = NeoHookean
            if self._params["elasticLaw"] == "Skalak":
                el = Skalak
            IBM_Triel_SetParams(self._bond_id,self._params["ind1"],self._params["ind2"],self._params["ind3"],self._params["maxDist"], el, self._params["k1"], self._params["k2"])
ELSE:
    class IBM_Triel(BondedInteractionNotDefined):
        name = "IBM_TRIEL"

# IBM tribend
IF IMMERSED_BOUNDARY == 1:
    class IBM_Tribend(BondedInteraction):

        def __init__(self, *args, **kwargs):
            """
            IBM_Tribend initializer. Used to instatiate an IBM_Tribend identifier
            with a given set of parameters.

            Parameters
            ----------
            indX : :obj:`int`
                  Specify first, second, third and fourth bonding partner. Used for initializing reference state
            kb : :obj:`float`
                  Specifies bending modulus
            refShape : :obj:`string`
                  Flat or Initial

            """
            super(IBM_Tribend, self).__init__(*args, **kwargs)

        def type_number(self):
            return BONDED_IA_IBM_TRIBEND

        def type_name(self):
            return "IBM_Tribend"

        def valid_keys(self):
            return "ind1", "ind2", "ind3", "ind4", "kb", "refShape"

        def required_keys(self):
            return "ind1", "ind2", "ind3", "ind4", "kb"

        def set_default_params(self):
            self._params = {"flat":False}

        def _get_params_from_es_core(self):
            return \
               {"kb": bonded_ia_params[self._bond_id].p.ibm_tribend.kb,
                "theta0": bonded_ia_params[self._bond_id].p.ibm_tribend.theta0}

        def _set_params_in_es_core(self):
            if self._params["refShape"] == "Flat":
                flat = True
            if self._params["refShape"] == "Initial":
                flat = False
            IBM_Tribend_SetParams(self._bond_id,self._params["ind1"],self._params["ind2"],self._params["ind3"],self._params["ind4"], self._params["kb"], flat)
ELSE:
    class IBM_Tribend(BondedInteractionNotDefined):
        name = "IBM_TRIBEND"

# IBM VolCons
IF IMMERSED_BOUNDARY == 1:
    class IBM_VolCons(BondedInteraction):

        def __init__(self, *args, **kwargs):
            """
            IBM_VolCons initializer. Used to instatiate an IBM_VolCons identifier
            with a given set of parameters.

            Parameters
            ----------
            softID : :obj:`int`
                  Used to identify the object to which this bond belongs. Each object (cell) needs its own ID
            kappaV : :obj:`float`
                  Modulus for volume force

            """
            super(IBM_VolCons, self).__init__(*args, **kwargs)

        def type_number(self):
            return BONDED_IA_IBM_VOLUME_CONSERVATION

        def type_name(self):
            return "IBM_VolCons"

        def valid_keys(self):
            return "softID", "kappaV"

        def required_keys(self):
            return "softID", "kappaV"

        def set_default_params(self):
            self._params = {}

        def _get_params_from_es_core(self):
            return \
              {"softID": bonded_ia_params[self._bond_id].p.ibmVolConsParameters.softID,
               "kappaV": bonded_ia_params[self._bond_id].p.ibmVolConsParameters.kappaV}

        def _set_params_in_es_core(self):
            immersed_boundaries.volume_conservation_set_params(self._bond_id,self._params["softID"],self._params["kappaV"])
ELSE:
    class IBM_VolCons(BondedInteractionNotDefined):
        name = "IBM_VOLCONS"

class OifGlobalForces(BondedInteraction):
    """
    Part of the :ref:`Object-in-fluid` method.

    """


    def type_number(self):
        return BONDED_IA_OIF_GLOBAL_FORCES

    def type_name(self):
        """Name of interaction type.

        """
        return "OIF_GLOBAL_FORCES"

    def valid_keys(self):
        """All parameters that can be set.

        """
        return "A0_g", "ka_g", "V0", "kv"

    def required_keys(self):
        """Parameters that have to be set.

        """
        return "A0_g", "ka_g", "V0", "kv"

    def set_default_params(self):
        """Sets parameters that are not required to their default value.

        """
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


class OifLocalForces(BondedInteraction):
    """
    Part of the :ref:`Object-in-fluid` method.

    """

    def type_number(self):
        return BONDED_IA_OIF_LOCAL_FORCES

    def type_name(self):
        """Name of interaction type.

        """
        return "OIF_LOCAL_FORCES"

    def valid_keys(self):
        """All parameters that can be set.

        """
        return "r0", "ks", "kslin", "phi0", "kb", "A01", "A02", "kal", "kvisc"

    def required_keys(self):
        """Parameters that have to be set.

        """
        return "r0", "ks", "kslin", "phi0", "kb", "A01", "A02", "kal", "kvisc"

    def set_default_params(self):
        """Sets parameters that are not required to their default value.

        """
        self._params = {"r0": 1., "ks": 0., "kslin": 0.,
                        "phi0": 0., "kb": 0., "A01": 0., "A02": 0., "kal": 0., "kvisc": 0.}

    def _get_params_from_es_core(self):
        return \
            {"r0": bonded_ia_params[self._bond_id].p.oif_local_forces.r0,
             "ks": bonded_ia_params[self._bond_id].p.oif_local_forces.ks,
             "kslin": bonded_ia_params[self._bond_id].p.oif_local_forces.kslin,
             "phi0": bonded_ia_params[self._bond_id].p.oif_local_forces.phi0,
             "kb": bonded_ia_params[self._bond_id].p.oif_local_forces.kb,
             "A01": bonded_ia_params[self._bond_id].p.oif_local_forces.A01,
             "A02": bonded_ia_params[self._bond_id].p.oif_local_forces.A02,
             "kal": bonded_ia_params[self._bond_id].p.oif_local_forces.kal,
             "kvisc": bonded_ia_params[self._bond_id].p.oif_local_forces.kvisc}

    def _set_params_in_es_core(self):
        oif_local_forces_set_params(
            self._bond_id, self._params["r0"], self._params["ks"], self._params["kslin"], self._params["phi0"], self._params["kb"], self._params["A01"], self._params["A02"], self._params["kal"], self._params["kvisc"])

IF MEMBRANE_COLLISION == 1:
    class OifOutDirection(BondedInteraction):
        def type_number(self):
            return BONDED_IA_OIF_OUT_DIRECTION

        def type_name(self):
            return "OIF_OUT_DIRECTION"

        def valid_keys(self):
            return ""

        def required_keys(self):
            return ""

        def set_default_params(self):
            self._params = {}

        def _get_params_from_es_core(self):
            return \
                {}

        def _set_params_in_es_core(self):
            oif_out_direction_set_params(
                self._bond_id)

ELSE:
    class OifOutDirection(BondedInteractionNotDefined):
        name = "OIF_OUT_DIRECTION"

bonded_interaction_classes = {
    int(BONDED_IA_FENE): FeneBond,
    int(BONDED_IA_HARMONIC): HarmonicBond,
    int(BONDED_IA_HARMONIC_DUMBBELL): HarmonicDumbbellBond,
    int(BONDED_IA_RIGID_BOND): RigidBond,
    int(BONDED_IA_DIHEDRAL): Dihedral,
    int(BONDED_IA_TABULATED): Tabulated,
    int(BONDED_IA_VIRTUAL_BOND): Virtual,
    int(BONDED_IA_ANGLE_HARMONIC): AngleHarmonic,
    int(BONDED_IA_ANGLE_COSINE): AngleCosine,
    int(BONDED_IA_ANGLE_COSSQUARE): AngleCossquare,
    int(BONDED_IA_OIF_GLOBAL_FORCES): OifGlobalForces,
    int(BONDED_IA_OIF_LOCAL_FORCES): OifLocalForces,
    int(BONDED_IA_OIF_OUT_DIRECTION): OifOutDirection,
    int(BONDED_IA_IBM_TRIEL): IBM_Triel,
    int(BONDED_IA_IBM_TRIBEND): IBM_Tribend,
    int(BONDED_IA_IBM_VOLUME_CONSERVATION): IBM_VolCons,
    int(BONDED_IA_THERMALIZED_DIST): ThermalizedBond
}
IF LENNARD_JONES:
    bonded_interaction_classes[int(BONDED_IA_SUBT_LJ)] = SubtLJ
IF ELECTROSTATICS:
    bonded_interaction_classes[int(BONDED_IA_BONDED_COULOMB)] = BondedCoulomb
IF P3M:
    bonded_interaction_classes[int(BONDED_IA_BONDED_COULOMB_P3M_SR)] = BondedCoulombP3MSRBond

class BondedInteractions(object):
    """Represents the bonded interactions.

    Individual interactions can be accessed using
    BondedInteractions[i], where i is the bond id. Will return a bonded interaction
    from bonded_interaction_classes"""

    def __getitem__(self, key):
        if not is_valid_type(key, int):
            raise ValueError(
                "Index to BondedInteractions[] has to be an integer referring to a bond id")

        # Find out the type of the interaction from Espresso
        if key >= bonded_ia_params.size():
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
        if not is_valid_type(key, int):
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
        return bonded_ia_params.size()

    # Support iteration over active bonded interactions
    def __iter__(self):
        for i in range(bonded_ia_params.size()):
            if bonded_ia_params[i].type != -1:
                yield self[i]

    def add(self, bonded_ia):
        """Add a bonded ia to the simulation>"""
        self[bonded_ia_params.size()] = bonded_ia

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
