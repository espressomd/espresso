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
from cython.operator cimport dereference
include "myconfig.pxi"
from espressomd cimport actors
from . import actors
cimport grid
cimport globals
import numpy as np
IF SCAFACOS == 1:
    from .scafacos import ScafacosConnector
    from . cimport scafacos
from espressomd.utils cimport handle_errors
from espressomd.utils import is_valid_type
from . cimport checks
from .c_analyze cimport partCfg, PartCfg
from .particle_data cimport particle


IF ELECTROSTATICS == 1:
    def check_neutrality(_params):
        if "check_neutrality" in _params:
            if(_params["check_neutrality"]):
                if not checks.check_charge_neutrality[PartCfg](partCfg()):
                    raise Exception("""
                    The system is not charge neutral. Please
                    neutralize the system before adding a new actor by adding
                    the corresponding counterions to the system. Alternatively
                    you can turn off the electroneutrality check by supplying
                    check_neutrality=False when creating the actor. In this
                    case you may be simulating a non-neutral system which will
                    affect physical observables like e.g. the pressure, the
                    chemical potentials of charged species or potential
                    energies of the system. Since simulations of non charge
                    neutral systems are special please make sure you know what
                    you are doing.
                    """)

    cdef class ElectrostaticInteraction(actors.Actor):
        def _tune(self):
            raise Exception(
                "Subclasses of ElectrostaticInteraction must define the _tune() \
                method or chosen method does not support tuning.")

        def _set_params_in_es_core(self):
            raise Exception(
                "Subclasses of ElectrostaticInteraction must define the \
                _set_params_in_es_core() method.")

        def _deactivate_method(self):
            deactivate_coulomb_method()
            handle_errors("Coulom method deactivation")

        def tune(self, **tune_params_subset):
            if tune_params_subset is not None:
                if all(k in self.valid_keys() for k in tune_params_subset):
                    self._params.update(tune_params_subset)
                else:
                    raise ValueError(
                        "Invalid parameter given to tune function.")
            self._tune()


IF ELECTROSTATICS:
    cdef class DH(ElectrostaticInteraction):
        """
        Solve electrostatics in the Debye-Hueckel framework see
        :ref:`Debye-Hückel potential` for more details.

        Parameters
        ----------
        prefactor : :obj:`float`
            Electrostatics prefactor (see :eq:`coulomb_prefactor`).
        kappa : :obj:`float`
            Inverse Debye screening length.
        r_cut : :obj:`float`
            Cut off radius for this interaction.

        """

        def validate_params(self):
            if (self._params["prefactor"] <= 0):
                raise ValueError(
                    "prefactor should be a positive float")
            if (self._params["kappa"] < 0):
                raise ValueError("kappa should be a non-negative double")
            if (self._params["r_cut"] < 0):
                raise ValueError("r_cut should be a non-negative double")

        def valid_keys(self):
            return "prefactor", "kappa", "r_cut", "check_neutrality"

        def required_keys(self):
            return "prefactor", "kappa", "r_cut"

        def _set_params_in_es_core(self):
            coulomb_set_prefactor(self._params["prefactor"])
            dh_set_params(self._params["kappa"], self._params["r_cut"])

        def _get_params_from_es_core(self):
            params = {}
            params.update(dh_params)
            return params

        def _activate_method(self):
            check_neutrality(self._params)
            coulomb.method = COULOMB_DH
            self._set_params_in_es_core()

        def default_params(self):
            return {"prefactor": -1,
                    "kappa": -1,
                    "r_cut": -1,
                    "check_neutrality": True}


IF P3M == 1:
    cdef class P3M(ElectrostaticInteraction):

        def __init__(self, *args, **kwargs):
            """
            P3M electrostatics solver.

            Particle–Particle-Particle–Mesh (P3M) is a Fourier-based Ewald
            summation method to calculate potentials in N-body simulation.

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
            epsilon : :obj:`str`, optional
                Use 'metallic' to set the dielectric constant of the
                surrounding medium to infinity (Default).
            epsilon : :obj:`float`, optional
                A positive number for the dielectric constant of the
                surrounding medium.
            mesh : :obj:`int`, optional
                The number of mesh points.
            mesh : array_like, optional
                The number of mesh points in x, y and z direction. This is
                relevant for noncubic boxes.
            r_cut : :obj:`float`, optional
                The real space cutoff.
            tune : :obj:`bool`, optional
                Used to activate/deactivate the tuning method on activation.
                Defaults to True.

            """
            super(type(self), self).__init__(*args, **kwargs)

        def validate_params(self):
            default_params = self.default_params()
            if not (self._params["prefactor"] > 0.0):
                raise ValueError("prefactor should be a positive float")

            if not (self._params["r_cut"] >= 0 or self._params["r_cut"] == default_params["r_cut"]):
                raise ValueError("P3M r_cut has to be >=0")

            if not (is_valid_type(self._params["mesh"], int) or len(self._params["mesh"])):
                raise ValueError(
                    "P3M mesh has to be an integer or integer list of length 3")

            if (isinstance(self._params["mesh"], basestring) and len(self._params["mesh"]) == 3):
                if (self._params["mesh"][0] % 2 != 0 and self._params["mesh"][0] != -1) or (self._params["mesh"][1] % 2 != 0 and self._params["mesh"][1] != -1) or (self._params["mesh"][2] % 2 != 0 and self._params["mesh"][2] != -1):
                    raise ValueError(
                        "P3M requires an even number of mesh points in all directions")

            if not (self._params["cao"] >= -1 and self._params["cao"] <= 7):
                raise ValueError(
                    "P3M cao has to be an integer between -1 and 7")

            if self._params["tune"] and not (self._params["accuracy"] >= 0):
                raise ValueError("P3M accuracy has to be positive")

            if self._params["epsilon"] == "metallic":
                self._params = 0.0

            if not (is_valid_type(self._params["epsilon"], float) or self._params["epsilon"] == "metallic"):
                raise ValueError("epsilon should be a double or 'metallic'")

            if not (self._params["inter"] == default_params["inter"] or self._params["inter"] > 0):
                raise ValueError("inter should be a positive integer")

            if not (self._params["mesh_off"] == default_params["mesh_off"] or len(self._params) != 3):
                raise ValueError(
                    "mesh_off should be a list of length 3 and values between 0.0 and 1.0")

            if not (self._params["alpha"] == default_params["alpha"] or self._params["alpha"] > 0):
                raise ValueError(
                    "alpha should be positive")

        def valid_keys(self):
            return "mesh", "cao", "accuracy", "epsilon", "alpha", "r_cut", "prefactor", "tune", "check_neutrality"

        def required_keys(self):
            return ["prefactor", "accuracy"]

        def default_params(self):
            return {"cao": 0,
                    "inter": -1,
                    "r_cut": -1,
                    "alpha": 0,
                    "accuracy": 0,
                    "mesh": [0, 0, 0],
                    "epsilon": 0.0,
                    "mesh_off": [-1, -1, -1],
                    "tune": True,
                    "check_neutrality": True}

        def _get_params_from_es_core(self):
            params = {}
            params.update(p3m.params)
            params["prefactor"] = coulomb.prefactor
            params["tune"] = self._params["tune"]
            return params

        def _set_params_in_es_core(self):
            #Sets lb, bcast, resets vars to zero if lb=0
            coulomb_set_prefactor(self._params["prefactor"])
            #Sets cdef vars and calls p3m_set_params() in core
            python_p3m_set_params(self._params["r_cut"],
                                  self._params["mesh"], self._params["cao"],
                                  self._params["alpha"], self._params["accuracy"])
            #p3m_set_params()  -> set r_cuts, mesh, cao, validates sanity, bcasts
            #Careful: bcast calls on_coulomb_change(), which calls p3m_init(),
            #         which resets r_cut if lb is zero. OK.
            #Sets eps, bcast
            p3m_set_eps(self._params["epsilon"])
            #Sets ninterpol, bcast
            p3m_set_ninterpol(self._params["inter"])
            python_p3m_set_mesh_offset(self._params["mesh_off"])

        def _tune(self):
            coulomb_set_prefactor(self._params["prefactor"])
            python_p3m_set_tune_params(self._params["r_cut"],
                                       self._params["mesh"],
                                       self._params["cao"],
                                       -1.0,
                                       self._params["accuracy"],
                                       self._params["inter"])
            resp = python_p3m_adaptive_tune()
            if resp:
                raise Exception(
                    "failed to tune P3M parameters to required accuracy")
            self._params.update(self._get_params_from_es_core())

        def _activate_method(self):
            check_neutrality(self._params)
            if self._params["tune"]:
                self._tune()
            self._set_params_in_es_core()

    IF CUDA:
        cdef class P3MGPU(ElectrostaticInteraction):

            def __init__(self, *args, **kwargs):
                """
                P3M electrostatics solver with GPU support.

                Particle–Particle-Particle–Mesh (P3M) is a Fourier-based Ewald
                summation method to calculate potentials in N-body simulation.

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
                epsilon : :obj:`str`, optional
                    Use 'metallic' to set the dielectric constant of the
                    surrounding medium to infinity (Default).
                epsilon : :obj:`float`, optional
                    A positive number for the dielectric constant of the
                    surrounding medium.
                mesh : :obj:`int`, optional
                    The number of mesh points.
                mesh : array_like, optional
                    The number of mesh points in x, y and z direction. This is
                    relevant for noncubic boxes.
                r_cut : :obj:`float`, optional
                    The real space cutoff
                tune : :obj:`bool`, optional
                    Used to activate/deactivate the tuning method on activation.
                    Defaults to True.

                """
                super(type(self), self).__init__(*args, **kwargs)

            def validate_params(self):
                default_params = self.default_params()

                if not (self._params["r_cut"] >= 0 or self._params["r_cut"] == default_params["r_cut"]):
                    raise ValueError("P3M r_cut has to be >=0")

                if not (is_valid_type(self._params["mesh"], int) or len(self._params["mesh"])):
                    raise ValueError(
                        "P3M mesh has to be an integer or integer list of length 3")

                if (isinstance(self._params["mesh"], basestring) and len(self._params["mesh"]) == 3):
                    if (self._params["mesh"][0] % 2 != 0 and self._params["mesh"][0] != -1) or (self._params["mesh"][1] % 2 != 0 and self._params["mesh"][1] != -1) or (self._params["mesh"][2] % 2 != 0 and self._params["mesh"][2] != -1):
                        raise ValueError(
                            "P3M requires an even number of mesh points in all directions")

                if not (self._params["cao"] >= -1 and self._params["cao"] <= 7):
                    raise ValueError(
                        "P3M cao has to be an integer between -1 and 7")

                if not (self._params["accuracy"] >= 0):
                    raise ValueError("P3M accuracy has to be positive")

                # if self._params["epsilon"] == "metallic":
                #  self._params = 0.0

                if not (is_valid_type(self._params["epsilon"], float) or self._params["epsilon"] == "metallic"):
                    raise ValueError(
                        "epsilon should be a double or 'metallic'")

                if not (self._params["inter"] == default_params["inter"] or self._params["inter"] > 0):
                    raise ValueError("inter should be a positive integer")

                if not (self._params["mesh_off"] == default_params["mesh_off"] or len(self._params) != 3):
                    raise ValueError(
                        "mesh_off should be a list of length 3 with values between 0.0 and 1.0")

            def valid_keys(self):
                return "mesh", "cao", "accuracy", "epsilon", "alpha", "r_cut", "prefactor", "tune", "check_neutrality"

            def required_keys(self):
                return ["prefactor", "accuracy"]

            def default_params(self):
                return {"cao": 0,
                        "inter": -1,
                        "r_cut": -1,
                        "alpha": 0,
                        "accuracy": 0,
                        "mesh": [0, 0, 0],
                        "epsilon": 0.0,
                        "mesh_off": [-1, -1, -1],
                        "tune": True,
                        "check_neutrality": True}

            def _get_params_from_es_core(self):
                params = {}
                params.update(p3m.params)
                params["prefactor"] = coulomb.prefactor
                params["tune"] = self._params["tune"]
                return params

            def _tune(self):
                coulomb_set_prefactor(self._params["prefactor"])
                python_p3m_set_tune_params(self._params["r_cut"],
                                           self._params["mesh"],
                                           self._params["cao"],
                                           -1.0,
                                           self._params["accuracy"],
                                           self._params["inter"])
                resp = python_p3m_adaptive_tune()
                if resp:
                    raise Exception(
                        "failed to tune P3M parameters to required accuracy")
                self._params.update(self._get_params_from_es_core())

            def _activate_method(self):
                check_neutrality(self._params)
                python_p3m_gpu_init(self._params)
                coulomb.method = COULOMB_P3M_GPU
                if self._params["tune"]:
                    self._tune()
                python_p3m_gpu_init(self._params)
                self._set_params_in_es_core()

            def _set_params_in_es_core(self):
                coulomb_set_prefactor(self._params["prefactor"])
                python_p3m_set_params(self._params["r_cut"],
                                      self._params["mesh"],
                                      self._params["cao"],
                                      self._params["alpha"],
                                      self._params["accuracy"])
                p3m_set_eps(self._params["epsilon"])
                p3m_set_ninterpol(self._params["inter"])
                python_p3m_set_mesh_offset(self._params["mesh_off"])
                handle_errors("p3m gpu init")

IF ELECTROSTATICS:
    cdef class MMM1D(ElectrostaticInteraction):
        """
        Electrostatics solver for Systems with one periodic direction.
        See :ref:`MMM1D Theory` for more details.

        Parameters
        ----------
        prefactor : :obj:`float`
            Electrostatics prefactor (see :eq:`coulomb_prefactor`).
        maxWPerror : :obj:`float`
            Maximal pairwise error.
        far_switch_radius : :obj:`float`, optional
            Radius where near-field and far-field calculation are switched.
        bessel_cutoff : :obj:`int`, optional
        tune : :obj:`bool`, optional
            Specify whether to automatically tune ore not. The default is True.
        """

        def validate_params(self):
            default_params = self.default_params()
            if self._params["prefactor"] <= 0.:
                raise ValueError("Prefactor should be a positive float")
            if self._params["maxPWerror"] < 0 and self._params["maxPWerror"] != default_params["maxPWerror"]:
                raise ValueError("maxPWerror should be a positive double")
            if self._params["far_switch_radius"] < 0 and self._params["far_switch_radius"] != default_params["far_switch_radius"]:
                raise ValueError("switch radius should be a positive double")
            if self._params["bessel_cutoff"] < 0 and self._params["bessel_cutoff"] != default_params["bessel_cutoff"]:
                raise ValueError("bessel_cutoff should be a positive integer")

        def default_params(self):
            return {"prefactor": -1,
                    "maxPWerror": -1,
                    "far_switch_radius": -1,
                    "bessel_cutoff": -1,
                    "tune": True,
                    "check_neutrality": True}

        def valid_keys(self):
            return "prefactor", "maxPWerror", "far_switch_radius", "bessel_cutoff", "tune", "check_neutrality"

        def required_keys(self):
            return ["prefactor", "maxPWerror"]

        def _get_params_from_es_core(self):
            params = {}
            params.update(mmm1d_params)
            params["far_switch_radius"] = np.sqrt(
                params["far_switch_radius_2"])
            del params["far_switch_radius_2"]
            params["prefactor"] = coulomb.prefactor
            params["tune"] = self._params["tune"]
            return params

        def _set_params_in_es_core(self):
            coulomb_set_prefactor(self._params["prefactor"])
            MMM1D_set_params(
                self._params["far_switch_radius"], self._params["maxPWerror"])

        def _tune(self):
            cdef int resp
            resp = pyMMM1D_tune()
            if resp:
                raise Exception("failed to tune mmm1d ")
            self._params.update(self._get_params_from_es_core())

        def _activate_method(self):
            check_neutrality(self._params)
            coulomb.method = COULOMB_MMM1D
            self._set_params_in_es_core()
            if self._params["tune"]:
                self._tune()

            self._set_params_in_es_core()

IF ELECTROSTATICS and MMM1D_GPU:
    cdef class MMM1DGPU(ElectrostaticInteraction):
        """
        Electrostatics solver for Systems with one periodic direction.
        See :ref:`MMM1D Theory` for more details.

        Parameters
        ----------
        prefactor : :obj:`float`
            Electrostatics prefactor (see :eq:`coulomb_prefactor`).
        maxWPerror : :obj:`float`
            Maximal pairwise error.
        far_switch_radius : :obj:`float`, optional
            Radius where near-field and far-field calculation are switched
        bessel_cutoff : :obj:`int`, optional
        tune : :obj:`bool`, optional
            Specify whether to automatically tune ore not. The default is True.
        """
        cdef Mmm1dgpuForce * thisptr
        cdef EspressoSystemInterface * interface
        cdef char * log
        cdef int resp

        def __cinit__(self):
            self.interface = EspressoSystemInterface._Instance()
            default_params = self.default_params()
            self.thisptr = new Mmm1dgpuForce(dereference(self.interface), 0.0, default_params["maxPWerror"])
            self.interface.update()
            self.interface.requestRGpu()

        def __dealloc__(self):
            del self.thisptr

        def validate_params(self):
            default_params = self.default_params()
            if self._params["prefactor"] <= 0:
                raise ValueError("prefactor should be a positive float")
            if self._params["maxPWerror"] < 0 and self._params["maxPWerror"] != default_params["maxPWerror"]:
                raise ValueError("maxPWerror should be a positive double")
            if self._params["far_switch_radius"] < 0 and self._params["far_switch_radius"] != default_params["far_switch_radius"]:
                raise ValueError("switch radius should be a positive double")
            if self._params["bessel_cutoff"] < 0 and self._params["bessel_cutoff"] != default_params["bessel_cutoff"]:
                raise ValueError("bessel_cutoff should be a positive integer")

        def default_params(self):
            return {"prefactor": -1,
                    "maxPWerror": -1.0,
                    "far_switch_radius": -1.0,
                    "bessel_cutoff": -1,
                    "tune": True,
                    "check_neutrality": True}

        def valid_keys(self):
            return "prefactor", "maxPWerror", "far_switch_radius", "bessel_cutoff", "tune", "check_neutrality"

        def required_keys(self):
            return ["prefactor", "maxPWerror"]

        def _get_params_from_es_core(self):
            params = {}
            params.update(mmm1d_params)
            params["prefactor"] = coulomb.prefactor
            return params

        def _set_params_in_es_core(self):
            coulomb_set_prefactor(self._params["prefactor"])
            default_params = self.default_params()

            self.thisptr.set_params(grid.box_l[2], coulomb.prefactor, self._params[
                                    "maxPWerror"], self._params["far_switch_radius"], self._params["bessel_cutoff"])

        def _tune(self):
            self.thisptr.setup(dereference(self.interface))
            self.thisptr.tune(dereference(self.interface), self._params[
                              "maxPWerror"], self._params["far_switch_radius"], self._params["bessel_cutoff"])

        def _activate_method(self):
            check_neutrality(self._params)
            self._set_params_in_es_core()
            coulomb.method = COULOMB_MMM1D_GPU
            if self._params["tune"]:
                self._tune()
            self._set_params_in_es_core()

IF ELECTROSTATICS:
    cdef class MMM2D(ElectrostaticInteraction):
        """
        Electrostatics solver for systems with two periodic dimensions.
        More detail are in the user guide :ref:`MMM2D Theory`

        Parameters
        ----------
        prefactor : :obj:`float`
            Electrostatics prefactor (see :eq:`coulomb_prefactor`).
        maxWPerror : :obj:`float`
            Maximal pairwise error.
        dielectric : :obj:`int`, optional
            Selector parameter for setting the dielectric constants manually
            (top, mid, bottom), mutually exclusive with dielectric-contrast
        top : :obj:`float`, optional
            If dielectric is specified this parameter sets the dielectric
            constant *above* the simulation box
            :math:`\\varepsilon_\\mathrm{top}`
        mid : :obj:`float`, optional
            If dielectric is specified this parameter sets the dielectric
            constant *in* the simulation box :math:`\\varepsilon_\\mathrm{mid}`.
        bottom : :obj:`float`, optional
            If dielectric is specified this parameter sets the dielectric
            constant *below* the simulation box
            :math:`\\varepsilon_\\mathrm{bot}`.
        dielectric_contrast_on : :obj:`int`, optional
            Selector parameter for setting a dielectric contrast between the
            upper simulation boundary and the simulation box, and between the
            lower simulation boundary and the simulation box, respectively.
        delta_mid_top : :obj:`float`, optional
            If dielectric-contrast mode is selected, then this parameter sets
            the dielectric contrast between the upper boundary and the
            simulation box :math:`\\Delta_t`.
        delta_mid_bottom : :obj:`float`, optional
            If dielectric-contrast mode is selected, then this parameter sets
            the dielectric contrast between the lower boundary and the
            simulation box :math:`\\Delta_b`.
        const_pot : :obj:`int`, optional
            Selector parameter for setting a constant electric potential
            between the top and bottom of the simulation box.
        pot_diff : :obj:`float`, optional
            If const_pot mode is selected this parameter controls the applied
            voltage.
        far_cut : :obj:`float`, optional
            Cut off radius, use with care, intended for testing purposes.

        """

        def validate_params(self):
            default_params = self.default_params()
            if self._params["prefactor"] <= 0:
                raise ValueError("prefactor should be a positive float")
            if self._params["maxPWerror"] < 0 and self._params["maxPWerror"] != default_params["maxPWerror"]:
                raise ValueError("maxPWerror should be a positive double")
            if self._params["dielectric"] == 1 and (self._params["top"] < 0 or self._params["mid"] < 0 or self._params["bot"] < 0):
                raise ValueError("Dielectric constants should be > 0!")
            if self._params["dielectric_contrast_on"] == 1 and (self._params["delta_mid_top"] == default_params["delta_mid_top"] or self._params["delta_mid_bot"] == default_params["delta_mid_bot"]):
                raise ValueError("Dielectric constrast not set!")
            if self._params["dielectric"] == 1 and self._params["dielectric_contrast_on"] == 1:
                raise ValueError(
                    "dielectric and dielectric_contrast are mutually exclusive!")

        def default_params(self):
            return {"prefactor": -1,
                    "maxPWerror": -1,
                    "far_cut": -1,
                    "top": 0,
                    "mid": 0,
                    "bot": 0,
                    "dielectric": 0,
                    "dielectric_contrast_on": 0,
                    "const_pot": 0,
                    "delta_mid_top": 0,
                    "delta_mid_bot": 0,
                    "pot_diff": 0,
                    "check_neutrality": True}

        def required_keys(self):
            return ["prefactor", "maxPWerror"]

        def valid_keys(self):
            return "prefactor", "maxPWerror", "top", "mid", "bot", "delta_mid_top", "delta_mid_bot", "pot_diff", "dielectric", "dielectric_contrast_on", "const_pot", "far_cut", "check_neutrality"

        def _get_params_from_es_core(self):
            params = {}
            params.update(mmm2d_params)
            params["prefactor"] = coulomb.prefactor
            if params["dielectric_contrast_on"] == 1 or params["const_pot"] == 1:
                params["dielectric"] = 1
            else:
                params["dielectric"] = 0
            return params

        def _set_params_in_es_core(self):
            coulomb_set_prefactor(self._params["prefactor"])
            if self._params["dielectric"]:
                self._params["delta_mid_top"] = (self._params[
                                                 "mid"] - self._params["top"]) / (self._params["mid"] + self._params["top"])
                self._params["delta_mid_bot"] = (self._params[
                                                 "mid"] - self._params["bot"]) / (self._params["mid"] + self._params["bot"])

            if self._params["const_pot"]:
                self._params["delta_mid_top"] = -1
                self._params["delta_mid_bot"] = -1

            res = MMM2D_set_params(self._params["maxPWerror"],
                                   self._params["far_cut"], self._params[
                                       "delta_mid_top"],
                                   self._params["delta_mid_bot"], self._params["const_pot"], self._params["pot_diff"])
            handle_errors("MMM2d setup")
            if res:
                raise Exception("MMM2D setup failed")

        def _activate_method(self):
            check_neutrality(self._params)
            coulomb.method = COULOMB_MMM2D
            self._set_params_in_es_core()
            MMM2D_init()
            handle_errors("MMM2d setup")
            res = MMM2D_sanity_checks()
            handle_errors("MMM2d setup")
            if res:
                raise Exception("MMM2D sanity checks failed.")
            mpi_bcast_coulomb_params()
            handle_errors("MMM2d setup")

    IF SCAFACOS == 1:
        class Scafacos(ScafacosConnector, ElectrostaticInteraction):

            """Calculates Coulomb interactions using method from the SCAFACOs library."""
            dipolar = False

            # Explicit constructor needed due to multiple inheritance
            def __init__(self, *args, **kwargs):
                actors.Actor.__init__(self, *args, **kwargs)

            def _activate_method(self):
                check_neutrality(self._params)
                coulomb_set_prefactor(self._params["prefactor"])
                self._set_params_in_es_core()
                mpi_bcast_coulomb_params()

            def default_params(self):
                return {}

            def _deactivate_method(self):
                super(Scafacos, self)._deactivate_method()
                scafacos.free_handle()
                mpi_bcast_coulomb_params()
