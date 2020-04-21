#
# Copyright (C) 2013-2019 The ESPResSo project
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
from cython.operator cimport dereference
include "myconfig.pxi"
from .actors cimport Actor
from .grid cimport box_geo
import numpy as np
IF SCAFACOS == 1:
    from .scafacos import ScafacosConnector
    from . cimport scafacos
from .utils cimport handle_errors
from .utils import is_valid_type, check_type_or_throw_except, to_str
from . cimport checks
from .analyze cimport partCfg, PartCfg
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

    cdef class ElectrostaticInteraction(Actor):
        def _tune(self):
            raise Exception(
                "Subclasses of ElectrostaticInteraction must define the "
                "_tune() method or chosen method does not support tuning.")

        def _set_params_in_es_core(self):
            raise Exception(
                "Subclasses of ElectrostaticInteraction must define the "
                "_set_params_in_es_core() method.")

        def _deactivate_method(self):
            deactivate_method()
            handle_errors("Coulomb method deactivation")

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
        Electrostatics solver based on the Debye-Hueckel framework. See
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
            if self._params["prefactor"] <= 0:
                raise ValueError("prefactor should be a positive float")
            if self._params["kappa"] < 0:
                raise ValueError("kappa should be a non-negative double")
            if self._params["r_cut"] < 0:
                raise ValueError("r_cut should be a non-negative double")

        def valid_keys(self):
            return ["prefactor", "kappa", "r_cut", "check_neutrality"]

        def required_keys(self):
            return ["prefactor", "kappa", "r_cut"]

        def _set_params_in_es_core(self):
            set_prefactor(self._params["prefactor"])
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

    cdef class ReactionField(ElectrostaticInteraction):
        """
        Electrostatics solver based on the Reaction-Field framework.

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
            Cut off radius for this interaction.

        """

        def validate_params(self):
            if self._params["prefactor"] <= 0:
                raise ValueError("prefactor should be a positive float")
            if self._params["kappa"] < 0:
                raise ValueError("kappa should be a non-negative double")
            if self._params["epsilon1"] < 0:
                raise ValueError("epsilon1 should be a non-negative double")
            if self._params["epsilon2"] < 0:
                raise ValueError("epsilon2 should be a non-negative double")
            if self._params["r_cut"] < 0:
                raise ValueError("r_cut should be a non-negative double")

        def valid_keys(self):
            return ["prefactor", "kappa", "epsilon1", "epsilon2", "r_cut",
                    "check_neutrality"]

        def required_keys(self):
            return ["prefactor", "kappa", "epsilon1", "epsilon2", "r_cut"]

        def _set_params_in_es_core(self):
            set_prefactor(self._params["prefactor"])
            rf_set_params(
                self._params["kappa"],
                self._params["epsilon1"],
                self._params["epsilon2"],
                self._params["r_cut"])

        def _get_params_from_es_core(self):
            params = {}
            params.update(rf_params)
            return params

        def _activate_method(self):
            check_neutrality(self._params)
            coulomb.method = COULOMB_RF
            self._set_params_in_es_core()

        def default_params(self):
            return {"prefactor": -1,
                    "kappa": -1,
                    "epsilon1": -1,
                    "epsilon2": -1,
                    "r_cut": -1,
                    "check_neutrality": True}


IF P3M == 1:
    cdef class P3M(ElectrostaticInteraction):
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
            The charge-assignment order, an integer between 0 and 7.
        epsilon : :obj:`float` or :obj:`str`, optional
            A positive number for the dielectric constant of the
            surrounding medium. Use ``'metallic'`` to set the dielectric
            constant of the surrounding medium to infinity (default).
        mesh : :obj:`int` or (3,) array_like of :obj:`int`, optional
            The number of mesh points in x, y and z direction. Use a single
            value for cubic boxes.
        r_cut : :obj:`float`, optional
            The real space cutoff.
        tune : :obj:`bool`, optional
            Used to activate/deactivate the tuning method on activation.
            Defaults to ``True``.
        check_neutrality : :obj:`bool`, optional
            Raise a warning if the system is not electrically neutral when
            set to ``True`` (default).

        """

        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)

        def validate_params(self):
            default_params = self.default_params()
            if not (self._params["prefactor"] > 0.0):
                raise ValueError("prefactor should be a positive float")

            if not (self._params["r_cut"] >= 0
                    or self._params["r_cut"] == default_params["r_cut"]):
                raise ValueError("P3M r_cut has to be >=0")

            if is_valid_type(self._params["mesh"], int):
                if self._params["mesh"] % 2 != 0 and self._params["mesh"] != -1:
                    raise ValueError(
                        "P3M requires an even number of mesh points in all directions")
            else:
                check_type_or_throw_except(self._params["mesh"], 3, int,
                                           "P3M mesh has to be an integer or integer list of length 3")
                if (self._params["mesh"][0] % 2 != 0 and self._params["mesh"][0] != -1) or \
                   (self._params["mesh"][1] % 2 != 0 and self._params["mesh"][1] != -1) or \
                   (self._params["mesh"][2] % 2 != 0 and self._params["mesh"][2] != -1):
                    raise ValueError(
                        "P3M requires an even number of mesh points in all directions")

            if not (self._params["cao"] >= -1 and self._params["cao"] <= 7):
                raise ValueError(
                    "P3M cao has to be an integer between -1 and 7")

            if self._params["tune"] and not (self._params["accuracy"] >= 0):
                raise ValueError("P3M accuracy has to be positive")

            if self._params["epsilon"] == "metallic":
                self._params = 0.0

            if not (is_valid_type(self._params["epsilon"], float)
                    or self._params["epsilon"] == "metallic"):
                raise ValueError("epsilon should be a double or 'metallic'")

            if self._params["mesh_off"] != default_params["mesh_off"]:
                check_type_or_throw_except(self._params["mesh_off"], 3, float,
                                           "mesh_off should be a (3,) array_like of values between 0.0 and 1.0")

            if not (self._params["alpha"] == default_params["alpha"]
                    or self._params["alpha"] > 0):
                raise ValueError(
                    "alpha should be positive")

        def valid_keys(self):
            return ["mesh", "cao", "accuracy", "epsilon", "alpha", "r_cut",
                    "prefactor", "tune", "check_neutrality"]

        def required_keys(self):
            return ["prefactor", "accuracy"]

        def default_params(self):
            return {"cao": 0,
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
            # Sets lb, bcast, resets vars to zero if lb=0
            set_prefactor(self._params["prefactor"])
            # Sets cdef vars and calls p3m_set_params() in core
            python_p3m_set_params(self._params["r_cut"],
                                  self._params["mesh"], self._params["cao"],
                                  self._params["alpha"], self._params["accuracy"])
            # p3m_set_params()  -> set r_cuts, mesh, cao, validates sanity, bcasts
            # Careful: bcast calls on_coulomb_change(), which calls p3m_init(),
            #         which resets r_cut if lb is zero. OK.
            # Sets eps, bcast
            p3m_set_eps(self._params["epsilon"])
            python_p3m_set_mesh_offset(self._params["mesh_off"])

        def _tune(self):
            set_prefactor(self._params["prefactor"])
            python_p3m_set_tune_params(self._params["r_cut"],
                                       self._params["mesh"],
                                       self._params["cao"],
                                       -1.0,
                                       self._params["accuracy"])
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
            r_cut : :obj:`float`, optional
                The real space cutoff
            tune : :obj:`bool`, optional
                Used to activate/deactivate the tuning method on activation.
                Defaults to ``True``.
            check_neutrality : :obj:`bool`, optional
                Raise a warning if the system is not electrically neutral when
                set to ``True`` (default).

            """

            def __init__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)

            def validate_params(self):
                default_params = self.default_params()

                if not (self._params["r_cut"] >= 0
                        or self._params["r_cut"] == default_params["r_cut"]):
                    raise ValueError("P3M r_cut has to be >=0")

                if is_valid_type(self._params["mesh"], int):
                    if self._params["mesh"] % 2 != 0 and self._params["mesh"] != -1:
                        raise ValueError(
                            "P3M requires an even number of mesh points in all directions")
                else:
                    check_type_or_throw_except(self._params["mesh"], 3, int,
                                               "P3M mesh has to be an integer or integer list of length 3")
                    if (self._params["mesh"][0] % 2 != 0 and self._params["mesh"][0] != -1) or \
                       (self._params["mesh"][1] % 2 != 0 and self._params["mesh"][1] != -1) or \
                       (self._params["mesh"][2] % 2 != 0 and self._params["mesh"][2] != -1):
                        raise ValueError(
                            "P3M requires an even number of mesh points in all directions")

                if not (self._params["cao"] >= -1
                        and self._params["cao"] <= 7):
                    raise ValueError(
                        "P3M cao has to be an integer between -1 and 7")

                if not (self._params["accuracy"] >= 0):
                    raise ValueError("P3M accuracy has to be positive")

                # if self._params["epsilon"] == "metallic":
                #  self._params = 0.0

                if not (is_valid_type(self._params["epsilon"], float)
                        or self._params["epsilon"] == "metallic"):
                    raise ValueError(
                        "epsilon should be a double or 'metallic'")

                if self._params["mesh_off"] != default_params["mesh_off"]:
                    check_type_or_throw_except(self._params["mesh_off"], 3, float,
                                               "mesh_off should be a (3,) array_like of values between 0.0 and 1.0")

            def valid_keys(self):
                return ["mesh", "cao", "accuracy", "epsilon", "alpha", "r_cut",
                        "prefactor", "tune", "check_neutrality"]

            def required_keys(self):
                return ["prefactor", "accuracy"]

            def default_params(self):
                return {"cao": 0,
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
                set_prefactor(self._params["prefactor"])
                python_p3m_set_tune_params(self._params["r_cut"],
                                           self._params["mesh"],
                                           self._params["cao"],
                                           -1.0,
                                           self._params["accuracy"])
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
                set_prefactor(self._params["prefactor"])
                python_p3m_set_params(self._params["r_cut"],
                                      self._params["mesh"],
                                      self._params["cao"],
                                      self._params["alpha"],
                                      self._params["accuracy"])
                p3m_set_eps(self._params["epsilon"])
                python_p3m_set_mesh_offset(self._params["mesh_off"])
                handle_errors("p3m gpu init")

IF ELECTROSTATICS:
    cdef class MMM1D(ElectrostaticInteraction):
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
        bessel_cutoff : :obj:`int`, optional
        tune : :obj:`bool`, optional
            Specify whether to automatically tune or not. Defaults to ``True``.

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
            return ["prefactor", "maxPWerror", "far_switch_radius",
                    "bessel_cutoff", "tune", "check_neutrality"]

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
            set_prefactor(self._params["prefactor"])
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
        tune : :obj:`bool`, optional
            Specify whether to automatically tune or not. Defaults to ``True``.
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
            dereference(self.thisptr).activate()

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
            return ["prefactor", "maxPWerror", "far_switch_radius",
                    "bessel_cutoff", "tune", "check_neutrality"]

        def required_keys(self):
            return ["prefactor", "maxPWerror"]

        def _get_params_from_es_core(self):
            params = {}
            params.update(mmm1d_params)
            params["prefactor"] = coulomb.prefactor
            return params

        def _set_params_in_es_core(self):
            set_prefactor(self._params["prefactor"])
            default_params = self.default_params()

            self.thisptr.set_params(
                box_geo.length()[2], coulomb.prefactor,
                self._params["maxPWerror"], self._params["far_switch_radius"],
                self._params["bessel_cutoff"])

        def _tune(self):
            self.thisptr.setup(dereference(self.interface))
            self.thisptr.tune(
                dereference(self.interface), self._params["maxPWerror"],
                self._params["far_switch_radius"], self._params["bessel_cutoff"])

        def _activate_method(self):
            check_neutrality(self._params)
            self._set_params_in_es_core()
            coulomb.method = COULOMB_MMM1D_GPU
            if self._params["tune"]:
                self._tune()
            self._set_params_in_es_core()

        def _deactivate_method(self):
            dereference(self.thisptr).deactivate()

IF ELECTROSTATICS:
    IF SCAFACOS == 1:
        class Scafacos(ScafacosConnector, ElectrostaticInteraction):

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
            """

            dipolar = False

            # Explicit constructor needed due to multiple inheritance
            def __init__(self, *args, **kwargs):
                Actor.__init__(self, *args, **kwargs)

            def _activate_method(self):
                check_neutrality(self._params)
                set_prefactor(self._params["prefactor"])
                self._set_params_in_es_core()
                mpi_bcast_coulomb_params()

            def default_params(self):
                return {}

            def _deactivate_method(self):
                super()._deactivate_method()
                scafacos.free_handle()
                mpi_bcast_coulomb_params()
