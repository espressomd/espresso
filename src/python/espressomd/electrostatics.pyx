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
from libcpp.memory cimport shared_ptr, make_shared
from cython.operator cimport dereference
include "myconfig.pxi"
from .actors cimport Actor
from .grid cimport box_geo
import numpy as np
IF SCAFACOS == 1:
    from .scafacos import ScafacosConnector
    from . cimport scafacos
from .utils import is_valid_type, check_type_or_throw_except, handle_errors
from .analyze cimport partCfg, PartCfg
from .particle_data cimport particle
import sys


IF ELECTROSTATICS == 1:
    def check_neutrality(_params):
        if _params.get("check_neutrality", False):
            if not check_charge_neutrality[PartCfg](partCfg()):
                raise Exception(' '.join("""
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
                    """.split()))

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
            Cutoff radius for this interaction.

        """

        def validate_params(self):
            pass

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

        def validate_params(self):
            pass

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
    cdef class _P3MBase(ElectrostaticInteraction):

        cdef _check_and_copy_mesh_size(self, int mesh[3], pmesh):
            if is_valid_type(pmesh, int):
                pmesh = 3 * [pmesh]
            else:
                check_type_or_throw_except(
                    pmesh, 3, int, "mesh size must be 3 ints")
            for i in range(3):
                mesh[i] = pmesh[i]

        def valid_keys(self):
            return ["mesh", "cao", "accuracy", "epsilon", "alpha", "r_cut",
                    "prefactor", "tune", "check_neutrality", "timings",
                    "verbose", "mesh_off"]

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
                    "timings": 10,
                    "check_neutrality": True,
                    "verbose": True}

        def _get_params_from_es_core(self):
            params = {}
            params.update(p3m.params)
            params["prefactor"] = coulomb.prefactor
            params["tune"] = self._params["tune"]
            params["timings"] = self._params["timings"]
            return params

        def _tune(self):
            cdef int mesh[3]
            self._check_and_copy_mesh_size(mesh, self._params["mesh"])

            set_prefactor(self._params["prefactor"])
            p3m_set_eps(self._params["epsilon"])
            p3m_set_tune_params(self._params["r_cut"], mesh,
                                self._params["cao"], self._params["accuracy"])
            tuning_error = p3m_adaptive_tune(
                self._params["timings"], self._params["verbose"])
            if tuning_error:
                handle_errors("P3M: tuning failed")
            self._params.update(self._get_params_from_es_core())

        def tune(self, **tune_params_subset):
            # update the three necessary parameters if not provided by the user
            default_params = self.default_params()
            for key in ["r_cut", "mesh", "cao"]:
                if key not in tune_params_subset:
                    tune_params_subset[key] = default_params[key]

            super().tune(**tune_params_subset)

        def _set_params_in_es_core(self):
            cdef int mesh[3]
            self._check_and_copy_mesh_size(mesh, self._params["mesh"])

            set_prefactor(self._params["prefactor"])
            # Sets p3m parameters
            # p3m_set_params() -> set parameters and bcasts
            # Careful: calls on_coulomb_change(), which calls p3m_init(),
            #          which resets r_cut if prefactor=0
            p3m_set_params(self._params["r_cut"], mesh, self._params["cao"],
                           self._params["alpha"], self._params["accuracy"])
            # Sets eps, bcast
            p3m_set_eps(self._params["epsilon"])
            p3m_set_mesh_offset(self._params["mesh_off"][0],
                                self._params["mesh_off"][1],
                                self._params["mesh_off"][2])

        def validate_params(self):
            default_params = self.default_params()
            if not (self._params["prefactor"] > 0.0):
                raise ValueError("prefactor should be a positive float")

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

            if self._params["epsilon"] == "metallic":
                self._params["epsilon"] = 0.0

            check_type_or_throw_except(
                self._params["epsilon"], 1, float,
                "epsilon should be a double or 'metallic'")

            if self._params["mesh_off"] != default_params["mesh_off"]:
                check_type_or_throw_except(self._params["mesh_off"], 3, float,
                                           "mesh_off should be a (3,) array_like of values between 0.0 and 1.0")

            if not is_valid_type(self._params["timings"], int):
                raise TypeError("DipolarP3M timings has to be an integer")
            if self._params["timings"] <= 0:
                raise ValueError("DipolarP3M timings must be > 0")

    cdef class P3M(_P3MBase):
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
        timings : :obj:`int`
            Number of force calculations during tuning.
        verbose : :obj:`bool`, optional
            If ``False``, disable log output during tuning.
        check_neutrality : :obj:`bool`, optional
            Raise a warning if the system is not electrically neutral when
            set to ``True`` (default).

        """

        def _activate_method(self):
            check_neutrality(self._params)
            if self._params["tune"]:
                self._tune()
            self._set_params_in_es_core()
            handle_errors("P3M: initialization failed")

    IF CUDA:
        cdef class P3MGPU(_P3MBase):
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
            timings : :obj:`int`
                Number of force calculations during tuning.
            verbose : :obj:`bool`, optional
                If ``False``, disable log output during tuning.
            check_neutrality : :obj:`bool`, optional
                Raise a warning if the system is not electrically neutral when
                set to ``True`` (default).

            """

            def _activate_method(self):
                cdef int mesh[3]
                self._check_and_copy_mesh_size(mesh, self._params["mesh"])

                check_neutrality(self._params)
                p3m_gpu_init(self._params["cao"], mesh, self._params["alpha"])
                handle_errors("P3M: tuning failed")
                coulomb.method = COULOMB_P3M_GPU
                if self._params["tune"]:
                    self._tune()
                p3m_gpu_init(self._params["cao"], mesh, self._params["alpha"])
                handle_errors("P3M: tuning failed")
                self._set_params_in_es_core()

            def _set_params_in_es_core(self):
                super()._set_params_in_es_core()
                handle_errors("P3M: initialization failed")

    cdef class ELC(ElectrostaticInteraction):
        """
        Electrostatics solver for systems with two periodic dimensions.
        See :ref:`Electrostatic Layer Correction (ELC)` for more details.

        Parameters
        ----------
        p3m_actor : :obj:`P3M`, required
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
            and the simulation box.
        delta_mid_bottom : :obj:`float`, optional
            Dielectric contrast :math:`\\Delta_b` between the lower boundary
            and the simulation box.
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

        def validate_params(self):
            # P3M
            if CUDA:
                if isinstance(self._params["p3m_actor"], P3MGPU):
                    raise ValueError(
                        "ELC is not set up to work with the GPU P3M")
            check_type_or_throw_except(
                self._params["p3m_actor"], 1, getattr(
                    sys.modules[__name__], "P3M"),
                "p3m_actor has to be a P3M solver")
            self._params["p3m_actor"]._params["epsilon"] = 0.0
            self._params["p3m_actor"].validate_params()
            # ELC
            check_type_or_throw_except(
                self._params["maxPWerror"], 1, float,
                "maxPWerror has to be a float")
            check_type_or_throw_except(self._params["gap_size"], 1, float,
                                       "gap_size has to be a float")
            check_type_or_throw_except(self._params["far_cut"], 1, float,
                                       "far_cut has to be a float")
            check_type_or_throw_except(
                self._params["neutralize"], 1, type(True),
                "neutralize has to be a bool")

        def valid_keys(self):
            return ["p3m_actor", "maxPWerror", "gap_size", "far_cut",
                    "neutralize", "delta_mid_top", "delta_mid_bot",
                    "const_pot", "pot_diff", "check_neutrality"]

        def required_keys(self):
            return ["p3m_actor", "maxPWerror", "gap_size"]

        def default_params(self):
            return {"maxPWerror": -1,
                    "gap_size": -1,
                    "far_cut": -1,
                    "delta_mid_top": 0,
                    "delta_mid_bot": 0,
                    "const_pot": False,
                    "pot_diff": 0.0,
                    "neutralize": True,
                    "check_neutrality": True}

        def _get_params_from_es_core(self):
            params = {}
            params.update(elc_params)
            params["p3m_actor"] = self._params["p3m_actor"]
            return params

        def _set_params_in_es_core(self):
            self._params["p3m_actor"]._set_params_in_es_core()
            if coulomb.method == COULOMB_P3M_GPU:
                raise Exception("ELC is not set up to work with the GPU P3M")

            if self._params["const_pot"]:
                self._params["delta_mid_top"] = -1
                self._params["delta_mid_bot"] = -1

            ELC_set_params(
                self._params["maxPWerror"],
                self._params["gap_size"],
                self._params["far_cut"],
                self._params["neutralize"],
                self._params["delta_mid_top"],
                self._params["delta_mid_bot"],
                self._params["const_pot"],
                self._params["pot_diff"])

        def tune(self, **tune_params_subset):
            self._params["p3m_actor"].tune(**tune_params_subset)

        def _activate_method(self):
            self._params["p3m_actor"]._activate_method()
            check_neutrality(self._params)
            self._set_params_in_es_core()

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
        timings : :obj:`int`
            Number of force calculations during tuning.

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
            if not is_valid_type(self._params["timings"], int):
                raise TypeError("DipolarP3M timings has to be an integer")
            if self._params["timings"] <= 0:
                raise ValueError("DipolarP3M timings must be > 0")

        def default_params(self):
            return {"prefactor": -1,
                    "maxPWerror": -1,
                    "far_switch_radius": -1,
                    "bessel_cutoff": -1,
                    "tune": True,
                    "timings": 1000,
                    "check_neutrality": True,
                    "verbose": True}

        def valid_keys(self):
            return ["prefactor", "maxPWerror", "far_switch_radius",
                    "bessel_cutoff", "tune", "check_neutrality", "timings",
                    "verbose"]

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
            params["timings"] = self._params["timings"]
            return params

        def _set_params_in_es_core(self):
            set_prefactor(self._params["prefactor"])
            MMM1D_set_params(
                self._params["far_switch_radius"], self._params["maxPWerror"])

        def _tune(self):
            resp = MMM1D_init()
            if resp:
                handle_errors("MMM1D: initialization failed")
            resp = mmm1d_tune(self._params["timings"], self._params["verbose"])
            if resp:
                handle_errors("MMM1D: tuning failed")
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
        cdef shared_ptr[Mmm1dgpuForce] sip
        cdef EspressoSystemInterface * interface

        def __cinit__(self):
            self.interface = &EspressoSystemInterface.Instance()
            self.sip = make_shared[Mmm1dgpuForce](dereference(self.interface))

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

            dereference(self.sip).set_params(
                < float > box_geo.length()[2], < float > coulomb.prefactor,
                self._params["maxPWerror"], self._params["far_switch_radius"],
                self._params["bessel_cutoff"])

        def _tune(self):
            dereference(self.sip).setup(dereference(self.interface))
            dereference(self.sip).tune(
                dereference(self.interface), self._params["maxPWerror"],
                self._params["far_switch_radius"], self._params["bessel_cutoff"])

        def _activate_method(self):
            check_neutrality(self._params)
            self._set_params_in_es_core()
            dereference(self.sip).activate()
            coulomb.method = COULOMB_MMM1D_GPU
            if self._params["tune"]:
                self._tune()
            self._set_params_in_es_core()

        def _deactivate_method(self):
            dereference(self.sip).deactivate()

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
                scafacos.free_handle(self.dipolar)
                mpi_bcast_coulomb_params()
