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
IF SCAFACOS == 1:
    from .scafacos import ScafacosConnector
    from . cimport scafacos

from .utils import handle_errors
from .utils import is_valid_type, check_type_or_throw_except

IF DIPOLES == 1:
    cdef class MagnetostaticInteraction(Actor):
        """Provide magnetostatic interactions.

        Parameters
        ----------
        prefactor : :obj:`float`
            Magnetostatics prefactor (:math:`\\mu_0/(4\\pi)`)

        """

        def validate_params(self):
            """Check validity of given parameters.
            """
            if not self._params["prefactor"] >= 0:
                raise ValueError("prefactor should be a positive float")

        def get_magnetostatics_prefactor(self):
            """
            Get the magnetostatics prefactor

            """
            return get_Dprefactor()

        def set_magnetostatics_prefactor(self):
            """
            Set the magnetostatics prefactor

            """
            set_Dprefactor(self._params["prefactor"])
            # also necessary on 1 CPU or GPU, does more than just broadcasting
            mpi_bcast_coulomb_params()

        def _deactivate_method(self):
            set_Dprefactor(0.0)
            disable_method_local()
            mpi_bcast_coulomb_params()

IF DP3M == 1:
    cdef class DipolarP3M(MagnetostaticInteraction):
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
            Charge-assignment order, an integer between -1 and 7.
        mesh : :obj:`int` or (3,) array_like of :obj:`int`
            The number of mesh points in x, y and z direction. Use a single
            value for cubic boxes.
        mesh_off : (3,) array_like of :obj:`float`
            Mesh offset.
        r_cut : :obj:`float`
            Real space cutoff.
        tune : :obj:`bool`, optional
            Activate/deactivate the tuning method on activation
            (default is ``True``, i.e., activated).
        timings : :obj:`int`
            Number of force calculations during tuning.

        """

        def validate_params(self):
            """Check validity of parameters.

            """
            super().validate_params()
            default_params = self.default_params()

            if is_valid_type(self._params["mesh"], int):
                pass
            else:
                check_type_or_throw_except(self._params["mesh"], 3, int,
                                           "DipolarP3M mesh has to be an integer or integer list of length 3")
                if (self._params["mesh"][0] != self._params["mesh"][1]) or \
                   (self._params["mesh"][0] != self._params["mesh"][2]):
                    raise ValueError(
                        "DipolarP3M requires a cubic box")

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

        def valid_keys(self):
            return ["prefactor", "alpha_L", "r_cut_iL", "mesh", "mesh_off",
                    "cao", "accuracy", "epsilon", "cao_cut", "a", "ai",
                    "alpha", "r_cut", "cao3", "tune", "timings", "verbose"]

        def required_keys(self):
            return ["accuracy", ]

        def default_params(self):
            return {"cao": -1,
                    "r_cut": -1,
                    "accuracy": -1,
                    "mesh": -1,
                    "epsilon": 0.0,
                    "mesh_off": [-1, -1, -1],
                    "tune": True,
                    "timings": 10,
                    "verbose": True}

        def _get_params_from_es_core(self):
            params = {}
            params.update(dp3m.params)
            params["prefactor"] = self.get_magnetostatics_prefactor()
            params["tune"] = self._params["tune"]
            params["timings"] = self._params["timings"]
            return params

        def _set_params_in_es_core(self):
            if hasattr(self._params["mesh"], "__getitem__"):
                mesh = self._params["mesh"][0]
            else:
                mesh = self._params["mesh"]

            self.set_magnetostatics_prefactor()
            dp3m_set_eps(self._params["epsilon"])
            dp3m_set_mesh_offset(self._params["mesh_off"][0],
                                 self._params["mesh_off"][1],
                                 self._params["mesh_off"][2])
            dp3m_set_params(self._params["r_cut"], mesh, self._params["cao"],
                            self._params["alpha"], self._params["accuracy"])

        def _tune(self):
            if hasattr(self._params["mesh"], "__getitem__"):
                mesh = self._params["mesh"][0]
            else:
                mesh = self._params["mesh"]

            self.set_magnetostatics_prefactor()
            dp3m_set_eps(self._params["epsilon"])
            dp3m_set_tune_params(self._params["r_cut"], mesh,
                                 self._params["cao"], self._params["accuracy"])
            tuning_error = dp3m_adaptive_tune(
                self._params["timings"], self._params["verbose"])
            if tuning_error:
                handle_errors("DipolarP3M: tuning failed")
            self._params.update(self._get_params_from_es_core())

        def _activate_method(self):
            if self._params["tune"]:
                self._tune()

            self._set_params_in_es_core()
            mpi_bcast_coulomb_params()

        def _deactivate_method(self):
            dp3m_deactivate()
            super()._deactivate_method()

IF DIPOLES == 1:
    cdef class DipolarDirectSumCpu(MagnetostaticInteraction):
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

        def default_params(self):
            return {}

        def required_keys(self):
            return ()

        def valid_keys(self):
            return ("prefactor",)

        def _get_params_from_es_core(self):
            return {"prefactor": self.get_magnetostatics_prefactor()}

        def _activate_method(self):
            self._set_params_in_es_core()
            mpi_bcast_coulomb_params()

        def _set_params_in_es_core(self):
            self.set_magnetostatics_prefactor()
            dawaanr_set_params()

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

        def default_params(self):
            return {}

        def required_keys(self):
            return ("n_replica",)

        def valid_keys(self):
            return ("prefactor", "n_replica")

        def _get_params_from_es_core(self):
            return {"prefactor": self.get_magnetostatics_prefactor(),
                    "n_replica": mdds_get_n_replica()}

        def _activate_method(self):
            self._set_params_in_es_core()
            mpi_bcast_coulomb_params()

        def _set_params_in_es_core(self):
            self.set_magnetostatics_prefactor()
            mdds_set_params(self._params["n_replica"])

    IF SCAFACOS_DIPOLES == 1:
        class Scafacos(ScafacosConnector, MagnetostaticInteraction):

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

            """

            dipolar = True

            # Explicit constructor needed due to multiple inheritance
            def __init__(self, *args, **kwargs):
                Actor.__init__(self, *args, **kwargs)

            def _activate_method(self):
                self.set_magnetostatics_prefactor()
                self._set_params_in_es_core()

            def _deactivate_method(self):
                scafacos.free_handle(self.dipolar)
                super()._deactivate_method()

            def default_params(self):
                return {}

    IF(CUDA == 1) and (DIPOLES == 1) and (ROTATION == 1):
        cdef class DipolarDirectSumGpu(MagnetostaticInteraction):
            """
            Calculate magnetostatic interactions by direct summation over all
            pairs. See :ref:`Dipolar direct sum` for more details.

            If the system has periodic boundaries, the minimum image convention
            is applied in the respective directions.

            This is the GPU version of :class:`espressomd.magnetostatics.DipolarDirectSumCpu`
            but uses floating point precision.

            Requires features ``DIPOLAR_DIRECT_SUM`` and ``CUDA``.

            Parameters
            ----------
            prefactor : :obj:`float`
                Magnetostatics prefactor (:math:`\\mu_0/(4\\pi)`)

            """
            cdef shared_ptr[DipolarDirectSum] sip

            def __cinit__(self):
                self.sip = make_shared[DipolarDirectSum](
                    EspressoSystemInterface.Instance())

            def default_params(self):
                return {}

            def required_keys(self):
                return ()

            def valid_keys(self):
                return ("prefactor",)

            def _get_params_from_es_core(self):
                return {"prefactor": self.get_magnetostatics_prefactor()}

            def _activate_method(self):
                self._set_params_in_es_core()
                dereference(self.sip).activate()

            def _deactivate_method(self):
                super()._deactivate_method()
                dereference(self.sip).deactivate()

            def _set_params_in_es_core(self):
                self.set_magnetostatics_prefactor()
                dereference(self.sip).set_params()

    IF(DIPOLAR_BARNES_HUT == 1):
        cdef class DipolarBarnesHutGpu(MagnetostaticInteraction):

            """
            Calculates magnetostatic interactions by direct summation over all
            pairs. See :ref:`Barnes-Hut octree sum on GPU` for more details.

            TODO: If the system has periodic boundaries, the minimum image
            convention is applied.
            """
            cdef shared_ptr[DipolarBarnesHut] sip

            def __cinit__(self):
                self.sip = make_shared[DipolarBarnesHut](
                    EspressoSystemInterface.Instance())

            def default_params(self):
                return {"epssq": 100.0,
                        "itolsq": 4.0}

            def required_keys(self):
                return ()

            def valid_keys(self):
                return ("prefactor", "epssq", "itolsq")

            def _get_params_from_es_core(self):
                return {"prefactor": self.get_magnetostatics_prefactor()}

            def _activate_method(self):
                self._set_params_in_es_core()
                dereference(self.sip).activate()

            def _deactivate_method(self):
                super()._deactivate_method()
                dereference(self.sip).deactivate()

            def _set_params_in_es_core(self):
                self.set_magnetostatics_prefactor()
                dereference(self.sip).set_params(
                    self._params["epssq"], self._params["itolsq"])
