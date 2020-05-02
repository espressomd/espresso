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
include "myconfig.pxi"
from .utils import requires_experimental_features
from .actors cimport Actor
IF SCAFACOS == 1:
    from .scafacos import ScafacosConnector
    from . cimport scafacos

from .utils cimport handle_errors
from .utils import is_valid_type, check_type_or_throw_except, to_str

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

        def set_magnetostatics_prefactor(self):
            """
            Set the magnetostatics prefactor

            """
            set_Dprefactor(self._params["prefactor"])
            handle_errors("Could not set magnetostatic prefactor")
            # also necessary on 1 CPU or GPU, does more than just broadcasting
            mpi_bcast_coulomb_params()

        def _get_active_method_from_es_core(self):
            return dipole.method

        def _deactivate_method(self):
            set_Dprefactor(0.0)
            dipole.method = DIPOLAR_NONE
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

        """

        def validate_params(self):
            """Check validity of parameters.

            """
            super().validate_params()
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

            if not (self._params["cao"] >= -1 and self._params["cao"] <= 7):
                raise ValueError(
                    "P3M cao has to be an integer between -1 and 7")

            if not (self._params["accuracy"] > 0):
                raise ValueError("P3M accuracy has to be positive")

            if self._params["epsilon"] == "metallic":
                self._params["epsilon"] = 0.0

            if not (is_valid_type(self._params["epsilon"], float)
                    or self._params["epsilon"] == "metallic"):
                raise ValueError("epsilon should be a double or 'metallic'")

            if not (self._params["inter"] == default_params["inter"]
                    or self._params["inter"] > 0):
                raise ValueError("inter should be a positive integer")

            if self._params["mesh_off"] != default_params["mesh_off"]:
                check_type_or_throw_except(self._params["mesh_off"], 3, float,
                                           "mesh_off should be a (3,) array_like of values between 0.0 and 1.0")

        def valid_keys(self):
            return ["prefactor", "alpha_L", "r_cut_iL", "mesh", "mesh_off",
                    "cao", "inter", "accuracy", "epsilon", "cao_cut", "a", "ai",
                    "alpha", "r_cut", "inter2", "cao3", "additional_mesh", "tune"]

        def required_keys(self):
            return ["accuracy", ]

        def default_params(self):
            return {"cao": -1,
                    "inter": -1,
                    "r_cut": -1,
                    "accuracy": -1,
                    "mesh": -1,
                    "epsilon": 0.0,
                    "mesh_off": [-1, -1, -1],
                    "tune": True}

        def _get_params_from_es_core(self):
            params = {}
            params.update(dp3m.params)
            params["prefactor"] = dipole.prefactor
            params["tune"] = self._params["tune"]
            return params

        def _set_params_in_es_core(self):
            self.set_magnetostatics_prefactor()
            dp3m_set_eps(self._params["epsilon"])
            self.python_dp3m_set_mesh_offset(self._params["mesh_off"])
            self.python_dp3m_set_params(
                self._params["r_cut"], self._params["mesh"],
                self._params["cao"], self._params["alpha"], self._params["accuracy"])

        def _tune(self):
            self.set_magnetostatics_prefactor()
            dp3m_set_eps(self._params["epsilon"])
            self.python_dp3m_set_tune_params(
                self._params["r_cut"], self._params["mesh"],
                self._params["cao"], -1., self._params["accuracy"], self._params["inter"])
            resp, log = self.python_dp3m_adaptive_tune()
            if resp:
                raise Exception(
                    "failed to tune dipolar P3M parameters to required accuracy")
            print(to_str(log))
            self._params.update(self._get_params_from_es_core())

        def _activate_method(self):
            if self._params["tune"]:
                self._tune()

            self._set_params_in_es_core()
            mpi_bcast_coulomb_params()

        def _deactivate_method(self):
            dp3m_deactivate()
            super()._deactivate_method()

        def python_dp3m_set_mesh_offset(self, mesh_off):
            cdef double mesh_offset[3]
            mesh_offset[0] = mesh_off[0]
            mesh_offset[1] = mesh_off[1]
            mesh_offset[2] = mesh_off[2]
            return dp3m_set_mesh_offset(
                mesh_offset[0], mesh_offset[1], mesh_offset[2])

        def python_dp3m_adaptive_tune(self):
            cdef char * log = NULL
            cdef int response
            response = dp3m_adaptive_tune(& log)
            handle_errors(
                "dipolar P3M_init: k-space cutoff is larger than half of box dimension")
            return response, log

        def python_dp3m_set_params(self, p_r_cut, p_mesh, p_cao, p_alpha,
                                   p_accuracy):
            cdef int mesh
            cdef double r_cut
            cdef int cao
            cdef double alpha
            cdef double accuracy
            r_cut = p_r_cut
            cao = p_cao
            alpha = p_alpha
            accuracy = p_accuracy
            if hasattr(p_mesh, "__getitem__"):
                mesh = p_mesh[0]
            else:
                mesh = p_mesh
            dp3m_set_params(r_cut, mesh, cao, alpha, accuracy)

        def python_dp3m_set_tune_params(self, p_r_cut, p_mesh, p_cao, p_alpha,
                                        p_accuracy, p_n_interpol):
            cdef int mesh
            cdef double r_cut
            cdef int cao
            cdef double alpha
            cdef double accuracy
            cdef int n_interpol
            r_cut = p_r_cut
            cao = p_cao
            alpha = p_alpha
            accuracy = p_accuracy
            n_interpol = p_n_interpol
            mesh = p_mesh
            dp3m_set_tune_params(r_cut, mesh, cao, alpha, accuracy, n_interpol)

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
            return {"prefactor": dipole.prefactor}

        def _activate_method(self):
            self._set_params_in_es_core()
            mpi_bcast_coulomb_params()

        def _set_params_in_es_core(self):
            self.set_magnetostatics_prefactor()
            dawaanr_set_params()
            handle_errors("Could not activate magnetostatics method "
                          + self.__class__.__name__)

    @requires_experimental_features("No test coverage")
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
            return {"prefactor": dipole.prefactor,
                    "n_replica": Ncut_off_magnetic_dipolar_direct_sum}

        def _activate_method(self):
            self._set_params_in_es_core()
            mpi_bcast_coulomb_params()

        def _set_params_in_es_core(self):
            self.set_magnetostatics_prefactor()
            mdds_set_params(self._params["n_replica"])
            handle_errors("Could not activate magnetostatics method "
                          + self.__class__.__name__)

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
                set_Dprefactor(self._params["prefactor"])
                self._set_params_in_es_core()

            def _deactivate_method(self):
                dipole.method = DIPOLAR_NONE
                scafacos.free_handle()
                mpi_bcast_coulomb_params()

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
                return {"prefactor": dipole.prefactor}

            def _activate_method(self):
                self._set_params_in_es_core()

            def _deactivate_method(self):
                super()._deactivate_method()
                deactivate_dipolar_direct_sum_gpu()

            def _set_params_in_es_core(self):
                self.set_magnetostatics_prefactor()
                activate_dipolar_direct_sum_gpu()

    IF(DIPOLAR_BARNES_HUT == 1):
        cdef class DipolarBarnesHutGpu(MagnetostaticInteraction):

            """
            Calculates magnetostatic interactions by direct summation over all
            pairs. See :ref:`Barnes-Hut octree sum on GPU` for more details.

            TODO: If the system has periodic boundaries, the minimum image
            convention is applied.
            """

            def default_params(self):
                return {"epssq": 100.0,
                        "itolsq": 4.0}

            def required_keys(self):
                return ()

            def valid_keys(self):
                return ("prefactor", "epssq", "itolsq")

            def _get_params_from_es_core(self):
                return {"prefactor": dipole.prefactor}

            def _activate_method(self):
                self._set_params_in_es_core()

            def _deactivate_method(self):
                super()._deactivate_method()
                deactivate_dipolar_barnes_hut()

            def _set_params_in_es_core(self):
                self.set_magnetostatics_prefactor()
                activate_dipolar_barnes_hut(
                    self._params["epssq"], self._params["itolsq"])
