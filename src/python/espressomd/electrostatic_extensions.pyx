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

from . cimport utils
include "myconfig.pxi"
from . cimport actors
from . import actors
import numpy as np
from .utils import handle_errors, array_locked
from .utils cimport check_type_or_throw_except, check_range_or_except, Vector3d, make_Vector3d, make_array_locked, make_array_locked_vector
from libcpp.vector cimport vector

IF ELECTROSTATICS and P3M:
    from espressomd.electrostatics import check_neutrality

    cdef class ElectrostaticExtensions(actors.Actor):
        pass

    cdef class ICC(ElectrostaticExtensions):
        """
        Interface to the induced charge calculation scheme for dielectric
        interfaces. See :ref:`Dielectric interfaces with the ICC algorithm`
        for more details.

        Parameters
        ----------
        n_icc : :obj:`int`
            Total number of ICC Particles.
        first_id : :obj:`int`, optional
            ID of the first ICC Particle.
        convergence : :obj:`float`, optional
            Abort criteria of the iteration. It corresponds to the maximum relative
            change of any of the interface particle's charge.
        relaxation : :obj:`float`, optional
            SOR relaxation parameter.
        ext_field : :obj:`float`, optional
            Homogeneous electric field added to the calculation of dielectric boundary forces.
        max_iterations : :obj:`int`, optional
            Maximal number of iterations.
        eps_out : :obj:`float`, optional
            Relative permittivity of the outer region (where the particles are).
        normals : (``n_icc``, 3) array_like :obj:`float`
            Normal vectors pointing into the outer region.
        areas : (``n_icc``, ) array_like :obj:`float`
            Areas of the discretized surface.
        sigmas : (``n_icc``, ) array_like :obj:`float`, optional
            Additional surface charge density in the absence of any charge
            induction.
        epsilons : (``n_icc``, ) array_like :obj:`float`
            Dielectric constant associated to the areas.

        """

        def validate_params(self):
            check_type_or_throw_except(self._params["n_icc"], 1, int, "")

            check_type_or_throw_except(
                self._params["first_id"], 1, int, "")

            check_type_or_throw_except(
                self._params["convergence"], 1, float, "")

            check_type_or_throw_except(
                self._params["relaxation"], 1, float, "")

            check_type_or_throw_except(
                self._params["ext_field"], 3, float, "")

            check_type_or_throw_except(
                self._params["max_iterations"], 1, int, "")

            check_type_or_throw_except(
                self._params["eps_out"], 1, float, "")

            n_icc = self._params["n_icc"]
            assert n_icc >= 0, "ICC: invalid number of particles"

            self._params["normals"] = np.array(self._params["normals"])
            if self._params["normals"].size != n_icc * 3:
                raise ValueError(
                    "Expecting normal list with " + str(n_icc * 3) + " entries.")
            check_type_or_throw_except(self._params["normals"], n_icc,
                                       np.ndarray, "Error in normal list.")

            check_type_or_throw_except(
                self._params["areas"], n_icc, float, "Error in area list.")

            if "sigmas" in self._params.keys():
                check_type_or_throw_except(
                    self._params["sigmas"], n_icc, float, "Error in sigma list.")
            else:
                self._params["sigmas"] = np.zeros(n_icc)

            check_type_or_throw_except(
                self._params["epsilons"], n_icc, float, "Error in epsilon list.")

        def valid_keys(self):
            return ["n_icc", "convergence", "relaxation", "ext_field",
                    "max_iterations", "first_id", "eps_out", "normals",
                    "areas", "sigmas", "epsilons", "check_neutrality"]

        def required_keys(self):
            return ["n_icc", "normals", "areas", "epsilons"]

        def default_params(self):
            return {"convergence": 1e-3,
                    "relaxation": 0.7,
                    "ext_field": [0, 0, 0],
                    "max_iterations": 100,
                    "first_id": 0,
                    "eps_out": 1,
                    "check_neutrality": True}

        def _get_params_from_es_core(self):
            params = {}
            params["n_icc"] = icc_cfg.n_icc
            params["first_id"] = icc_cfg.first_id
            params["max_iterations"] = icc_cfg.num_iteration
            params["convergence"] = icc_cfg.convergence
            params["relaxation"] = icc_cfg.relax
            params["eps_out"] = icc_cfg.eout
            params["normals"] = make_array_locked_vector(icc_cfg.normals)
            params["areas"] = array_locked(icc_cfg.areas)
            params["epsilons"] = array_locked(icc_cfg.ein)
            params["sigmas"] = array_locked(icc_cfg.sigma)
            params["ext_field"] = make_array_locked(icc_cfg.ext_field)

            return params

        def _set_params_in_es_core(self):
            cdef Vector3d ext_field = make_Vector3d(self._params["ext_field"])
            cdef vector[double] areas, e_in, sigma
            cdef vector[Vector3d] normals
            areas.resize(self._params["n_icc"])
            e_in.resize(self._params["n_icc"])
            sigma.resize(self._params["n_icc"])
            normals.resize(self._params["n_icc"])

            for i in range(self._params["n_icc"]):
                areas[i] = self._params["areas"][i]
                e_in[i] = self._params["epsilons"][i]
                sigma[i] = self._params["sigmas"][i]

                for j in range(3):
                    normals[i][j] = self._params["normals"][i][j]

            icc_set_params(self._params["n_icc"],
                           self._params["convergence"],
                           self._params["relaxation"],
                           ext_field,
                           self._params["max_iterations"],
                           self._params["first_id"],
                           self._params["eps_out"],
                           areas,
                           e_in,
                           sigma,
                           normals)

        def _activate_method(self):
            check_neutrality(self._params)
            self._set_params_in_es_core()

        def _deactivate_method(self):
            icc_deactivate()

        def last_iterations(self):
            """
            Number of iterations needed in last relaxation to
            reach the convergence criterion.

            Returns
            -------
            iterations : :obj:`int`
                Number of iterations

            """
            return icc_cfg.citeration
