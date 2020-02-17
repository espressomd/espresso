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
from .utils cimport handle_errors, check_type_or_throw_except, check_range_or_except

IF ELECTROSTATICS and P3M:
    from espressomd.electrostatics import check_neutrality

    cdef class ElectrostaticExtensions(actors.Actor):
        pass

    cdef class ELC(ElectrostaticExtensions):
        """
        Electrostatics solver for systems with two periodic dimensions.
        See :ref:`Electrostatic Layer Correction (ELC)` for more details.

        Parameters
        ----------
        gap_size : :obj:`float`, required
            The gap size gives the height :math:`h` of the empty region between
            the system box and the neighboring artificial images. |es| does not
            make sure that the gap is actually empty, this is the user's
            responsibility. The method will run even if the condition is not
            fulfilled, however, the error bound will not be reached. Therefore
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
            default_params = self.default_params()
            check_type_or_throw_except(
                self._params["maxPWerror"], 1, float, "")
            check_range_or_except(
                self._params, "maxPWerror", 0, False, "inf", True)
            check_type_or_throw_except(self._params["gap_size"], 1, float, "")
            check_range_or_except(
                self._params, "gap_size", 0, False, "inf", True)
            check_type_or_throw_except(self._params["far_cut"], 1, float, "")
            check_type_or_throw_except(
                self._params["neutralize"], 1, type(True), "")

        def valid_keys(self):
            return ["maxPWerror", "gap_size", "far_cut", "neutralize",
                    "delta_mid_top", "delta_mid_bot", "const_pot", "pot_diff",
                    "check_neutrality"]

        def required_keys(self):
            return ["maxPWerror", "gap_size"]

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
            return params

        def _set_params_in_es_core(self):
            if coulomb.method == COULOMB_P3M_GPU:
                raise Exception(
                    "ELC tuning failed, ELC is not set up to work with the GPU P3M")

            if self._params["const_pot"]:
                self._params["delta_mid_top"] = -1
                self._params["delta_mid_bot"] = -1

            if ELC_set_params(
                self._params["maxPWerror"],
                self._params["gap_size"],
                self._params["far_cut"],
                self._params["neutralize"],
                self._params["delta_mid_top"],
                self._params["delta_mid_bot"],
                self._params["const_pot"],
                    self._params["pot_diff"]):
                handle_errors(
                    "ELC tuning failed, ELC is not set up to work with the GPU P3M")

        def _activate_method(self):
            check_neutrality(self._params)
            self._set_params_in_es_core()

        def _deactivate_method(self):
            raise Exception(
                "Unable to remove ELC as the state of the underlying electrostatics method will remain unclear.")

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
        epsilons : (``n_icc``, ) array_like :obj:`float`, optional
            Dielectric constant associated to the areas.

        """

        def validate_params(self):
            default_params = self.default_params()

            check_type_or_throw_except(self._params["n_icc"], 1, int, "")
            check_range_or_except(
                self._params, "n_icc", 1, True, "inf", True)

            check_type_or_throw_except(
                self._params["convergence"], 1, float, "")
            check_range_or_except(
                self._params, "convergence", 0, False, "inf", True)

            check_type_or_throw_except(
                self._params["relaxation"], 1, float, "")
            check_range_or_except(
                self._params, "relaxation", 0, False, "inf", True)

            check_type_or_throw_except(
                self._params["ext_field"], 3, float, "")

            check_type_or_throw_except(
                self._params["max_iterations"], 1, int, "")
            check_range_or_except(
                self._params, "max_iterations", 0, False, "inf", True)

            check_type_or_throw_except(
                self._params["first_id"], 1, int, "")
            check_range_or_except(
                self._params, "first_id", 0, True, "inf", True)

            check_type_or_throw_except(
                self._params["eps_out"], 1, float, "")

            n_icc = self._params["n_icc"]

            # Required list input
            self._params["normals"] = np.array(self._params["normals"])
            if self._params["normals"].size != n_icc * 3:
                raise ValueError(
                    "Expecting normal list with " + str(n_icc * 3) + " entries.")
            check_type_or_throw_except(self._params["normals"], n_icc,
                                       np.ndarray, "Error in normal list.")

            check_type_or_throw_except(
                self._params["areas"], n_icc, float, "Error in area list.")

            # Not Required
            if "sigmas" in self._params.keys():
                check_type_or_throw_except(
                    self._params["sigmas"], n_icc, float, "Error in sigma list.")
            else:
                self._params["sigmas"] = np.zeros(n_icc)

            if "epsilons" in self._params.keys():
                check_type_or_throw_except(
                    self._params["epsilons"], n_icc, float, "Error in epsilon list.")
            else:
                self._params["epsilons"] = np.zeros(n_icc)

        def valid_keys(self):
            return ["n_icc", "convergence", "relaxation", "ext_field",
                    "max_iterations", "first_id", "eps_out", "normals",
                    "areas", "sigmas", "epsilons", "check_neutrality"]

        def required_keys(self):
            return ["n_icc", "normals", "areas"]

        def default_params(self):
            return {"n_icc": 0,
                    "convergence": 1e-3,
                    "relaxation": 0.7,
                    "ext_field": [0, 0, 0],
                    "max_iterations": 100,
                    "first_id": 0,
                    "esp_out": 1,
                    "normals": [],
                    "areas": [],
                    "sigmas": [],
                    "epsilons": [],
                    "check_neutrality": True}

        def _get_params_from_es_core(self):
            params = {}
            params["n_icc"] = iccp3m_cfg.n_ic

            # Fill Lists
            normals = []
            areas = []
            sigmas = []
            epsilons = []
            for i in range(iccp3m_cfg.n_ic):
                normals.append([iccp3m_cfg.normals[i][0], iccp3m_cfg.normals[
                               i][1], iccp3m_cfg.normals[i][2]])
                areas.append(iccp3m_cfg.areas[i])
                epsilons.append(iccp3m_cfg.ein[i])
                sigmas.append(iccp3m_cfg.sigma[i])

            params["normals"] = normals
            params["areas"] = areas
            params["epsilons"] = epsilons
            params["sigmas"] = sigmas

            params["ext_field"] = [iccp3m_cfg.ext_field[0],
                                   iccp3m_cfg.ext_field[1], iccp3m_cfg.ext_field[2]]
            params["first_id"] = iccp3m_cfg.first_id
            params["max_iterations"] = iccp3m_cfg.num_iteration
            params["convergence"] = iccp3m_cfg.convergence
            params["relaxation"] = iccp3m_cfg.relax
            params["eps_out"] = iccp3m_cfg.eout

            return params

        def _set_params_in_es_core(self):
            # First set number of icc particles
            iccp3m_cfg.n_ic = self._params["n_icc"]
            # Allocate ICC lists
            iccp3m_alloc_lists()

            # Fill Lists
            for i in range(iccp3m_cfg.n_ic):
                iccp3m_cfg.normals[i][0] = self._params["normals"][i][0]
                iccp3m_cfg.normals[i][1] = self._params["normals"][i][1]
                iccp3m_cfg.normals[i][2] = self._params["normals"][i][2]

                iccp3m_cfg.areas[i] = self._params["areas"][i]
                iccp3m_cfg.ein[i] = self._params["epsilons"][i]
                iccp3m_cfg.sigma[i] = self._params["sigmas"][i]

            iccp3m_cfg.ext_field[0] = self._params["ext_field"][0]
            iccp3m_cfg.ext_field[1] = self._params["ext_field"][1]
            iccp3m_cfg.ext_field[2] = self._params["ext_field"][2]
            iccp3m_cfg.first_id = self._params["first_id"]
            iccp3m_cfg.num_iteration = self._params["max_iterations"]
            iccp3m_cfg.convergence = self._params["convergence"]
            iccp3m_cfg.relax = self._params["relaxation"]
            iccp3m_cfg.eout = self._params["eps_out"]

            # Broadcasts vars
            mpi_iccp3m_init()

        def _activate_method(self):
            check_neutrality(self._params)
            self._set_params_in_es_core()

        def _deactivate_method(self):
            iccp3m_cfg.n_ic = 0
            # Allocate ICC lists
            iccp3m_alloc_lists()

            # Broadcasts vars
            mpi_iccp3m_init()

        def last_iterations(self):
            """
            Number of iterations needed in last relaxation to
            reach the convergence criterion.

            Returns
            -------
            iterations : :obj:`int`
                Number of iterations

            """
            return iccp3m_cfg.citeration
